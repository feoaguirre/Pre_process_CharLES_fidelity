"""
Module: mesh_planner
Description: Calculates Voronoi mesh resolution parameters for CharLES Stitch.
             Estimates near-wall cell sizing based on user-defined y+ constraints
             and computes the required binary refinement levels.
"""

import math
from dataclasses import dataclass
from typing import Any

@dataclass
class MeshState:
    """Stores the calculated mesh parameters required by stitch.in"""
    hcp_delta: float
    max_refinement_level: int
    nlayers: int
    nsmooth: int
    refinement_height: float
    half_w: float
    target_y_spacing: float


class MeshPlanner:
    """
    Translates fluid mechanics constraints into grid discretization parameters.
    """
    def __init__(self, state: Any, physics_state: Any):
        """
        Args:
            state: SimulationState (parsed from YAML)
            physics_state: PhysicsState (from physics_engine)
        """
        self.state = state
        self.physics = physics_state

    def _estimate_wall_spacing(self) -> float:
        """
        Calculates the non-dimensional wall spacing (Delta y) corresponding 
        to the target y+. 
        Strictly enforces that the unit length in the mesh is L_ref = delta* = 1.0.
        """
        # 1. Non-dimensional reference values (Strictly mapped to L_ref = delta* = 1.0)
        mach = self.state.flow_physics.mach_number
        gamma = self.state.flow_physics.gamma
        prandtl = getattr(self.state.flow_physics, 'pr_lam', 0.72)
        re_delta_star = self.state.flow_physics.target_reynolds
        
        rho_inf = 1.0
        u_inf = 1.0
        # Since Re = (rho * U * L) / mu and L = delta* = 1.0
        mu_inf = (rho_inf * u_inf * 1.0) / re_delta_star 
        
        # 2. Wall thermodynamic properties 
        wall_bc = self.state.boundary_conditions.wall_bc
        wall_T = getattr(self.state.boundary_conditions, 'wall_T', -1.0)
        
        r = math.sqrt(prandtl)
        t_aw = 1.0 + r * ((gamma - 1.0) / 2.0) * (mach ** 2) # t_inf = 1.0
        t_w = t_aw if wall_bc.upper() == "ADIABATIC" else (wall_T if wall_T > 0 else t_aw)
        t_star = 1.0 + 0.5 * (t_w - 1.0) + 0.22 * (t_aw - 1.0)
        
        n = getattr(self.state.flow_physics, 'mu_power_law', 0.76)
        mu_w = mu_inf * (t_w ** n)
        rho_w = rho_inf * (1.0 / t_w) # Constant static pressure across BL
        
        # 3. Local Skin Friction at the leading edge of the gap (x = L2)
        c_star = t_star ** (n - 1.0)
        
        # We use the already calculated L2 length from the physics state!
        l2_length = self.physics.pregap_noslip_length
        target_re_x = re_delta_star * l2_length
        
        # Compressible Skin Friction Coefficient
        cf_inc = 0.664 / math.sqrt(target_re_x)
        cf_comp = cf_inc * math.sqrt(c_star)
        
        tau_w = 0.5 * rho_inf * (u_inf ** 2) * cf_comp
        u_tau = math.sqrt(tau_w / rho_w)
        
        # 4. Strict Non-Dimensional Delta y mapped to the mesh geometry
        target_y_plus = self.state.mesh_control.target_y_plus
        delta_y = (target_y_plus * mu_w) / (rho_w * u_tau)
        
        return delta_y

    def calculate_mesh_parameters(self) -> MeshState:
        """
        Calculates the hierarchical refinement tree for CharLES Stitch.
        Returns a MeshState object to be injected into the template writer,
        respecting user boundaries defined in the YAML configuration.
        """
        # Calculate required near-wall spacing based on delta* = 1.0
        target_delta_y = self._estimate_wall_spacing()
        
        # User Configuration Overrides
        user_hcp_delta = getattr(self.state.mesh_control, 'hcp_delta', None)
        user_max_level = getattr(self.state.mesh_control, 'max_refinement_level', None)
        
        # Determine HCP_DELTA (Base Cell Size)
        if user_hcp_delta and user_hcp_delta > 0:
            hcp_delta = user_hcp_delta
        else:
            hcp_delta = min(5.0 * self.physics.estimated_delta99, self.physics.domain_height / 10.0)
            
        # Determine MAX_REFINEMENT_LEVEL
        if target_delta_y >= hcp_delta:
            max_level = 0
        else:
            max_level = math.ceil(math.log2(hcp_delta / target_delta_y))
            
        # Apply strict caps to avoid RAM explosion
        if user_max_level is not None:
            max_level = min(max_level, user_max_level)
        else:
            max_level = min(max_level, 9) # Hard safety net
        
        smallest_cell_size = hcp_delta / (2 ** max_level)
        
        # Mandatory Quasi-2D Coupling for CharLES Stitch
        domain_z = self.state.geometry.domain_z
        if domain_z > 0.0:
            half_w = domain_z / 2.0
        else:
            half_w = smallest_cell_size / 2.0
            
        refinement_height = 1.5 * self.physics.estimated_delta99
        
        print("\n=== MESH PLANNER DEBUG ===")
        print(f"Target Y+: {self.state.mesh_control.target_y_plus}")
        print(f"Required Delta Y (delta* = 1): {target_delta_y:.6f}")
        print(f"HCP_DELTA (Far-field): {hcp_delta:.4f}")
        print(f"MAX_LEVEL applied: {max_level}")
        print(f"Actual Wall Cell Size: {smallest_cell_size:.6f}")
        print(f"Will target Y+ be reached? {'YES' if smallest_cell_size <= target_delta_y else 'NO (Limited by Max Level)'}")
        print("==========================\n")
        
        return MeshState(
            hcp_delta=round(hcp_delta, 4),
            max_refinement_level=max_level,
            nlayers=getattr(self.state.mesh_control, 'nlayers', 8),
            nsmooth=getattr(self.state.mesh_control, 'nsmooth', 20),
            refinement_height=round(refinement_height, 4),
            half_w=round(half_w, 8),
            target_y_spacing=target_delta_y
        )