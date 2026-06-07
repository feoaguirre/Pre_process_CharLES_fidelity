"""
Module: mesh_planner
Description: Computes the Voronoi volumetric mesh resolution parameters required 
             by the CharLES Stitch generator. It translates the aerodynamic requirements 
             (like near-wall y+ spacing) into explicit binary octree refinement levels 
             and safely caps memory allocation to prevent Out-Of-Memory (OOM) failures.
"""

import math
from dataclasses import dataclass
from typing import Any

@dataclass
class MeshState:
    """Data container encapsulating the structural parameters injected into stitch.in"""
    hcp_delta: float             # Base element size in the far-field (unrefined region)
    max_refinement_level: int    # Number of successive octree bisections applied near the wall
    nlayers: int                 # Number of continuous elements maintaining the highest refinement level
    nsmooth: int                 # Mesh transition smoothing factor
    refinement_height: float     # Absolute Y-coordinate up to which the max refinement is enforced
    half_w: float                # Extrusion distance applied for quasi-2D isotropic bounding
    target_y_spacing: float      # The mathematical physical wall spacing (Delta Y) required to achieve target y+


class MeshPlanner:
    """
    Acts as the translation layer between theoretical fluid mechanics constraints 
    and the discrete grid parameters utilized by the volumetric mesher.
    """
    def __init__(self, state: Any, physics_state: Any):
        """
        Args:
            state: SimulationState (parsed securely from YAML)
            physics_state: PhysicsState (computed properties from the physics_engine)
        """
        self.state = state
        self.physics = physics_state

    def _estimate_wall_spacing(self) -> float:
        """
        Calculates the absolute physical wall spacing (Delta y) required to satisfy 
        the target y+ condition at the leading edge of the structure.
        
        Strictly enforces the assumption that the non-dimensional base length 
        in the mesh coordinate system is L_ref = delta* = 1.0.
        """
        # 1. Non-dimensional reference values mapped to L_ref = delta* = 1.0
        mach = self.state.flow_physics.mach_number
        gamma = self.state.flow_physics.gamma
        prandtl = getattr(self.state.flow_physics, 'pr_lam', 0.72)
        re_delta_star = self.state.flow_physics.target_reynolds
        
        rho_inf = 1.0
        u_inf = 1.0
        
        # Dynamic viscosity evaluated at freestream conditions
        # mu_inf = (rho_inf * u_inf * L_ref) / Re
        mu_inf = (rho_inf * u_inf * 1.0) / re_delta_star 
        
        # 2. Evaluate thermodynamic properties at the wall (Adiabatic vs Isothermal)
        wall_bc = self.state.boundary_conditions.wall_bc
        wall_T = getattr(self.state.boundary_conditions, 'wall_T', -1.0)
        
        r = math.sqrt(prandtl)
        t_aw = 1.0 + r * ((gamma - 1.0) / 2.0) * (mach ** 2) # Assuming T_inf = 1.0
        
        # Assign wall temperature based on thermal boundary condition
        t_w = t_aw if wall_bc.upper() == "ADIABATIC" else (wall_T if wall_T > 0 else t_aw)
        
        # Eckert's reference temperature methodology
        t_star = 1.0 + 0.5 * (t_w - 1.0) + 0.22 * (t_aw - 1.0)
        
        # Wall fluid properties utilizing the power-law viscosity model
        n = getattr(self.state.flow_physics, 'mu_power_law', 0.76)
        mu_w = mu_inf * (t_w ** n)
        rho_w = rho_inf * (1.0 / t_w) # Derived assuming constant static pressure across the boundary layer
        
        # 3. Local Skin Friction evaluation at the structural leading edge (x = L2)
        c_star = t_star ** (n - 1.0)
        
        # Extract the absolute development length (L2) pre-calculated by the physics engine
        l2_length = self.physics.pregap_noslip_length
        target_re_x = re_delta_star * l2_length
        
        # Compressible Blasius Skin Friction Coefficient
        cf_inc = 0.664 / math.sqrt(target_re_x)
        cf_comp = cf_inc * math.sqrt(c_star)
        
        # Wall shear stress and friction velocity (u_tau)
        tau_w = 0.5 * rho_inf * (u_inf ** 2) * cf_comp
        u_tau = math.sqrt(tau_w / rho_w)
        
        # 4. Strict dimensional scaling for the required mesh element height
        target_y_plus = self.state.mesh_control.target_y_plus
        delta_y = (target_y_plus * mu_w) / (rho_w * u_tau)
        
        return delta_y

    def calculate_mesh_parameters(self) -> MeshState:
        """
        Constructs the hierarchical octree refinement architecture for CharLES Stitch.
        Applies necessary mathematical constraints and user-defined caps to generate 
        a safe, functional MeshState ready for template injection.
        """
        # Determine the physical wall distance required to achieve the target y+
        target_delta_y = self._estimate_wall_spacing()
        
        # Check for explicit user overrides in the YAML configuration
        user_hcp_delta = getattr(self.state.mesh_control, 'hcp_delta', None)
        user_max_level = getattr(self.state.mesh_control, 'max_refinement_level', None)
        
        # Define the base Hexagonal Close-Packed (HCP) far-field resolution
        if user_hcp_delta and user_hcp_delta > 0:
            hcp_delta = user_hcp_delta
        else:
            # Automatic heuristic: ensures at least 5 cells resolve the boundary layer thickness
            hcp_delta = min(5.0 * self.physics.estimated_delta99, self.physics.domain_height / 10.0)
            
        # Calculate the required binary bisection levels (Octree Refinement)
        if target_delta_y >= hcp_delta:
            max_level = 0
        else:
            max_level = math.ceil(math.log2(hcp_delta / target_delta_y))
            
        # Apply strict memory caps to prevent HPC node failure
        if user_max_level is not None:
            max_level = min(max_level, user_max_level)
        else:
            max_level = min(max_level, 9) # Absolute safety limit
        
        # Calculate the actual physical height of the smallest cell generated near the wall
        smallest_cell_size = hcp_delta / (2 ** max_level)
        
        # Enforce Quasi-2D volumetric constraints for CharLES Stitch
        # The generator requires a strictly isotropic spanwise boundary to process 2D
        domain_z = self.state.geometry.domain_z
        if domain_z > 0.0:
            half_w = domain_z / 2.0
        else:
            # Force the span to perfectly match one cell width to create a true quasi-2D slice
            half_w = smallest_cell_size / 2.0
            
        # Define the vertical extent where maximum refinement is strictly enforced
        refinement_height = 1.5 * self.physics.estimated_delta99
        
        print("\n=== Mesh Refinement Log ===")
        print(f"  Target Wall y+:            {self.state.mesh_control.target_y_plus:.2f}")
        print(f"  Required Spacing (Delta):  {target_delta_y:.6e}")
        print(f"  Far-field Cell Size:       {hcp_delta:.4f}")
        print(f"  Octree Refinement Level:   {max_level}")
        print(f"  Resulting Wall Cell Size:  {smallest_cell_size:.6e}")
        
        if smallest_cell_size <= target_delta_y:
            print("  Status: TARGET Y+ ACHIEVED")
        else:
            print("  Status: Y+ LIMITED BY MAX_REFINEMENT_LEVEL CONSTRAINT")
        print("===========================\n")
        
        return MeshState(
            hcp_delta=round(hcp_delta, 4),
            max_refinement_level=max_level,
            nlayers=getattr(self.state.mesh_control, 'nlayers', 8),
            nsmooth=getattr(self.state.mesh_control, 'nsmooth', 20),
            refinement_height=round(refinement_height, 4),
            half_w=round(half_w, 8),
            target_y_spacing=target_delta_y
        )