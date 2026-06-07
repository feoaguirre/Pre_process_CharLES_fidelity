"""
Module: physics_engine
Description: Handles fluid mechanics calculations, scaling transformations, and 
             boundary layer estimations. Implements multiple analytical and numerical 
             strategies (Eckert's Reference Temperature, Compressible Similarity Solutions)
             to accurately predict boundary layer growth and position the geometric structure.
"""

import math
import numpy as np
from scipy.integrate import solve_bvp
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any


@dataclass
class PhysicsState:
    """Data container for the calculated physical properties of the flow."""
    mu_inf: float
    p_inf: float
    t_aw: float                
    l_char: float              
    pregap_slip_length: float  
    pregap_noslip_length: float 
    postgap_length: float      
    domain_height: float       
    estimated_delta99: float   


class BoundaryLayerStrategy(ABC):
    """Abstract base class for calculating boundary layer development."""
    
    @abstractmethod
    def calculate_development_length(self, target_re: float, mach: float, gamma: float, 
                                     prandtl: float, t_inf: float, wall_bc: str, wall_T: float, state: Any) -> float:
        """Calculates the physical X-coordinate (L2) required to reach the target Reynolds number."""
        pass

    @abstractmethod
    def estimate_delta99(self, target_re: float, mach: float, gamma: float, 
                         prandtl: float, t_inf: float, wall_bc: str, wall_T: float, state: Any) -> float:
        """Estimates the 99% boundary layer thickness at the domain exit."""
        pass


class EckertAnalyticalStrategy(BoundaryLayerStrategy):
    """Calculates boundary layer growth using Eckert's Reference Temperature methodology."""
    
    def _get_t_star(self, mach, gamma, prandtl, t_inf, wall_bc, wall_T) -> tuple:
        """Calculates the adiabatic wall temperature and Eckert's reference temperature (T*)."""
        r = math.sqrt(prandtl)
        t_aw = t_inf * (1.0 + r * ((gamma - 1.0) / 2.0) * (mach ** 2))
        t_w = t_aw if wall_bc.upper() == "ADIABATIC" else (wall_T if wall_T > 0 else t_aw)
        
        # Eckert's Reference Temperature equation
        t_star = t_inf + 0.5 * (t_w - t_inf) + 0.22 * (t_aw - t_inf)
        return t_star, t_w

    def calculate_development_length(self, target_re, mach, gamma, prandtl, t_inf, wall_bc, wall_T, state) -> float:
        t_star, t_w = self._get_t_star(mach, gamma, prandtl, t_inf, wall_bc, wall_T)
        n = getattr(state.flow_physics, 'mu_power_law', 0.76)
        
        # Chapman-Rubesin constant based on the reference temperature
        c_star = (t_star / t_inf) ** (n - 1.0)
        
        # Analytical approximation for the Compressible Shape Factor (H)
        r = math.sqrt(prandtl)
        h_inc = 2.59 # Shape factor for incompressible Blasius boundary layer
        h_comp = (t_w / t_inf) * h_inc + r * ((gamma - 1.0) / 2.0) * (mach ** 2)
        
        # Exact relation for displacement thickness: delta* / x = I1 / sqrt(Re_x)
        i1 = h_comp * 0.664 * math.sqrt(c_star)
        
        # Since delta* = 1.0 (non-dimensionalization constraint), target_re represents Re_delta*
        # target_re = I1 * sqrt(Re_x) => Re_x = (target_re / I1)^2
        target_re_x = (target_re / i1) ** 2
        
        # Returns the development length (L2) in non-dimensional units of delta*
        return target_re_x / target_re

    def estimate_delta99(self, target_re, mach, gamma, prandtl, t_inf, wall_bc, wall_T, state) -> float:
        # Extract the position of the cavity's leading edge
        l2 = self.calculate_development_length(target_re, mach, gamma, prandtl, t_inf, wall_bc, wall_T, state)
        
        t_star, t_w = self._get_t_star(mach, gamma, prandtl, t_inf, wall_bc, wall_T)
        r = math.sqrt(prandtl)
        h_inc = 2.59
        h_comp = (t_w / t_inf) * h_inc + r * ((gamma - 1.0) / 2.0) * (mach ** 2)
        
        # Robust estimation of the physical ratio between delta99 and delta* at the leading edge
        ratio_99_to_star = (5.0 / 0.664) / h_comp * (t_star / t_inf)
        delta99_start = ratio_99_to_star * 1.0 # Since delta* = 1.0 by definition
        
        # Total length up to the physical domain exit
        L = state.geometry.structure_length
        # Fetches L3 dynamically; if it fails, assumes 20*L as a standard safety margin
        postgap_mult = getattr(state.domain_sizing, 'postgap_length_multiplier', 20.0)
        L3 = L * postgap_mult
        
        x_end = l2 + L + L3
        
        # Scales delta99 considering the Blasius boundary layer parabolic growth (~ sqrt(x))
        delta99_end = delta99_start * math.sqrt(x_end / l2)
        return delta99_end


class SimilaritySolutionStrategy(BoundaryLayerStrategy):
    """
    Calculates exact boundary layer growth by numerically solving the compressible 
    similarity equations utilizing the Illingworth-Stewartson transformation.
    """
    def _solve_compressible_blasius(self, mach, gamma, prandtl, wall_bc, wall_T_ratio):
        def ode_system(eta, y):
            f, fp, fpp, g, gp = y
            # Momentum equation: Exact transformed Blasius equation
            fppp = -0.5 * f * fpp
            # Energy equation: Includes the compressible viscous dissipation term
            gpp = -0.5 * prandtl * f * gp - prandtl * (gamma - 1.0) * (mach**2) * (fpp**2)
            return np.vstack((fp, fpp, fppp, gp, gpp))

        def bc(ya, yb):
            f0 = ya[0]
            fp0 = ya[1]
            if wall_bc.upper() == "ADIABATIC":
                gp0 = ya[4]
                wall_cond = gp0
            else:
                g0 = ya[3] - wall_T_ratio
                wall_cond = g0
                
            fp_inf = yb[1] - 1.0
            g_inf = yb[3] - 1.0
            
            return np.array([f0, fp0, wall_cond, fp_inf, g_inf])

        # Highly resolved grid (eta up to 15 ensures complete boundary layer convergence)
        eta = np.linspace(0, 15, 1000)
        y_guess = np.zeros((5, eta.size))
        y_guess[1, :] = 1.0 - np.exp(-eta)
        y_guess[2, :] = np.exp(-eta)
        
        # Initial guess based on the Crocco-Busemann relation (stabilizes the ODE solver)
        y_guess[3, :] = wall_T_ratio + (1.0 - wall_T_ratio + 0.5 * (gamma - 1.0) * mach**2) * y_guess[1, :] - 0.5 * (gamma - 1.0) * mach**2 * (y_guess[1, :]**2)
        
        res = solve_bvp(ode_system, bc, eta, y_guess, tol=1e-5, max_nodes=5000)
        if not res.success:
            return None
            
        fp = res.y[1]
        g = res.y[3]
        
        # Integral of the physical density displacement (g - f')
        integrand_star = g - fp
        i1 = np.trapezoid(integrand_star, x=res.x)
        
        # Calculation of the point where velocity reaches 99% of freestream
        idx_99 = np.searchsorted(fp, 0.99)
        if idx_99 < len(res.x):
            eta_99 = res.x[idx_99]
            i99 = np.trapezoid(g[:idx_99], x=res.x[:idx_99])
        else:
            i99 = np.trapezoid(g, x=res.x)
            
        # Returns the transformation integrals and the reached T_wall
        return i1, i99, g[0] 

    def calculate_development_length(self, target_re, mach, gamma, prandtl, t_inf, wall_bc, wall_T, state) -> float:
        wall_T_ratio = wall_T / t_inf if wall_T > 0 else 1.0
        sol = self._solve_compressible_blasius(mach, gamma, prandtl, wall_bc, wall_T_ratio)
        
        if sol is None:
            print("[WARNING] Similarity ODE solver failed. Falling back to Eckert analytical method.")
            return EckertAnalyticalStrategy().calculate_development_length(target_re, mach, gamma, prandtl, t_inf, wall_bc, wall_T, state)
            
        i1, _, t_w_ratio = sol
        n = getattr(state.flow_physics, 'mu_power_law', 0.76)
        
        # Calculates the exact Chapman-Rubesin C* from the ODE outputs
        r = math.sqrt(prandtl)
        t_aw_ratio = 1.0 + r * ((gamma - 1.0) / 2.0) * (mach ** 2)
        t_star_ratio = 0.5 * (1.0 + t_w_ratio) + 0.22 * (t_aw_ratio - 1.0)
        c_star = t_star_ratio ** (n - 1.0)
        
        # Exact Compressible Blasius relation
        target_re_x = (target_re / (i1 * math.sqrt(c_star))) ** 2
        return target_re_x / target_re

    def estimate_delta99(self, target_re, mach, gamma, prandtl, t_inf, wall_bc, wall_T, state) -> float:
        wall_T_ratio = wall_T / t_inf if wall_T > 0 else 1.0
        sol = self._solve_compressible_blasius(mach, gamma, prandtl, wall_bc, wall_T_ratio)
        
        if sol is None:
            return EckertAnalyticalStrategy().estimate_delta99(target_re, mach, gamma, prandtl, t_inf, wall_bc, wall_T, state)
            
        i1, i99, _ = sol
        
        # Since the problem is non-dimensionalized by delta* = 1.0, 
        # the physical ratio delta99 / delta* is exactly I_99 / I_1 in the transformation.
        delta99_start = (i99 / i1) * 1.0 
        
        l2 = self.calculate_development_length(target_re, mach, gamma, prandtl, t_inf, wall_bc, wall_T, state)
        
        L = state.geometry.structure_length
        postgap_mult = getattr(state.domain_sizing, 'postgap_length_multiplier', 20.0)
        L3 = L * postgap_mult
        x_end = l2 + L + L3
        
        # Strict scaling of the boundary layer growth (parabolic in x)
        delta99_end = delta99_start * math.sqrt(x_end / l2)
        return delta99_end


class FixedLengthStrategy(BoundaryLayerStrategy):
    """Bypasses physics calculations and uses a user-defined L2 length directly."""
    def calculate_development_length(self, target_re, mach, gamma, prandtl, t_inf, wall_bc, wall_T, state) -> float:
        return state.boundary_layer_setup.fixed_l2_length

    def estimate_delta99(self, target_re, mach, gamma, prandtl, t_inf, wall_bc, wall_T, state) -> float:
        # If L2 is fixed, we still need to estimate delta99 for the domain height (H).
        # We use the analytical method as a placeholder for the fixed domain.
        return EckertAnalyticalStrategy().estimate_delta99(target_re, mach, gamma, prandtl, t_inf, wall_bc, wall_T, state)


# =============================================================================
# FUTURE IMPLEMENTATION: Local DNS / RANS Strategy
# =============================================================================
class LocalDNSStrategy(BoundaryLayerStrategy):
    """
    [PLACEHOLDER]
    This class will dynamically generate a 2D mesh, run a fast CharLES/SU2 RANS or DNS
    simulation of a flat plate, extract the boundary layer profile, and iteratively
    find the exact physical X location where Re_delta_star matches the target.
    """
    def calculate_development_length(self, target_re, mach, gamma, prandtl, t_inf, wall_bc, wall_T, state) -> float:
        print("[INFO] Initializing Local DNS/RANS solver to determine L2 length...")
        # TODO 1: Generate 2D flat plate mesh up to Re_x_max
        # TODO 2: Create temporary charLES.in for 2D run
        # TODO 3: Execute solver using subprocess (mpirun -np X charles_ig.exe)
        # TODO 4: Load results, compute integral delta_star(x) along the plate
        # TODO 5: Interpolate to find X where Re_delta_star(X) == target_re
        raise NotImplementedError("Local DNS/RANS boundary layer calculation is not yet implemented.")

    def estimate_delta99(self, target_re, mach, gamma, prandtl, t_inf, wall_bc, wall_T, state) -> float:
        # TODO: Extract delta99 directly from the velocity profile obtained in calculate_development_length
        raise NotImplementedError("Local DNS/RANS boundary layer calculation is not yet implemented.")
# =============================================================================


class PhysicsEngine:
    """Master controller that initializes the requested boundary layer physics strategy."""
    def __init__(self, state: Any):
        self.state = state
        self.bl_strategy = self._select_strategy()

    def _select_strategy(self) -> BoundaryLayerStrategy:
        """Factory method to instantiate the chosen boundary layer strategy from YAML."""
        method = self.state.boundary_layer_setup.method.lower()
        if method == "similarity_solution":
            return SimilaritySolutionStrategy()
        elif method == "fixed_length":
            return FixedLengthStrategy()
        elif method == "local_dns":
            return LocalDNSStrategy()
        else:
            return EckertAnalyticalStrategy() # Default fallback

    def calculate_derived_properties(self) -> PhysicsState:
        """Executes the chosen strategy and constructs the physical state container."""
        phys = self.state.flow_physics
        geom = self.state.geometry
        bc = self.state.boundary_conditions
        ds = self.state.domain_sizing
        
        p_inf = (phys.inflow_rho * (phys.inflow_u ** 2)) / (phys.gamma * (phys.mach_number ** 2))
        mu_inf = (phys.inflow_rho * phys.inflow_u * 1.0) / phys.target_reynolds
        
        r = math.sqrt(phys.prandtl)
        t_aw = phys.inflow_T * (1.0 + r * ((phys.gamma - 1.0) / 2.0) * (phys.mach_number ** 2))
        
        # Dynamically calculates L2 and delta99 based on the chosen strategy in the YAML
        l2_distance = self.bl_strategy.calculate_development_length(
            phys.target_reynolds, phys.mach_number, phys.gamma, phys.prandtl, 
            phys.inflow_T, bc.wall_bc, bc.wall_T, self.state
        )
        
        delta99 = self.bl_strategy.estimate_delta99(
            phys.target_reynolds, phys.mach_number, phys.gamma, phys.prandtl, 
            phys.inflow_T, bc.wall_bc, bc.wall_T, self.state
        )
        
        l1_distance = ds.l1_fixed_length
        
        if ds.l3_definition == "multiple_of_L":
            l3_distance = ds.l3_value * geom.structure_length
        elif ds.l3_definition == "fixed_end_x":
            l3_distance = ds.l3_value - geom.structure_length
        else:
            raise ValueError(f"Unknown l3_definition: {ds.l3_definition}")
            
        if ds.h_definition == "multiple_of_delta99":
            domain_height = ds.h_value * delta99
        elif ds.h_definition == "fixed_end_y":
            domain_height = ds.h_value
        else:
            raise ValueError(f"Unknown h_definition: {ds.h_definition}")
        
        scale_factor = 1.0 
            
        return PhysicsState(
            mu_inf=mu_inf * scale_factor,
            p_inf=p_inf,
            t_aw=t_aw,
            l_char=scale_factor,
            pregap_slip_length=l1_distance * scale_factor,
            pregap_noslip_length=l2_distance * scale_factor,
            postgap_length=l3_distance * scale_factor,
            domain_height=domain_height * scale_factor,
            estimated_delta99=delta99 * scale_factor
        )