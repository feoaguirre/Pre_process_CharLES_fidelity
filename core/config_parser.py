"""
Module: config_parser
Description: Parses the simulation_config.yaml file into strongly typed, 
             immutable Python dataclasses. Handles automated 2D grid validation
             and nested probe spacing configurations. 
             
             By using @dataclass(frozen=True), we ensure that once the configuration
             is loaded, no other part of the code can accidentally change a value 
             (e.g., changing Mach number halfway through the mesh generation).
"""

import yaml
from dataclasses import dataclass, field
from typing import Dict, List, Optional

# =============================================================================
# INDIVIDUAL CONFIGURATION BLOCKS
# =============================================================================

@dataclass(frozen=True)
class IdentityConfig:
    """Basic identification and core solver decisions."""
    run_name: str                # Name of the output folder and PBS job
    solver_type: str             # "DNS" or "LES" (influences mesh y+ target and SGS models)
    flow_regime: str             # "subsonic", "supersonic", or "hypersonic"
    dimensionalization_base: str # E.g., "delta_star". Used by the physics engine to scale outputs.

@dataclass(frozen=True)
class GeometryConfig:
    """Physical dimensions of the structure being simulated."""
    structure_type: str          # "symmetric_cavity", "step", etc.
    structure_length: float      # L
    structure_depth: float       # D
    span_z: float                # Span of the structure (0.0 for 2D)
    domain_z: float              # Span of the entire domain. If 0.0, parser assumes quasi-2D.

@dataclass(frozen=True)
class DomainSizingConfig:
    """Rules for how the computational domain boundaries are calculated."""
    l1_fixed_length: float       # Physical length of the slip wall before the no-slip plate
    l3_definition: str           # Rule for the wake region ("multiple_of_L" or "fixed_end_x")
    l3_value: float              # Multiplier or exact coordinate based on l3_definition
    h_definition: str            # Rule for domain height ("multiple_of_delta99" or "fixed_end_y")
    h_value: float               # Multiplier or exact coordinate based on h_definition

@dataclass
class MeshControlConfig:
    target_y_plus: float
    # Add these new lines:
    hcp_delta: float
    max_refinement_level: int
    # Keep your existing fields:
    nlayers: int
    nsmooth: int

@dataclass(frozen=True)
class FlowPhysicsConfig:
    """Thermodynamic and aerodynamic properties of the freestream flow."""
    target_reynolds: float       # Desired Re at the structure's leading edge
    mach_number: float           # Freestream Mach (M_inf)
    inflow_u: float              # Non-dimensional velocity (usually 1.0)
    inflow_rho: float            # Non-dimensional density (usually 1.0)
    inflow_T: float              # Non-dimensional temperature (usually 1.0)
    gamma: float                 # Heat capacity ratio (Cp/Cv, typically 1.4 for air)
    prandtl: float               # Prandtl number (typically 0.72 for air)
    mu_power_law: float          # Viscosity power-law exponent (e.g., 0.76)

@dataclass(frozen=True)
class BoundaryLayerSetupConfig:
    """Configuration for calculating boundary layer growth to position the structure."""
    method: str                  # "eckert", "similarity_solution", "fixed_length", or "local_dns"
    fixed_l2_length: float       # Absolute distance L2 (used only if method == "fixed_length")

@dataclass(frozen=True)
class BoundaryConditionsConfig:
    """Thermal conditions applied to the solid walls."""
    wall_bc: str                 # "ADIABATIC" (q_w = 0) or "ISOTHERMAL" (T_w = constant)
    wall_T: float                # Wall temperature ratio (ignored if wall_bc is ADIABATIC)
    z_boundaries: str            # "PERIODIC" or "SYMMETRY" (used for 3D spans)

@dataclass(frozen=True)
class SimulationControlConfig:
    """Solver time-stepping and high-performance computing (HPC) settings."""
    pbs_queue: str               # Target Zeus cluster queue (e.g., mafat_new_q)
    transient_simtime_ftt: float # Flow Through Times to flush initial transients
    steady_simtime_ftt: float    # Flow Through Times for statistical data gathering
    cfl: float                   # Courant–Friedrichs–Lewy limit for time stepping

# =============================================================================
# PROBES & DATA EXTRACTION BLOCKS
# =============================================================================

@dataclass(frozen=True)
class Spacing1DConfig:
    """Defines how points are distributed along a single axis (1D)."""
    type: str                                # "uniform", "logarithmic", "exponential", "mesh_like", "custom"
    points: Optional[int] = None             # Number of points (used by uniform/log/exp)
    value: Optional[float] = None            # Multiplier factor (used by mesh_like)
    custom_vector: List[float] = field(default_factory=list) # Explicit coordinate array (if type == "custom")

@dataclass(frozen=True)
class SpaceProbeRegionConfig:
    """Defines a 3D bounding box and internal spacing for spatial probes (snapshots)."""
    x_bounds: List[float]        # [start_x, end_x]. If start==end, the axis is fixed.
    y_bounds: List[float]        # [start_y, end_y]
    z_bounds: List[float]        # [start_z, end_z]
    x_spacing: Spacing1DConfig   # Point distribution rule for X
    y_spacing: Spacing1DConfig   # Point distribution rule for Y
    z_spacing: Spacing1DConfig   # Point distribution rule for Z

@dataclass(frozen=True)
class TimeProbeRegionConfig:
    """Defines spatial bounds and explicit Z-planes for high-frequency time probes."""
    x_bounds: List[float]
    y_bounds: List[float]
    z_planes: List[float]        # Explicit list of Z coordinates where the 2D grid will be cloned
    x_spacing: Spacing1DConfig
    y_spacing: Spacing1DConfig

@dataclass(frozen=True)
class IOAndProbesConfig:
    """Master controller for simulation outputs and sensor intervals."""
    check_interval_steps: int
    image_interval_steps: int
    space_probes_write_interval: int
    time_probes_write_interval: int
    space_probes: Dict[str, SpaceProbeRegionConfig]
    time_probes: Dict[str, TimeProbeRegionConfig]

# =============================================================================
# MASTER STATE & PARSING LOGIC
# =============================================================================

@dataclass(frozen=True)
class SimulationState:
    """
    Master data container for the current simulation configuration.
    This object is passed to all other modules (Physics, Mesh, Templates)
    to serve as the single source of truth.
    """
    identity: IdentityConfig
    geometry: GeometryConfig
    domain_sizing: DomainSizingConfig
    mesh_control: MeshControlConfig          # <-- DEVE ESTAR AQUI
    flow_physics: FlowPhysicsConfig
    boundary_layer_setup: BoundaryLayerSetupConfig
    boundary_conditions: BoundaryConditionsConfig
    simulation_control: SimulationControlConfig
    io_and_probes: IOAndProbesConfig

    @classmethod
    def from_yaml(cls, filepath: str) -> 'SimulationState':
        """
        Loads the YAML file, validates edge cases (like 2D assumptions),
        and strongly types every parameter into nested dataclasses.
        """
        # Open and load the YAML file into a standard Python dictionary
        with open(filepath, 'r') as file:
            raw_data = yaml.safe_load(file)

        # Map simple top-level dictionaries to their respective dataclasses
        identity = IdentityConfig(**raw_data.get('identity', {}))
        
        # Geometry Processing (Includes automated Quasi-2D validation)
        geom_dict = raw_data.get('geometry', {})
        if geom_dict.get('domain_z', 0.0) == 0.0:
            print("[INFO] domain_z is set to 0.0. Processing as a quasi-2D simulation setup.")
        geometry = GeometryConfig(**geom_dict)
        
        # Load mid-level dictionaries
        domain_sizing = DomainSizingConfig(**raw_data.get('domain_sizing', {}))
        mesh_control = MeshControlConfig(**raw_data.get('mesh_control', {})) # <-- CORRIGIDO AQUI
        flow_physics = FlowPhysicsConfig(**raw_data.get('flow_physics', {}))
        boundary_layer_setup = BoundaryLayerSetupConfig(**raw_data.get('boundary_layer_setup', {}))
        boundary_conditions = BoundaryConditionsConfig(**raw_data.get('boundary_conditions', {}))
        simulation_control = SimulationControlConfig(**raw_data.get('simulation_control', {}))
        
        # ---------------------------------------------------------------------
        # Unpacking Nested Probe Configurations
        # ---------------------------------------------------------------------
        io_data = raw_data.get('io_and_probes', {})
        
        # Parse Space Probes
        space_probes_parsed = {}
        for region_name, p_data in io_data.get('space_probes', {}).items():
            space_probes_parsed[region_name] = SpaceProbeRegionConfig(
                x_bounds=p_data.get('x_bounds', [0.0, 0.0]),
                y_bounds=p_data.get('y_bounds', [0.0, 0.0]),
                z_bounds=p_data.get('z_bounds', [0.0, 0.0]),
                x_spacing=Spacing1DConfig(**p_data.get('x_spacing', {})),
                y_spacing=Spacing1DConfig(**p_data.get('y_spacing', {})),
                z_spacing=Spacing1DConfig(**p_data.get('z_spacing', {}))
            )
            
        # Parse Time Probes
        time_probes_parsed = {}
        for region_name, p_data in io_data.get('time_probes', {}).items():
            time_probes_parsed[region_name] = TimeProbeRegionConfig(
                x_bounds=p_data.get('x_bounds', [0.0, 0.0]),
                y_bounds=p_data.get('y_bounds', [0.0, 0.0]),
                z_planes=p_data.get('z_planes', [0.0]),
                x_spacing=Spacing1DConfig(**p_data.get('x_spacing', {})),
                y_spacing=Spacing1DConfig(**p_data.get('y_spacing', {}))
            )

        # Consolidate IO and Probes
        io_and_probes = IOAndProbesConfig(
            check_interval_steps=io_data.get('check_interval_steps', 1000),
            image_interval_steps=io_data.get('image_interval_steps', 3000),
            space_probes_write_interval=io_data.get('space_probes_write_interval', 10000),
            time_probes_write_interval=io_data.get('time_probes_write_interval', 100),
            space_probes=space_probes_parsed,
            time_probes=time_probes_parsed
        )

        # Return the fully constructed, immutable master state
        return cls(
            identity=identity,
            geometry=geometry,
            domain_sizing=domain_sizing,
            mesh_control=mesh_control,       # <-- CORRIGIDO AQUI
            flow_physics=flow_physics,
            boundary_layer_setup=boundary_layer_setup,
            boundary_conditions=boundary_conditions,
            simulation_control=simulation_control,
            io_and_probes=io_and_probes
        )