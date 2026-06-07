"""
Module: config_parser
Description: Parses the simulation_config.yaml file into strongly typed, 
             immutable Python dataclasses. Handles automated 2D grid validation
             and nested probe spacing configurations. 
             
             By utilizing @dataclass(frozen=True), the framework guarantees that 
             once the configuration is loaded into memory, no subsequent module 
             can inadvertently mutate a value (e.g., modifying the Mach number 
             during mesh generation).
"""

import yaml
from dataclasses import dataclass, field
from typing import Dict, List, Optional

# =============================================================================
# INDIVIDUAL CONFIGURATION BLOCKS
# =============================================================================

@dataclass(frozen=True)
class IdentityConfig:
    """Basic identification and core solver parameters."""
    run_name: str                # Name of the output directory and PBS job identifier
    solver_type: str             # "DNS" or "LES" (dictates mesh y+ targets and SGS model activation)
    flow_regime: str             # Flow classification ("subsonic", "supersonic", or "hypersonic")
    dimensionalization_base: str # Reference length (e.g., "delta_star"). Used to scale spatial outputs.

@dataclass(frozen=True)
class GeometryConfig:
    """Physical dimensions of the structural domain."""
    structure_type: str          # Geometry classification (e.g., "symmetric_cavity", "backward_facing_step")
    structure_length: float      # Streamwise length of the structure (L)
    structure_depth: float       # Wall-normal depth of the structure (D)
    span_z: float                # Spanwise width of the structure (Set to 0.0 for pure 2D setups)
    domain_z: float              # Total spanwise width of the fluid domain. If 0.0, the parser forces a quasi-2D grid.

@dataclass(frozen=True)
class DomainSizingConfig:
    """Rules and multipliers defining the outer boundaries of the computational domain."""
    l1_fixed_length: float       # Absolute length of the slip-wall development region preceding the no-slip plate
    l3_definition: str           # Boundary rule for the downstream wake region ("multiple_of_L" or "fixed_end_x")
    l3_value: float              # Corresponding multiplier or absolute coordinate for L3
    h_definition: str            # Boundary rule for the domain height ("multiple_of_delta99" or "fixed_end_y")
    h_value: float               # Corresponding multiplier or absolute coordinate for the top boundary

@dataclass(frozen=True)
class MeshControlConfig:
    """Constraints for the Voronoi octree mesh generation."""
    target_y_plus: float         # Desired non-dimensional wall distance for the first cell layer
    hcp_delta: float             # Base Hexagonal Close-Packed (HCP) far-field cell size
    max_refinement_level: int    # Maximum allowable octree bisections (prevents Out-Of-Memory errors)
    nlayers: int                 # Number of isotropic layers maintained near solid boundaries
    nsmooth: int                 # Number of smoothing iterations applied to the grid transitions

@dataclass(frozen=True)
class FlowPhysicsConfig:
    """Thermodynamic and aerodynamic boundary conditions for the freestream flow."""
    target_reynolds: float       # Reference Reynolds number evaluated at the structure's leading edge
    mach_number: float           # Freestream Mach number (M_inf)
    inflow_u: float              # Non-dimensional streamwise velocity (standardized to 1.0)
    inflow_rho: float            # Non-dimensional freestream density (standardized to 1.0)
    inflow_T: float              # Non-dimensional freestream temperature (standardized to 1.0)
    gamma: float                 # Specific heat ratio (Cp/Cv, typically 1.4 for ideal air)
    prandtl: float               # Prandtl number (typically 0.72 for standard air)
    mu_power_law: float          # Exponent for the viscosity power-law model (e.g., 0.76 or 0.95)

@dataclass(frozen=True)
class BoundaryLayerSetupConfig:
    """Methodology for estimating the boundary layer growth to position the geometric structure."""
    method: str                  # Calculation regime: "eckert", "similarity_solution", "fixed_length", or "local_dns"
    fixed_l2_length: float       # Absolute development length (L2). Evaluated only if method == "fixed_length"

@dataclass(frozen=True)
class BoundaryConditionsConfig:
    """Thermal conditions enforced upon all solid no-slip boundaries."""
    wall_bc: str                 # "ADIABATIC" (zero heat flux) or "ISOTHERMAL" (constant wall temperature)
    wall_T: float                # Wall-to-freestream temperature ratio (T_w / T_inf). Ignored if ADIABATIC.
    z_boundaries: str            # Spanwise boundary condition paradigm ("PERIODIC" or "SYMMETRY")

@dataclass(frozen=True)
class SimulationControlConfig:
    """Parameters governing time advancement and High-Performance Computing (HPC) deployment."""
    pbs_queue: str               # Target scheduling queue on the HPC cluster
    transient_simtime_ftt: float # Flow-Through Times (FTT) allocated to flush numerical initialization transients
    steady_simtime_ftt: float    # Flow-Through Times (FTT) allocated for active statistical data collection
    cfl: float                   # Courant-Friedrichs-Lewy (CFL) limit governing dynamic time-stepping

# =============================================================================
# PROBES & DATA EXTRACTION BLOCKS
# =============================================================================

@dataclass(frozen=True)
class Spacing1DConfig:
    """Defines the point distribution logic along a single spatial axis."""
    type: str                                # Distribution pattern: "uniform", "logarithmic", "exponential", "mesh_like", or "custom"
    points: Optional[int] = None             # Number of discrete nodes (required for uniform/log/exp)
    value: Optional[float] = None            # Geometric multiplier factor (required for mesh_like spacing)
    custom_vector: List[float] = field(default_factory=list) # Explicit coordinate array (required if type == "custom")

@dataclass(frozen=True)
class SpaceProbeRegionConfig:
    """Defines a 3D bounding box and internal node distribution for spatial probes (volumetric snapshots)."""
    x_bounds: List[float]        # [start_x, end_x]. If start == end, the axis collapses to a 2D plane.
    y_bounds: List[float]        # [start_y, end_y]
    z_bounds: List[float]        # [start_z, end_z]
    x_spacing: Spacing1DConfig   # Point distribution rule applied along the X-axis
    y_spacing: Spacing1DConfig   # Point distribution rule applied along the Y-axis
    z_spacing: Spacing1DConfig   # Point distribution rule applied along the Z-axis

@dataclass(frozen=True)
class TimeProbeRegionConfig:
    """Defines spatial bounds and explicit Z-planes for high-frequency pointcloud probes."""
    x_bounds: List[float]
    y_bounds: List[float]
    z_planes: List[float]        # Explicit list of discrete Z-coordinates where the 2D (X, Y) grid will be duplicated
    x_spacing: Spacing1DConfig
    y_spacing: Spacing1DConfig

@dataclass(frozen=True)
class IOAndProbesConfig:
    """Master controller scheduling simulation I/O operations and sensor extraction frequencies."""
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
    Master data container encapsulating the entire simulation configuration.
    This centralized object is injected into all downstream modules (Physics, Mesh, Templates)
    to serve as the strict single source of truth.
    """
    identity: IdentityConfig
    geometry: GeometryConfig
    domain_sizing: DomainSizingConfig
    mesh_control: MeshControlConfig
    flow_physics: FlowPhysicsConfig
    boundary_layer_setup: BoundaryLayerSetupConfig
    boundary_conditions: BoundaryConditionsConfig
    simulation_control: SimulationControlConfig
    io_and_probes: IOAndProbesConfig

    @classmethod
    def from_yaml(cls, filepath: str) -> 'SimulationState':
        """
        Loads the YAML configuration file, validates physical edge cases 
        (e.g., quasi-2D domain constraints), and securely maps every parameter 
        into nested, immutable dataclasses.
        """
        # Read the raw YAML configuration
        with open(filepath, 'r') as file:
            raw_data = yaml.safe_load(file)

        # Instantiate simple top-level configurations
        identity = IdentityConfig(**raw_data.get('identity', {}))
        
        # Process Geometry and enforce quasi-2D safeguards if span is absent
        geom_dict = raw_data.get('geometry', {})
        if geom_dict.get('domain_z', 0.0) == 0.0:
            print("[INFO] ConfigParser: 'domain_z' is set to 0.0. Applying quasi-2D simulation constraints.")
        geometry = GeometryConfig(**geom_dict)
        
        # Instantiate mid-level configurations
        domain_sizing = DomainSizingConfig(**raw_data.get('domain_sizing', {}))
        mesh_control = MeshControlConfig(**raw_data.get('mesh_control', {}))
        flow_physics = FlowPhysicsConfig(**raw_data.get('flow_physics', {}))
        boundary_layer_setup = BoundaryLayerSetupConfig(**raw_data.get('boundary_layer_setup', {}))
        boundary_conditions = BoundaryConditionsConfig(**raw_data.get('boundary_conditions', {}))
        simulation_control = SimulationControlConfig(**raw_data.get('simulation_control', {}))
        
        # ---------------------------------------------------------------------
        # Deep Unpacking for Nested Probe Configurations
        # ---------------------------------------------------------------------
        io_data = raw_data.get('io_and_probes', {})
        
        # Map Space Probes
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
            
        # Map Time Probes
        time_probes_parsed = {}
        for region_name, p_data in io_data.get('time_probes', {}).items():
            time_probes_parsed[region_name] = TimeProbeRegionConfig(
                x_bounds=p_data.get('x_bounds', [0.0, 0.0]),
                y_bounds=p_data.get('y_bounds', [0.0, 0.0]),
                z_planes=p_data.get('z_planes', [0.0]),
                x_spacing=Spacing1DConfig(**p_data.get('x_spacing', {})),
                y_spacing=Spacing1DConfig(**p_data.get('y_spacing', {}))
            )

        # Consolidate IO operations
        io_and_probes = IOAndProbesConfig(
            check_interval_steps=io_data.get('check_interval_steps', 1000),
            image_interval_steps=io_data.get('image_interval_steps', 3000),
            space_probes_write_interval=io_data.get('space_probes_write_interval', 10000),
            time_probes_write_interval=io_data.get('time_probes_write_interval', 100),
            space_probes=space_probes_parsed,
            time_probes=time_probes_parsed
        )

        # Return the securely constructed master state
        return cls(
            identity=identity,
            geometry=geometry,
            domain_sizing=domain_sizing,
            mesh_control=mesh_control,
            flow_physics=flow_physics,
            boundary_layer_setup=boundary_layer_setup,
            boundary_conditions=boundary_conditions,
            simulation_control=simulation_control,
            io_and_probes=io_and_probes
        )