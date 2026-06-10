"""
Module: config_parser
Description: Parses the simulation_config.yaml file into strongly typed, 
             immutable Python dataclasses. Handles automated 2D grid validation,
             dynamic SGS model allocation for LES/DNS, granular probe spacing, 
             variable selection, image generation, and extraction frequencies.
             
             By utilizing @dataclass(frozen=True), the framework guarantees that 
             once the configuration is loaded into memory, no subsequent module 
             can inadvertently mutate a core value.
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
    run_name: str
    solver_type: str             # "DNS" or "LES"
    flow_regime: str             # "subsonic", "supersonic", or "hypersonic"
    dimensionalization_base: str 
    sgs_model: str               # Sub-Grid Scale model (Enforced to "NONE" if DNS)

@dataclass(frozen=True)
class GeometryConfig:
    """Physical dimensions of the structural domain."""
    structure_type: str
    structure_length: float
    structure_depth: float
    span_z: float
    domain_z: float              # 0.0 triggers automatic quasi-2D enforcement

@dataclass(frozen=True)
class DomainSizingConfig:
    """Rules and multipliers defining the outer boundaries of the computational domain."""
    l1_fixed_length: float
    l3_definition: str
    l3_value: float
    h_definition: str
    h_value: float

@dataclass(frozen=True)
class MeshControlConfig:
    """Constraints for the Voronoi octree mesh generation."""
    target_y_plus: float
    hcp_delta: float
    max_refinement_level: int
    nlayers: int
    nsmooth: int

@dataclass(frozen=True)
class FlowPhysicsConfig:
    """Thermodynamic and aerodynamic boundary conditions for the freestream flow."""
    target_reynolds: float
    mach_number: float
    inflow_u: float
    inflow_rho: float
    inflow_T: float
    gamma: float
    prandtl: float
    mu_power_law: float

@dataclass(frozen=True)
class BoundaryLayerSetupConfig:
    """Methodology for estimating the boundary layer growth."""
    method: str
    fixed_l2_length: float

@dataclass(frozen=True)
class BoundaryConditionsConfig:
    """Thermal conditions enforced upon all solid no-slip boundaries."""
    wall_bc: str                 # "ADIABATIC" or "ISOTHERMAL"
    wall_T: float
    z_boundaries: str

@dataclass(frozen=True)
class SimulationControlConfig:
    """Parameters governing time advancement and High-Performance Computing (HPC) deployment."""
    pbs_queue: str
    transient_simtime_ftt: float
    steady_simtime_ftt: float
    cfl: float

# =============================================================================
# ADVANCED PROBING & I/O ARCHITECTURE
# =============================================================================

@dataclass(frozen=True)
class Spacing1DConfig:
    """Defines the point distribution logic along a single spatial axis."""
    type: str                                # "uniform", "logarithmic", "exponential", "mesh_multiple", etc.
    points: Optional[int] = None
    value: Optional[float] = None
    custom_vector: List[float] = field(default_factory=list)

@dataclass(frozen=True)
class BaseProbeConfig:
    """Standard spatial/temporal probe bounding box with granular extraction controls."""
    write_interval: int
    variables: List[str]
    x_bounds: List[float]; y_bounds: List[float]
    x_spacing: Spacing1DConfig; y_spacing: Spacing1DConfig

@dataclass(frozen=True)
class SpaceProbeRegionConfig(BaseProbeConfig):
    z_bounds: List[float]
    z_spacing: Spacing1DConfig

@dataclass(frozen=True)
class TimeProbeRegionConfig(BaseProbeConfig):
    z_planes: List[float]

@dataclass(frozen=True)
class MeanFlowCfConfig(BaseProbeConfig):
    """Skin Friction (Cf) probes tracking du/dy gradients near Y~0."""
    z_planes: List[float]

@dataclass(frozen=True)
class MeanFlowBLConfig(BaseProbeConfig):
    """Boundary Layer profiles extracting vertical lines to compute integral thicknesses."""
    z_planes: List[float]

@dataclass(frozen=True)
class LargeDataBlockConfig:
    """SPOD/POD volumetric blocks featuring split transient/steady output frequencies."""
    transient_write_interval: int
    steady_write_interval: int
    variables: List[str]
    x_bounds: List[float]; y_bounds: List[float]; z_bounds: List[float]
    x_spacing: Spacing1DConfig; y_spacing: Spacing1DConfig; z_spacing: Spacing1DConfig

@dataclass(frozen=True)
class ImageOutputConfig:
    """Defines native CharLES 2D slice image exports for real-time visualization."""
    variable: str
    write_interval: int

@dataclass(frozen=True)
class IOAndProbesConfig:
    """Master controller scheduling simulation IO and complex sensor registries."""
    check_interval_steps: int
    save_full_mesh: bool
    full_mesh_write_interval: int
    
    image_outputs: Dict[str, ImageOutputConfig]
    space_probes: Dict[str, SpaceProbeRegionConfig]
    time_probes: Dict[str, TimeProbeRegionConfig]
    cf_probes: Dict[str, MeanFlowCfConfig]
    boundary_layer_probes: Dict[str, MeanFlowBLConfig]
    large_data_probes: Dict[str, LargeDataBlockConfig]

# =============================================================================
# MASTER STATE & LOGIC LAYER
# =============================================================================

@dataclass(frozen=True)
class SimulationState:
    """Master immutable data container encapsulating the entire simulation architecture."""
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
        """Parses the YAML, applies core aerospace logic rules (e.g., DNS vs LES SGS handling)."""
        with open(filepath, 'r') as file:
            raw_data = yaml.safe_load(file)

        # ---------------------------------------------------------------------
        # 1. Identity & Dynamic SGS Model Logic
        # ---------------------------------------------------------------------
        identity_raw = raw_data.get('identity', {})
        solver_type = identity_raw.get('solver_type', 'DNS').upper()
        user_sgs = identity_raw.get('sgs_model')

        if solver_type == "LES":
            # Business Rule: If LES is requested but no model is provided, default to VREMAN
            assigned_sgs = str(user_sgs).upper() if user_sgs else "VREMAN"
            print(f"[INFO] ConfigParser: LES Solver active. Applied SGS Model: {assigned_sgs}")
        else:
            # Business Rule: DNS natively requires SGS to be disabled
            assigned_sgs = "NONE"
            if user_sgs and str(user_sgs).upper() != "NONE":
                print(f"[WARNING] ConfigParser: DNS Solver active. Forcing SGS Model to NONE (Ignored '{user_sgs}').")

        identity = IdentityConfig(
            run_name=identity_raw.get('run_name', 'default_run'),
            solver_type=solver_type,
            flow_regime=identity_raw.get('flow_regime', 'hypersonic'),
            dimensionalization_base=identity_raw.get('dimensionalization_base', 'delta_star'),
            sgs_model=assigned_sgs
        )
        
        # ---------------------------------------------------------------------
        # 2. Geometry & Domain Mapping
        # ---------------------------------------------------------------------
        geom_dict = raw_data.get('geometry', {})
        if geom_dict.get('domain_z', 0.0) == 0.0:
            print("[INFO] ConfigParser: 'domain_z' is 0.0. Processing as a Quasi-2D setup.")
        geometry = GeometryConfig(**geom_dict)
        
        domain_sizing = DomainSizingConfig(**raw_data.get('domain_sizing', {}))
        mesh_control = MeshControlConfig(**raw_data.get('mesh_control', {}))
        flow_physics = FlowPhysicsConfig(**raw_data.get('flow_physics', {}))
        boundary_layer_setup = BoundaryLayerSetupConfig(**raw_data.get('boundary_layer_setup', {}))
        boundary_conditions = BoundaryConditionsConfig(**raw_data.get('boundary_conditions', {}))
        simulation_control = SimulationControlConfig(**raw_data.get('simulation_control', {}))
        
        # ---------------------------------------------------------------------
        # 3. Granular Probe & Image Architecture Unpacking
        # ---------------------------------------------------------------------
        io_data = raw_data.get('io_and_probes', {})
        
        def _parse_spacing(sp_data: dict) -> Spacing1DConfig:
            return Spacing1DConfig(**sp_data) if sp_data else Spacing1DConfig(type="uniform", points=1)

        image_outputs_parsed = {}
        for name, p_data in io_data.get('image_outputs', {}).items():
            image_outputs_parsed[name] = ImageOutputConfig(
                variable=p_data.get('variable', 'rho'),
                write_interval=p_data.get('write_interval', 10000)
            )

        space_probes_parsed = {}
        for name, p_data in io_data.get('space_probes', {}).items():
            space_probes_parsed[name] = SpaceProbeRegionConfig(
                write_interval=p_data.get('write_interval', 10000),
                variables=p_data.get('variables', ["comp(u,0)", "comp(u,1)", "p", "rho"]),
                x_bounds=p_data.get('x_bounds', [0.0, 0.0]), y_bounds=p_data.get('y_bounds', [0.0, 0.0]), z_bounds=p_data.get('z_bounds', [0.0, 0.0]),
                x_spacing=_parse_spacing(p_data.get('x_spacing')), y_spacing=_parse_spacing(p_data.get('y_spacing')), z_spacing=_parse_spacing(p_data.get('z_spacing'))
            )
            
        time_probes_parsed = {}
        for name, p_data in io_data.get('time_probes', {}).items():
            time_probes_parsed[name] = TimeProbeRegionConfig(
                write_interval=p_data.get('write_interval', 1000),
                variables=p_data.get('variables', ["comp(u,0)", "comp(u,1)", "p", "rho"]),
                x_bounds=p_data.get('x_bounds', [0.0, 0.0]), y_bounds=p_data.get('y_bounds', [0.0, 0.0]), z_planes=p_data.get('z_planes', [0.0]),
                x_spacing=_parse_spacing(p_data.get('x_spacing')), y_spacing=_parse_spacing(p_data.get('y_spacing'))
            )

        cf_probes_parsed = {}
        for name, p_data in io_data.get('cf_probes', {}).items():
            cf_probes_parsed[name] = MeanFlowCfConfig(
                write_interval=p_data.get('write_interval', 50000),
                variables=p_data.get('variables', ["comp(u,0)", "p", "rho", "T"]),
                x_bounds=p_data.get('x_bounds', [0.0, 0.0]), y_bounds=p_data.get('y_bounds', [0.0, 0.0]), z_planes=p_data.get('z_planes', [0.0]),
                x_spacing=_parse_spacing(p_data.get('x_spacing')), y_spacing=_parse_spacing(p_data.get('y_spacing'))
            )

        bl_probes_parsed = {}
        for name, p_data in io_data.get('boundary_layer_probes', {}).items():
            bl_probes_parsed[name] = MeanFlowBLConfig(
                write_interval=p_data.get('write_interval', 50000),
                variables=p_data.get('variables', ["comp(u,0)", "p", "rho", "T"]),
                x_bounds=p_data.get('x_bounds', [0.0, 0.0]), y_bounds=p_data.get('y_bounds', [0.0, 0.0]), z_planes=p_data.get('z_planes', [0.0]),
                x_spacing=_parse_spacing(p_data.get('x_spacing')), y_spacing=_parse_spacing(p_data.get('y_spacing'))
            )

        large_data_probes_parsed = {}
        for name, p_data in io_data.get('large_data_probes', {}).items():
            large_data_probes_parsed[name] = LargeDataBlockConfig(
                transient_write_interval=p_data.get('transient_write_interval', 100000),
                steady_write_interval=p_data.get('steady_write_interval', 500),
                variables=p_data.get('variables', ["comp(u,0)", "comp(u,1)", "p", "rho"]),
                x_bounds=p_data.get('x_bounds', [0.0, 0.0]), y_bounds=p_data.get('y_bounds', [0.0, 0.0]), z_bounds=p_data.get('z_bounds', [0.0, 0.0]),
                x_spacing=_parse_spacing(p_data.get('x_spacing')), y_spacing=_parse_spacing(p_data.get('y_spacing')), z_spacing=_parse_spacing(p_data.get('z_spacing'))
            )

        io_and_probes = IOAndProbesConfig(
            check_interval_steps=io_data.get('check_interval_steps', 1000),
            save_full_mesh=io_data.get('save_full_mesh', False),
            full_mesh_write_interval=io_data.get('full_mesh_write_interval', 10000),
            image_outputs=image_outputs_parsed,
            space_probes=space_probes_parsed,
            time_probes=time_probes_parsed,
            cf_probes=cf_probes_parsed,
            boundary_layer_probes=bl_probes_parsed,
            large_data_probes=large_data_probes_parsed
        )

        return cls(
            identity=identity, geometry=geometry, domain_sizing=domain_sizing,
            mesh_control=mesh_control, flow_physics=flow_physics, boundary_layer_setup=boundary_layer_setup,
            boundary_conditions=boundary_conditions, simulation_control=simulation_control, io_and_probes=io_and_probes
        )