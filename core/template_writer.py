"""
Module: template_writer
Description: 
    Responsible for the secure injection of calculated physical and structural 
    parameters into the solver templates (.in and .sh).
    
    Bypasses native template libraries (like Jinja) to prevent syntax conflicts 
    with HPC cluster environment variables (PBS) and CharLES internal commands.
    Manages the output directory hierarchy and dynamically allocates hardware 
    resources based on the cluster queues defined in 'templates/queues.yaml'.
"""

import os
import re
import yaml
from typing import Any, Dict

class TemplateWriter:
    """
    Translates the calculated simulation states into executable CharLES scripts,
    enforcing strict syntax rules and cluster compatibility.
    """
    
    def __init__(self, state: Any, physics_state: Any, mesh_state: Any, template_dir: str = "templates"):
        """
        Args:
            state: SimulationState (parsed securely from YAML)
            physics_state: PhysicsState (computed fluid mechanics properties)
            mesh_state: MeshState (computed Voronoi octree discretization parameters)
            template_dir: Path to the directory containing base .in and .sh templates
        """
        self.state = state
        self.physics = physics_state
        self.mesh = mesh_state
        self.template_dir = template_dir
        self.run_name = self.state.identity.run_name
        self.output_dir = os.path.join("output_simulations", self.run_name)
        
        # Load hardware/queue definitions for the HPC cluster
        self.queues_path = os.path.join(self.template_dir, "queues.yaml")
        self.queue_data = self._load_queues()
        
        # Build the master dictionary for token replacement
        self.replacements = self._build_replacement_dict()

    def _load_queues(self) -> Dict[str, Any]:
        """Loads HPC hardware configurations from the queues YAML file."""
        if not os.path.exists(self.queues_path):
            raise FileNotFoundError(f"[ERROR] Queue configuration file not found at {self.queues_path}")
        with open(self.queues_path, 'r') as file:
            return yaml.safe_load(file)

    def _build_replacement_dict(self) -> Dict[str, Any]:
        """
        Constructs a key-value mapping between the DEFINE tokens expected in the 
        .in templates and the physical/numerical data calculated by the framework.
        """
        phys = self.state.flow_physics
        geom = self.state.geometry
        ctrl = self.state.simulation_control
        io = self.state.io_and_probes
        
        # Format the thermal boundary condition to match strict CharLES syntax
        raw_wall_bc = self.state.boundary_conditions.wall_bc.upper()
        wall_t = getattr(self.state.boundary_conditions, 'wall_T', 1.0)
        
        if raw_wall_bc == "ISOTHERMAL":
            formatted_wall_bc = f"WALL_ISOTHERMAL T_WALL {wall_t}"
        elif raw_wall_bc == "ADIABATIC":
            formatted_wall_bc = "WALL_ADIABATIC"
        else:
            formatted_wall_bc = raw_wall_bc
            
        # Calculate Total Domain Length and convert user FTTs to Adimensional Simulation Time
        total_domain_length = self.physics.pregap_slip_length + self.physics.pregap_noslip_length + geom.structure_length + self.physics.postgap_length
        u_inf = phys.inflow_u
        
        # 1 Flow-Through Time (FTT) = Time required for the freestream to traverse the entire domain
        ftt_duration = total_domain_length / u_inf
        
        # Convert user-defined YAML FTT limits into actual CharLES solver time
        actual_transient_time = ctrl.transient_simtime_ftt * ftt_duration
        actual_steady_time = ctrl.steady_simtime_ftt * ftt_duration
        
        return {
            # --- Geometric Parameters ---
            "L": geom.structure_length,
            "D": geom.structure_depth,
            "HALF_W": self.mesh.half_w,
            "L1": self.physics.pregap_slip_length,
            "L2": self.physics.pregap_noslip_length,
            "L3": self.physics.postgap_length,
            "H": self.physics.domain_height,
            "TOTAL_DOMAIN_LENGTH": total_domain_length,
            
            # --- Volumetric Mesh Parameters (Stitch) ---
            "HCP_DELTA_VAL": self.mesh.hcp_delta,
            "MAX_REFINEMENT_LEVEL": self.mesh.max_refinement_level,
            "NLAYERS_VAL": self.mesh.nlayers,
            "NSMOOTH_VAL": self.mesh.nsmooth,
            "REFINEMENT_HEIGHT_ABOVE_PLATE": self.mesh.refinement_height,
            
            # --- Freestream Fluid Physics ---
            "RE_DELTA_STAR": phys.target_reynolds,
            "MACH": phys.mach_number,
            "U_INF": phys.inflow_u,
            "RHO_INF": phys.inflow_rho,
            "T_INF": phys.inflow_T,
            "G": phys.gamma,
            "PR_LAM_VAL": phys.prandtl,
            "MU_POWER_LAW_VAL": phys.mu_power_law,
            
            # --- Boundary Conditions & Solver Paradigms ---
            "WALL_BC": formatted_wall_bc,
            "SGS_MODEL_VAL": "NONE" if self.state.identity.solver_type.upper() == "DNS" else "VREMAN",
            
            # --- Simulation Control (Converted from FTT to non-dimensional time) ---
            "TRANSIENT_SIMTIME_VAL": round(actual_transient_time, 4),
            "STEADY_SIMTIME_VAL": round(actual_steady_time, 4),
            "CFL_VAL": ctrl.cfl,
            
            # --- I/O and Data Extraction Intervals ---
            "CHECK_INTERVAL_VAL": io.check_interval_steps,
            "IMAGE_WRITE_INTERVAL": io.image_interval_steps,
            "SPACE_PROBES_WRITE_INTERVAL": io.space_probes_write_interval,
            "TIME_PROBES_WRITE_INTERVAL": io.time_probes_write_interval
        }

    def _setup_directories(self):
        """Builds the strict directory hierarchy required for simulation results."""
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "helping_files"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "logs"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "results", "images"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "results", "mesh_images"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "probe_coordinates"), exist_ok=True)

    def _generate_probe_lines(self) -> str:
        """
        Generates the standardized POINTCLOUD_PROBE command blocks.
        Enforces the use of pointclouds to fundamentally prevent native solver memory leaks.
        """
        lines = []
        io = self.state.io_and_probes
        
        # Variables extracted for high-frequency Time Probes (Instantaneous fields only)
        time_vars = "comp(u,0) comp(u,1) comp(u,2) p rho T"
        
        # Variables extracted for low-frequency Space Probes (Instantaneous + Statistical Averages)
        space_vars = "comp(u,0) comp(u,1) comp(u,2) p rho T comp(avg(u),0) comp(avg(u),1) comp(avg(u),2) avg(p) avg(rho) avg(T)"
        
        # --- SPACE PROBES ---
        for name in io.space_probes.keys():
            probe_dir = os.path.join(self.output_dir, "results", "space_probes", name)
            os.makedirs(probe_dir, exist_ok=True)
            
            csv_path = f"probe_coordinates/SpaceProbe_{name}.csv"
            interval = io.space_probes_write_interval
            
            lines.append(f"POINTCLOUD_PROBE NAME=./results/space_probes/{name}/data INTERVAL {interval} GEOM=FILE {csv_path} VARS={space_vars}")
            
        # --- TIME PROBES ---
        for name in io.time_probes.keys():
            probe_dir = os.path.join(self.output_dir, "results", "time_probes", name)
            os.makedirs(probe_dir, exist_ok=True)
            
            csv_path = f"probe_coordinates/TimeProbe_{name}.csv"
            interval = io.time_probes_write_interval
            
            # Mandatorily uses POINTCLOUD_PROBE for hardware stability
            lines.append(f"POINTCLOUD_PROBE NAME=./results/time_probes/{name}/data INTERVAL {interval} GEOM=FILE {csv_path} VARS={time_vars}")
            
        return "\n".join(lines)
    
    def _process_in_file(self, template_path: str, output_name: str):
        """Reads a CharLES .in template, safely injects DEFINE variables, and expands probe blocks."""
        output_path = os.path.join(self.output_dir, output_name)
        with open(template_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                # 1. Expand dynamic probe block
                if "{{PROBES_INJECTION}}" in line:
                    outfile.write(self._generate_probe_lines() + "\n")
                    continue
                
                # 2. Inject mapped DEFINE variables
                if line.startswith("DEFINE "):
                    parts = line.split("=")
                    if len(parts) == 2:
                        var_name = parts[0].replace("DEFINE", "").strip()
                        if var_name in self.replacements:
                            val = self.replacements[var_name]
                            outfile.write(f"DEFINE {var_name} = {val}\n")
                            continue
                outfile.write(line)
        print(f"[INFO] Generated input script: {output_name}")

    def _process_sh_file(self, template_path: str, output_name: str, step_name: str):
        """Processes PBS bash scripts, injecting the correct job name, cluster queue, and core counts."""
        output_path = os.path.join(self.output_dir, output_name)
        job_name = f"{self.run_name}_{step_name}"
        
        q_name = self.state.simulation_control.pbs_queue
        if q_name not in self.queue_data:
            raise ValueError(f"[ERROR] Queue '{q_name}' is not defined in queues.yaml!")
            
        max_cores = self.queue_data[q_name]['ncpus']
        
        with open(template_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                # Inject PBS Job Name
                if line.startswith("#PBS -N"):
                    outfile.write(f"#PBS -N {job_name}\n")
                    continue
                
                # Inject PBS Queue Designation
                if line.startswith("#PBS -q"):
                    outfile.write(f"#PBS -q {q_name}\n")
                    continue
                
                # Dynamically allocate MPI constraints based on yaml definitions
                if step_name in ["steady", "transient"]:
                    if line.startswith("#PBS -l select="):
                        outfile.write(f"#PBS -l select=1:ncpus={max_cores}:mpiprocs={max_cores}\n")
                        continue
                    if "charles_launch.sh" in line:
                        line = re.sub(r"-np \d+", f"-np {max_cores}", line)
                        line = re.sub(r"-perhost \d+", f"-perhost {max_cores}", line)
                        outfile.write(line)
                        continue
                outfile.write(line)
        
        # Grant executable permissions to the generated bash script
        os.chmod(output_path, 0o755)
        print(f"[INFO] Generated bash script: {output_name} (Allocated Queue: {q_name})")

    def execute(self):
        """Main orchestrator executing the translation of all base templates into run-ready scripts."""
        print("--- Generating Simulation Files ---")
        self._setup_directories()
        
        in_files = {
            "template-surfer.in": "surfer.in",
            "template-stitch.in": "stitch.in",
            "template-transient_charles_ig.in": "transient_charles_ig.in",
            "template-steady_charles_ig.in": "steady_charles_ig.in"
        }
        
        sh_files = {
            "template-run-surfer.sh": ("run-surfer.sh", "surfer"),
            "template-run-stitch.sh": ("run-stitch.sh", "stitch"),
            "template-run-transient-charles.sh": ("run-transient-charles.sh", "transient"),
            "template-run-steady-charles.sh": ("run-steady-charles.sh", "steady"),
            "template-move_logs.sh": ("move_logs.sh", "logs")
        }

        for tpl_file, out_file in in_files.items():
            tpl_path = os.path.join(self.template_dir, tpl_file)
            if os.path.exists(tpl_path):
                self._process_in_file(tpl_path, out_file)
            else:
                print(f"[WARNING] Template {tpl_file} not found in {self.template_dir}")

        for tpl_file, (out_file, step_name) in sh_files.items():
            tpl_path = os.path.join(self.template_dir, tpl_file)
            if os.path.exists(tpl_path):
                self._process_sh_file(tpl_path, out_file, step_name)
            else:
                print(f"[WARNING] Template {tpl_file} not found in {self.template_dir}")

        print(f"--- Simulation Setup Complete: {self.run_name} ---")