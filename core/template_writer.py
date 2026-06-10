"""
Module: template_writer
Description: 
    Responsible for the secure injection of calculated physical and structural 
    parameters into the unified solver template (.in) and batch scripts (.sh).
    
    Dynamically maps Space, Time, Mean Flow (Cf/BL), and Large Data (SPOD) 
    probes. It parses the unified template twice (Transient and Steady) to inject
    phase-specific simulation times, restart logic, and extraction intervals.
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
        self.state = state
        self.physics = physics_state
        self.mesh = mesh_state
        self.template_dir = template_dir
        self.run_name = self.state.identity.run_name
        self.output_dir = os.path.join("output_simulations", self.run_name)
        
        self.queues_path = os.path.join(self.template_dir, "queues.yaml")
        self.queue_data = self._load_queues()

    def _load_queues(self) -> Dict[str, Any]:
        """Loads HPC hardware configurations from the queues YAML file."""
        if not os.path.exists(self.queues_path):
            raise FileNotFoundError(f"[ERROR] Queue configuration file not found at {self.queues_path}")
        with open(self.queues_path, 'r') as file:
            return yaml.safe_load(file)

    def _build_replacement_dict(self, step: str) -> Dict[str, Any]:
        """
        Constructs a key-value mapping between the DEFINE tokens expected in the 
        .in templates and the physical/numerical data.
        """
        phys = self.state.flow_physics
        geom = self.state.geometry
        ctrl = self.state.simulation_control
        io = self.state.io_and_probes
        
        # Thermal boundary condition string formatting
        raw_wall_bc = self.state.boundary_conditions.wall_bc.upper()
        wall_t = getattr(self.state.boundary_conditions, 'wall_T', 1.0)
        formatted_wall_bc = f"WALL_ISOTHERMAL T_WALL {wall_t}" if raw_wall_bc == "ISOTHERMAL" else ("WALL_ADIABATIC" if raw_wall_bc == "ADIABATIC" else raw_wall_bc)
            
        # Total Domain Length & Adimensional Simulation Time 
        total_domain_length = self.physics.pregap_slip_length + self.physics.pregap_noslip_length + geom.structure_length + self.physics.postgap_length
        ftt_duration = total_domain_length / phys.inflow_u
        
        # Phase-Specific Logic (Transient vs Steady)
        if step == "transient":
            simtime = ctrl.transient_simtime_ftt * ftt_duration
            restart_cmd = "RESTART ./helping_files/restart.mles"
            p_inf_calc = f"({phys.inflow_u}*{phys.inflow_u}*{phys.inflow_rho}*1.0/({phys.gamma}*{phys.mach_number}*{phys.mach_number}))"
            # Inicia fluido em repouso para evitar impulsive start
            init_cmd = f"INIT_RUP {phys.inflow_rho} 0 0 0 {p_inf_calc}"
        else:
            simtime = ctrl.steady_simtime_ftt * ftt_duration
            restart_cmd = "RESTART ./helping_files/restart.mles sles"
            init_cmd = "# Flow already developed. Initialization bypassed in steady state."

        return {
            "TOTAL_DOMAIN_LENGTH": total_domain_length,
            "P_INF": f"({phys.inflow_u}*{phys.inflow_u}*{phys.inflow_rho}*1.0/({phys.gamma}*{phys.mach_number}*{phys.mach_number}))",
            "MU_INF": f"({phys.inflow_rho}*{phys.inflow_u}*1.0/{phys.target_reynolds})",
            "RESTART_CMD": restart_cmd,
            "PR_LAM_VAL": phys.prandtl,
            "MU_POWER_LAW_VAL": phys.mu_power_law,
            "RHO_INF": phys.inflow_rho,
            "T_INF": phys.inflow_T,
            "U_INF": phys.inflow_u,
            "G": phys.gamma,
            "CFL_VAL": ctrl.cfl,
            "SIMTIME_VAL": round(simtime, 4),
            "CHECK_INTERVAL_VAL": io.check_interval_steps,
            "SGS_MODEL_VAL": self.state.identity.sgs_model.upper(),
            "INIT_CMD": init_cmd,
            "WALL_BC": formatted_wall_bc,
            
            # --- NOVAS VARIÁVEIS DE METADADOS (YAML) ---
            "GEOM_STRUCTURE_TYPE": geom.structure_type,
            "GEOM_SPAN_Z": geom.span_z,
            "GEOM_DOMAIN_Z": geom.domain_z,
            "MESH_TARGET_Y_PLUS": self.state.mesh_control.target_y_plus,
            "PHYSICS_TARGET_REYNOLDS": phys.target_reynolds,
            "PHYSICS_MACH_NUMBER": phys.mach_number,
            "WALL_BC_TYPE": self.state.boundary_conditions.wall_bc,
            "WALL_T": wall_t,
            "Z_BOUNDARIES": self.state.boundary_conditions.z_boundaries,
            # -------------------------------------------
            
            # Legacy Stitch Replacements
            "HCP_DELTA_VAL": self.mesh.hcp_delta,
            "MAX_REFINEMENT_LEVEL": self.mesh.max_refinement_level,
            "NLAYERS_VAL": self.mesh.nlayers,
            "NSMOOTH_VAL": self.mesh.nsmooth,
            "REFINEMENT_HEIGHT_ABOVE_PLATE": self.mesh.refinement_height,
            "L": geom.structure_length,
            "D": geom.structure_depth,
            "HALF_W": self.mesh.half_w,
            "L1": self.physics.pregap_slip_length,
            "L2": self.physics.pregap_noslip_length,
            "L3": self.physics.postgap_length,
            "H": self.physics.domain_height,
        }
    
    
    def _setup_directories(self):
        """Builds the strict directory hierarchy required for simulation results."""
        dirs_to_create = [
            "helping_files", "logs", "probe_coordinates",
            "results/images", "results/mesh_images",
            "results/space_probes", "results/time_probes",
            "results/cf_probes", "results/boundary_layer_probes", 
            "results/large_data_probes", "results/full_mesh"
        ]
        for directory in dirs_to_create:
            os.makedirs(os.path.join(self.output_dir, directory), exist_ok=True)

    def _generate_probe_lines(self, step: str) -> str:
        """
        Generates the standardized POINTCLOUD_PROBE command blocks for all categories.
        Reads custom variables and independent writing intervals directly from the YAML.
        """
        lines = []
        io = self.state.io_and_probes
        
        def append_probe(registry: dict, category_name: str, prefix: str, is_large_data: bool = False):
            for name, config in registry.items():
                
                # ========================================================
                # CORREÇÃO: Criar as subpastas específicas de cada probe!
                # ========================================================
                probe_dir = os.path.join(self.output_dir, "results", category_name, name)
                os.makedirs(probe_dir, exist_ok=True)
                
                csv_path = f"probe_coordinates/{prefix}_{name}.csv"
                out_path = f"./results/{category_name}/{name}/data"
                
                # Assign interval based on probe type and simulation phase (Transient/Steady)
                if is_large_data:
                    interval = config.transient_write_interval if step == "transient" else config.steady_write_interval
                else:   
                    interval = config.write_interval
                
                # Format variables explicitly requested by the user
                vars_str = " ".join(config.variables)
                lines.append(f"POINTCLOUD_PROBE NAME={out_path} INTERVAL {interval} GEOM=FILE {csv_path} VARS={vars_str}")

        # Dispatch probe injection blocks
        append_probe(io.space_probes, "space_probes", "SpaceProbe")
        append_probe(io.time_probes, "time_probes", "TimeProbe")
        append_probe(io.cf_probes, "cf_probes", "CfProbe")
        append_probe(io.boundary_layer_probes, "boundary_layer_probes", "BLProbe")
        append_probe(io.large_data_probes, "large_data_probes", "LargeData", is_large_data=True)
        
        # Inject Full Mesh Dump if enabled
        if getattr(io, 'save_full_mesh', False):
            lines.append(f"WRITE_DATA NAME=./results/full_mesh/domain INTERVAL {io.full_mesh_write_interval}")
            
        return "\n".join(lines)
        
    def _generate_image_lines(self) -> str:
        """
        Generates native WRITE_IMAGE blocks based on user-defined configurations.
        Defaults the slice to the domain mid-plane (Z_PLANE_FRAC 0.5) to capture 2D physics.
        """
        lines = []
        for name, config in self.state.io_and_probes.image_outputs.items():
            out_path = f"./results/images/{name}"
            lines.append(f"WRITE_IMAGE NAME={out_path} INTERVAL {config.write_interval} GEOM=Z_PLANE_FRAC 0.5 VAR={config.variable} SIZE 4000 4000 UP=0 1 0")
        return "\n".join(lines)
    
    def _process_in_file(self, template_path: str, output_name: str, step: str = "generic"):
        """Reads a CharLES .in template, safely injects variables, and expands dynamically generated blocks."""
        output_path = os.path.join(self.output_dir, output_name)
        
        # Fetch replacement mappings adapted to the current step (Transient/Steady)
        replacements = self._build_replacement_dict(step)
        
        with open(template_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                # 1. Expand dynamic probe block
                if "{{PROBES_INJECTION}}" in line:
                    outfile.write(self._generate_probe_lines(step) + "\n")
                    continue
                
                # 2. Expand dynamic image visualization block
                if "{{IMAGES_INJECTION}}" in line:
                    outfile.write(self._generate_image_lines() + "\n")
                    continue
                
                # 3. Replace DEFINE macros
                if line.startswith("DEFINE "):
                    parts = line.split("=")
                    if len(parts) == 2:
                        var_name = parts[0].replace("DEFINE", "").strip()
                        if var_name in replacements:
                            val = replacements[var_name]
                            outfile.write(f"DEFINE {var_name} = {val}\n")
                            continue
                
                # 4. Replace dynamically injected $(CMD) macros
                for key, val in replacements.items():
                    if f"$({key})" in line:
                        line = line.replace(f"$({key})", str(val))
                        
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
                if line.startswith("#PBS -N"):
                    outfile.write(f"#PBS -N {job_name}\n")
                    continue
                if line.startswith("#PBS -q"):
                    outfile.write(f"#PBS -q {q_name}\n")
                    continue
                
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
        
        os.chmod(output_path, 0o755)
        print(f"[INFO] Generated bash script: {output_name} (Allocated Queue: {q_name})")

    def execute(self):
        """Main orchestrator executing the translation of templates into run-ready scripts."""
        print("--- Generating Simulation Files ---")
        self._setup_directories()
        
        # The writer pulls the unified template to build both the Transient and Steady files
        unified_solver_template = os.path.join(self.template_dir, "template-charles.in")
        stitch_template = os.path.join(self.template_dir, "template-stitch.in")
        surfer_template = os.path.join(self.template_dir, "template-surfer.in")

        if os.path.exists(unified_solver_template):
            self._process_in_file(unified_solver_template, "transient_charles_ig.in", step="transient")
            self._process_in_file(unified_solver_template, "steady_charles_ig.in", step="steady")
        else:
            print(f"[ERROR] Unified solver template not found at {unified_solver_template}")

        if os.path.exists(stitch_template):
            self._process_in_file(stitch_template, "stitch.in", step="generic")
            
        if os.path.exists(surfer_template):
            self._process_in_file(surfer_template, "surfer.in", step="generic")

        # Process PBS submission shell scripts (Standardized names without '_2')
        sh_files = {
            "template-run-surfer.sh": ("run-surfer.sh", "surfer"),
            "template-run-stitch.sh": ("run-stitch.sh", "stitch"),
            "template-run-transient-charles.sh": ("run-transient-charles.sh", "transient"),
            "template-run-steady-charles.sh": ("run-steady-charles.sh", "steady")
        }

        for tpl_file, (out_file, step_name) in sh_files.items():
            tpl_path = os.path.join(self.template_dir, tpl_file)
            if os.path.exists(tpl_path):
                self._process_sh_file(tpl_path, out_file, step_name)
            else:
                print(f"[WARNING] Template {tpl_file} not found in {self.template_dir}")

        print(f"--- Simulation Setup Complete: {self.run_name} ---")