"""
Module: template_writer
Description: 
    Responsável pela injeção segura de parâmetros calculados nos templates (.in e .sh).
    Evita o uso de bibliotecas de template nativas para prevenir conflitos de sintaxe
    com variáveis de ambiente do cluster (PBS) e variáveis internas do CharLES.
    
    Gerencia a hierarquia de diretórios de saída e a alocação dinâmica de hardware
    baseada nas filas configuradas em 'templates/queues.yaml'.
"""

import os
import re
import yaml
from typing import Any, Dict

class TemplateWriter:
    
    def __init__(self, state: Any, physics_state: Any, mesh_state: Any, template_dir: str = "templates"):
        self.state = state
        self.physics = physics_state
        self.mesh = mesh_state
        self.template_dir = template_dir
        self.run_name = self.state.identity.run_name
        self.output_dir = os.path.join("output_simulations", self.run_name)
        
        # Carrega a definição de filas/hardware do cluster Zeus
        self.queues_path = os.path.join(self.template_dir, "queues.yaml")
        self.queue_data = self._load_queues()
        
        # Cria o dicionário de mapeamento para substituição de tokens
        self.replacements = self._build_replacement_dict()

    def _load_queues(self) -> Dict[str, Any]:
        """Carrega configurações de hardware das filas HPC."""
        if not os.path.exists(self.queues_path):
            raise FileNotFoundError(f"[ERROR] Arquivo de filas não encontrado em {self.queues_path}")
        with open(self.queues_path, 'r') as file:
            return yaml.safe_load(file)

    def _build_replacement_dict(self) -> Dict[str, Any]:
        """
        Constrói um mapeamento (chave -> valor) entre os tokens esperados nos 
        arquivos .in e os dados calculados pelo framework.
        """
        phys = self.state.flow_physics
        geom = self.state.geometry
        ctrl = self.state.simulation_control
        io = self.state.io_and_probes
        
        # Formata a string de condição de contorno da parede para a sintaxe estrita do CharLES
        raw_wall_bc = self.state.boundary_conditions.wall_bc.upper()
        wall_t = getattr(self.state.boundary_conditions, 'wall_T', 1.0)
        
        if raw_wall_bc == "ISOTHERMAL":
            formatted_wall_bc = f"WALL_ISOTHERMAL T_WALL {wall_t}"
        elif raw_wall_bc == "ADIABATIC":
            formatted_wall_bc = "WALL_ADIABATIC"
        else:
            formatted_wall_bc = raw_wall_bc
            
        # Cálculo de Comprimento Total e Conversão de FTT para Tempo Adimensional
        total_domain_length = self.physics.pregap_slip_length + self.physics.pregap_noslip_length + geom.structure_length + self.physics.postgap_length
        u_inf = phys.inflow_u
        
        # 1 FTT = Tempo necessário para o escoamento cruzar todo o domínio
        ftt_duration = total_domain_length / u_inf
        
        # Conversão dos tempos do YAML (em FTT) para tempo de simulação do CharLES
        actual_transient_time = ctrl.transient_simtime_ftt * ftt_duration
        actual_steady_time = ctrl.steady_simtime_ftt * ftt_duration
        
        return {
            # --- Parâmetros Geométricos ---
            "L": geom.structure_length,
            "D": geom.structure_depth,
            "HALF_W": self.mesh.half_w,
            "L1": self.physics.pregap_slip_length,
            "L2": self.physics.pregap_noslip_length,
            "L3": self.physics.postgap_length,
            "H": self.physics.domain_height,
            "TOTAL_DOMAIN_LENGTH": total_domain_length,
            
            # --- Parâmetros de Malha (Stitch) ---
            "HCP_DELTA_VAL": self.mesh.hcp_delta,
            "MAX_REFINEMENT_LEVEL": self.mesh.max_refinement_level,
            "NLAYERS_VAL": self.mesh.nlayers,
            "NSMOOTH_VAL": self.mesh.nsmooth,
            "REFINEMENT_HEIGHT_ABOVE_PLATE": self.mesh.refinement_height,
            
            # --- Parâmetros de Física de Escoamento ---
            "RE_DELTA_STAR": phys.target_reynolds,
            "MACH": phys.mach_number,
            "U_INF": phys.inflow_u,
            "RHO_INF": phys.inflow_rho,
            "T_INF": phys.inflow_T,
            "G": phys.gamma,
            "PR_LAM_VAL": phys.prandtl,
            "MU_POWER_LAW_VAL": phys.mu_power_law,
            
            # --- Condições de Contorno e Solver ---
            "WALL_BC": formatted_wall_bc,
            "SGS_MODEL_VAL": "NONE" if self.state.identity.solver_type.upper() == "DNS" else "VREMAN",
            
            # --- Controle de Simulação (Agora convertidos de FTT para tempo real) ---
            "TRANSIENT_SIMTIME_VAL": round(actual_transient_time, 4),
            "STEADY_SIMTIME_VAL": round(actual_steady_time, 4),
            "CFL_VAL": ctrl.cfl,
            
            # --- Intervalos de Escrita (IO) ---
            "CHECK_INTERVAL_VAL": io.check_interval_steps,
            "IMAGE_WRITE_INTERVAL": io.image_interval_steps,
            "SPACE_PROBES_WRITE_INTERVAL": io.space_probes_write_interval,
            "TIME_PROBES_WRITE_INTERVAL": io.time_probes_write_interval
        }
    def _setup_directories(self):
        """Cria a estrutura de pastas necessária para os resultados da simulação."""
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "helping_files"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "logs"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "results", "images"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "results", "mesh_images"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "probe_coordinates"), exist_ok=True)

    def _generate_probe_lines(self) -> str:
        """Gera o bloco de texto padronizado usando POINTCLOUD_PROBE para evitar memory leak no solver."""
        lines = []
        io = self.state.io_and_probes
        
        # Variáveis exclusivas para as Time Probes (apenas instantâneas, conforme solicitado)
        time_vars = "comp(u,0) comp(u,1) comp(u,2) p rho T"
        
        # Variáveis para as Space Probes (instantâneas + médias temporais)
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
            
            # Alterado obrigatoriamente para POINTCLOUD_PROBE
            lines.append(f"POINTCLOUD_PROBE NAME=./results/time_probes/{name}/data INTERVAL {interval} GEOM=FILE {csv_path} VARS={time_vars}")
            
        return "\n".join(lines)
    
    def _process_in_file(self, template_path: str, output_name: str):
        """Lê um template .in, injeta variáveis 'DEFINE' e expande blocos dinâmicos."""
        output_path = os.path.join(self.output_dir, output_name)
        with open(template_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                # 1. Expandir bloco dinâmico de sondas
                if "{{PROBES_INJECTION}}" in line:
                    outfile.write(self._generate_probe_lines() + "\n")
                    continue
                
                # 2. Injeção de variáveis DEFINE
                if line.startswith("DEFINE "):
                    parts = line.split("=")
                    if len(parts) == 2:
                        var_name = parts[0].replace("DEFINE", "").strip()
                        if var_name in self.replacements:
                            val = self.replacements[var_name]
                            outfile.write(f"DEFINE {var_name} = {val}\n")
                            continue
                outfile.write(line)
        print(f"[INFO] Gerado arquivo de entrada: {output_name}")

    def _process_sh_file(self, template_path: str, output_name: str, step_name: str):
        """Processa scripts PBS, injetando o nome do job, a fila correta e recursos."""
        output_path = os.path.join(self.output_dir, output_name)
        job_name = f"{self.run_name}_{step_name}"
        
        q_name = self.state.simulation_control.pbs_queue
        if q_name not in self.queue_data:
            raise ValueError(f"[ERROR] Fila '{q_name}' não definida em queues.yaml!")
            
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
        print(f"[INFO] Gerado script bash: {output_name} (Fila: {q_name})")

    def execute(self):
        """Método principal para executar o processamento de todos os templates."""
        print("--- Gerando Arquivos de Simulação ---")
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
                print(f"[WARNING] Template {tpl_file} não encontrado")

        for tpl_file, (out_file, step_name) in sh_files.items():
            tpl_path = os.path.join(self.template_dir, tpl_file)
            if os.path.exists(tpl_path):
                self._process_sh_file(tpl_path, out_file, step_name)
            else:
                print(f"[WARNING] Template {tpl_file} não encontrado")

        print(f"--- Setup de Simulação Finalizado: {self.run_name} ---")