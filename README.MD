# CharLES Hypersonic Cavity Simulation Pipeline

This repository contains a robust, automated Python-based framework designed to prepare, execute, and post-process high-fidelity 3D Direct Numerical Simulations (DNS) and Large Eddy Simulations (LES) for hypersonic boundary layers and cavity flows using the Cascade Technologies CharLES solver.

## 📌 Core Architecture & Features

* **Unified Solver Template:** Eliminates redundancy by utilizing a single, dynamically injected `template-charles.in` for both Transient and Steady states. The Python pipeline automatically manages initialization commands (avoiding impulsive start thermal shocks), restart flags, and extraction configurations.
* **Dynamic Sub-Grid Scale (SGS) Logic:** Strictly enforces fluid dynamics principles. If `DNS` is selected, SGS models are rigorously disabled (`NONE`). If `LES` is selected, it automatically applies the user-defined model (defaulting to `VREMAN` if omitted).
* **Advanced Probing Architecture:** Bypasses native memory leaks by exclusively using vectorized `POINTCLOUD_PROBE` topologies. Introduces specific mathematical domains for:
  * *Mean Flow Analysis:* Skin Friction (`Cf`) arrays capturing wall-normal gradients (du/dy), and Boundary Layer (`BL`) profiles for integral thicknesses.
  * *Modal Decomposition:* Massive volumetric blocks (`Large Data`) optimized for SPOD/POD and Resolvent Analysis, featuring independent writing frequencies for Transient (low-frequency) and Steady (high-frequency) phases.
* **Native Volumetric & Image Diagnostics:** Granular control directly from the YAML to export 2D `.png` slices (e.g., density, dilatation, velocity) and full-domain mesh outputs without altering solver scripts.
* **Intelligent Mesh Sizing:** Translates theoretical boundary layer physics (y+ targets) into absolute Voronoi Hexagonal Close Packed (HCP) constraints. It dynamically detects Quasi-2D vs. 3D domains based on `domain_z`, safely capping memory allocation to prevent Out-Of-Memory (OOM) failures.

## 💻 Environment & Requirements

### HPC Cluster Requirements
* **Solver:** CharLES `fidelityCharles_2024.2`
* **MPI Implementation:** Intel MPI (`impi2021.10.0` or equivalent)
* **Compute:** High-memory nodes are required for Voronoi octree generation (RAM usage scales heavily with the `MAX_REFINEMENT_LEVEL`).

### Data Engineering & Post-Processing Requirements
* **Languages:** Python 3.x and/or MATLAB.
* **Compute Restrictions:** All post-processing scripts must execute exclusively on the CPU. GPU acceleration and parallelization/multiprocessing features must be strictly disabled.
* **File Management:** All automated data writing operations must enforce a strict 10-step Safe-Save Protocol to guarantee absolute protection against data corruption.

## 📂 Project Structure

```Ensure the working directory strictly follows this hierarchy:

├── core/
│   ├── config_parser.py          # YAML architecture and strict logic typing
│   ├── physics_engine.py         # Boundary layer and hypersonic scaling
│   ├── mesh_planner.py           # Physics-to-discretization logic
│   ├── probe_generator.py        # Vectorized coordinate generator
│   └── template_writer.py        # YAML-to-CharLES syntax injection
├── templates/
│   ├── queues.yaml               # HPC cluster job scheduling parameters
│   ├── template-surfer.in        # Stereolithography surface generator
│   ├── template-stitch.in        # Voronoi volume generator
│   └── template-charles.in       # Unified fluid solver template
├── configs/
│   └── simulation_config.yaml    # Master control file (User Input)
└── run_preprocess.py             # Main pipeline initialization script
```

## ⚙️ Execution Workflow

### 1. Initialization
Define your fluid properties, structural geometry, extraction intervals, and target constraints in `configs/simulation_config.yaml`. Initialize the pipeline:

> python run_preprocess.py configs/simulation_config.yaml

This creates a dedicated `output_simulations/<run_name>/` directory natively populated with subdirectories, generated coordinate arrays, and cluster submission scripts perfectly mapped to the available HPC queues.

### 2. Mesh Construction
Navigate to your specific simulation directory and execute the meshing sequence:

> cd output_simulations/<run_name>
> sbatch run-surfer.sh  

Wait for Surfer to finish, then build the Voronoi octree volume:

> sbatch run-stitch.sh  

### 3. Flow Solver Execution
Once `stitch` successfully outputs the `restart.mles` file, submit the transient fluid solver to initialize the flow and flush the domain:

> sbatch run-transient-charles.sh

After the transient phase reaches a stable development, launch the steady-state data collection:

> sbatch run-steady-charles.sh

### 4. Post-Processing & Consolidation
Upon completion, the solver will populate the `results/` directories (`time_probes`, `space_probes`, `cf_probes`, `boundary_layer_probes`, `large_data_probes`, and `images`). 

Run your dedicated Python or MATLAB scripts to sequentially parse these files into continuous, time-resolved datasets for Modal and Mean Flow analyses.

## ⚠️ Troubleshooting & Critical Guidelines

* **Probe File Formatting:** The CharLES parser will fatally crash if coordinate files contain text headers or commas. The framework handles this via strict Pandas export rules; do not manually modify the `.csv` files inside `probe_coordinates/`.
* **MPI Over-partitioning (Signal 9):** If testing small meshes (under 2 million cells) in 2D, reduce the MPI ranks significantly (e.g., `-np 64`). Using 256+ cores on a small mesh creates an enormous MPI communication buffer overhead that will trigger the OS OOM Killer during probe allocation.
* **Transient Initialization:** The framework enforces an initial condition of u=0 to prevent infinite velocity gradients at the wall during high-Mach setups. Allow sufficient Flow-Through Times (FTT) for the freestream to fully establish the boundary layer before switching to steady state.
