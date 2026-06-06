# CharLES Hypersonic Cavity Simulation Pipeline

This repository contains a robust, automated Python-based framework designed to prepare, execute, and post-process high-fidelity 3D Direct Numerical Simulations (DNS) and Large Eddy Simulations (LES) for hypersonic boundary layers and cavity flows using the Cascade Technologies CharLES solver.

## 📌 Core Architecture & Features

* **Intelligent Mesh Sizing:** Automatically calculates Voronoi Hexagonal Close Packed (HCP) base parameters and octree refinement levels. It translates theoretical boundary layer physics and target wall spacings into structural limits, preventing Out-Of-Memory (OOM) crashes during generation.
* **Dynamic Template Injection:** Reads YAML configuration files and securely maps parameters to strict CharLES `.in` syntax. It handles complex implementations, such as isothermal wall boundaries, without breaking solver parsing rules.
* **Stable Volumetric Probing:** Implements `POINTCLOUD_PROBE` linked to KD-Tree structures for all spatial and temporal data extraction. This fundamentally bypasses the severe memory leaks and OOM crashes natively triggered when using standard `PROBE` commands with heavy external coordinate files.
* **Automated Time Conversions:** Translates user-defined Flow-Through Times (FTT) into the strict non-dimensional simulation times required by the solver, dynamically adjusting to changes in the domain length.

## 💻 Environment & Requirements

### HPC Cluster Requirements
* **Solver:** CharLES `fidelityCharles_2024.2`
* **MPI Implementation:** Intel MPI (`impi2021.10.0` or equivalent)
* **Compute:** High-memory nodes are required for Voronoi octree generation (RAM usage scales heavily with the `MAX_REFINEMENT_LEVEL`).

### Data Engineering & Post-Processing Requirements
* **Languages:** Python 3.x and/or MATLAB.
* **Compute Restrictions:** All post-processing and data engineering scripts must execute exclusively on the CPU. GPU acceleration and parallelization/multiprocessing features must be strictly disabled.
* **File Management:** All automated data writing operations must enforce a strict 10-step Safe-Save Protocol (Open -> Create backup -> Close -> Check equality -> etc.) to guarantee absolute protection against data corruption.

## 📂 Project Structure

Ensure the working directory strictly follows this hierarchy:
```text
├── core/
│   ├── mesh_planner.py           # Physics-to-discretization logic
│   └── template_writer.py        # YAML-to-CharLES syntax injection
├── templates/
│   ├── queues.yaml               # HPC cluster job scheduling parameters
│   ├── template-surfer.in        # Stereolithography surface generator
│   ├── template-stitch.in        # Voronoi volume generator
│   ├── template-steady_charles_ig.in
│   └── template-transient_charles_ig.in
├── probe_coordinates/            # Pre-calculated coordinate arrays (.csv)
└── run_preprocess.py             # Main pipeline initialization script
```

## ⚙️ Execution Workflow

### 1. Initialization
Define your fluid properties, structural geometry, and target constraints in the master YAML configuration, then initialize the pipeline:

> python run_preprocess.py

This creates a dedicated `output_simulations/<run_name>/` directory populated with dynamically generated input files and cluster submission scripts perfectly mapped to the available HPC queues.

### 2. Mesh Construction
Navigate to your specific simulation directory and execute the meshing sequence:

> cd output_simulations/<run_name>
> sbatch run-surfer.sh  

Wait for the surface generation to complete successfully, then build the volume:

> sbatch run-stitch.sh  

### 3. Flow Solver Execution
Once `stitch` successfully outputs the `restart.mles` file, submit the main fluid solver:

> sbatch run-transient-charles.sh

### 4. Post-Processing & Consolidation
Upon completion, the solver will populate `results/time_probes/` and `results/space_probes/` with static `.pxyz` coordinate files and multiple iterative `.pcd` data files. 

Run your dedicated Python or MATLAB scripts to sequentially parse these files and consolidate them into continuous, time-resolved datasets. Ensure the data engineering compute restrictions (CPU-only, single-thread, 10-step Safe-Save) are actively enforced during this phase.

## ⚠️ Troubleshooting & Critical Guidelines

* **Probe File Formatting:** The CharLES parser will crash if coordinate files contain text. All `.csv` files inside `probe_coordinates/` must be stripped of headers (e.g., no "x, y, z") and must use a single space as the delimiter, not commas.
* **MPI Over-partitioning (Signal 9):** If testing small meshes (under 2 million cells), reduce the MPI ranks significantly (e.g., `-np 64`). Using 256+ cores on a small mesh creates an enormous MPI communication buffer overhead that will drain the node's RAM and trigger the OS OOM Killer during probe allocation.
* **Isotropic Constraints:** The `stitch` generator creates strictly isotropic Voronoi cells. Modifying resolution limits will affect all dimensions equally. Manage the `HCP_DELTA` cautiously to prevent severe cell count explosions in the spanwise dimension.
