"""
Main entry point for the CharLES Pre-Processing Pipeline.
"""
import os
import sys
import argparse

# Clean and professional imports using the new __init__.py API
from core import (
    SimulationState,
    PhysicsEngine,
    MeshPlanner,
    ProbeGenerator,
    TemplateWriter
)

def parse_args():
    """Parses command line arguments."""
    parser = argparse.ArgumentParser(description="CharLES CFD Pre-Processing Pipeline")
    parser.add_argument(
        "--config", 
        type=str, 
        required=True, 
        help="Path to the simulation configuration YAML file."
    )
    return parser.parse_args()

def main():
    print("=======================================================")
    print("      CharLES Pre-Processing Pipeline Initialized      ")
    print("=======================================================")

    args = parse_args()
    yaml_path = args.config

    if not os.path.exists(yaml_path):
        print(f"[ERROR] Configuration file not found at {yaml_path}")
        sys.exit(1)

    try:
        # 1. Parse Configuration
        print("[1/5] Loading and validating configuration...")
        state = SimulationState.from_yaml(yaml_path)
        print(f"      -> Run Name: {state.identity.run_name} | Solver: {state.identity.solver_type}")

        # 2. Compute Fluid Physics and Domain Sizing
        print("[2/5] Calculating fluid properties and domain sizing...")
        physics_engine = PhysicsEngine(state)
        physics_state = physics_engine.calculate_derived_properties()

        # 3. Calculate Mesh Resolution for Stitch
        print("[3/5] Planning Voronoi mesh parameters...")
        mesh_planner = MeshPlanner(state, physics_state)
        mesh_state = mesh_planner.calculate_mesh_parameters()

        # 4. Generate Spatial and Temporal Probes
        print("[4/5] Generating spatial and temporal probe clouds...")
        probe_gen = ProbeGenerator(state, physics_state)
        probe_gen.execute()

        # 5. Write Templates and Finalize Directory
        print("[5/5] Injecting variables into CharLES templates...")
        writer = TemplateWriter(state, physics_state, mesh_state, template_dir="templates")
        writer.execute()

        print("=======================================================")
        print(f"[SUCCESS] Simulation setup complete!")
        print(f"Output directory: output_simulations/{state.identity.run_name}/")
        print("Ready to transfer to Zeus.")
        print("=======================================================")

    except Exception as e:
        print(f"\n[FATAL ERROR] Pipeline failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()