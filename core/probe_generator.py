"""
Module: probe_generator
Description: Generates highly customized 3D spatial and temporal point clouds 
             using vectorized NumPy operations. Outputs strict, CharLES-compatible 
             .csv coordinate files (no headers, space-separated).
"""

import os
import numpy as np
import pandas as pd
from typing import Any

class ProbeGenerator:
    """
    Translates user-defined bounds and spacing laws into raw X, Y, Z coordinate matrices
    for CharLES temporal and spatial POINTCLOUD_PROBE arrays.
    """
    def __init__(self, state: Any, physics_state: Any):
        """
        Args:
            state: SimulationState (parsed securely from YAML)
            physics_state: PhysicsState (computed properties, specifically delta99)
        """
        self.state = state
        self.physics = physics_state
        self.output_dir = os.path.join(f"output_simulations/{state.identity.run_name}/probe_coordinates")
        os.makedirs(self.output_dir, exist_ok=True)

    def _generate_1d_array(self, bounds: list, config: Any) -> np.ndarray:
        """
        Core mathematical engine for generating 1D point distributions based on user spacing rules.
        """
        start, end = bounds[0], bounds[1]
        
        # Fixed axis bypass: If start equals end, collapse the array to a single 2D plane coordinate
        if np.isclose(start, end):
            return np.array([start])

        pts = config.points if config.points is not None else int(config.value)
        
        if config.type == "uniform":
            # Linear, equidistant point distribution
            return np.linspace(start, end, pts)
            
        elif config.type == "logarithmic":
            # Clusters points near the 'start' boundary (e.g., highly resolved near a solid wall)
            # Implements a base-10 logarithmic space mapping, adding an epsilon to prevent log(0)
            space = np.geomspace(1e-5, 1.0, pts) 
            normalized_space = (space - space.min()) / (space.max() - space.min())
            return start + normalized_space * (end - start)
            
        elif config.type == "exponential":
            # Clusters points near the 'end' boundary (e.g., highly resolved in the freestream/shear layer)
            space = np.exp(np.linspace(0, 5, pts))
            normalized_space = (space - space.min()) / (space.max() - space.min())
            return start + normalized_space * (end - start)
            
        elif config.type == "linear":
            # Arithmetic progression (constant sequential increase in distance delta)
            # sum(d + i*d) = total_length
            normalized_space = np.cumsum(np.linspace(0, 1, pts))
            normalized_space = (normalized_space - normalized_space.min()) / (normalized_space.max() - normalized_space.min())
            return start + normalized_space * (end - start)
            
        elif config.type == "custom":
            # Injects the explicit coordinate array provided by the user in the YAML
            return np.array(config.custom_vector)
            
        elif config.type in ["mesh_like", "mesh_multiple"]:
            # Leverages the physics_engine's estimated boundary layer thickness (delta99) 
            # to dynamically approximate the required physical clustering near solid surfaces.
            # Currently defaults to a robust exponential-like clustering methodology.
            multiplier = config.value if config.value else 1.0
            estimated_pts = int(20 * multiplier * self.physics.estimated_delta99)
            estimated_pts = max(5, estimated_pts) # Ensure at least 5 resolution points exist
            space = np.geomspace(1e-4, 1.0, estimated_pts)
            normalized_space = (space - space.min()) / (space.max() - space.min())
            return start + normalized_space * (end - start)
            
        else:
            raise ValueError(f"Unknown spacing type: {config.type}")

    def generate_space_probes(self):
        """Generates bounding coordinate files for spatial probes (volumetric snapshots)."""
        for region_name, config in self.state.io_and_probes.space_probes.items():
            
            x_arr = self._generate_1d_array(config.x_bounds, config.x_spacing)
            y_arr = self._generate_1d_array(config.y_bounds, config.y_spacing)
            z_arr = self._generate_1d_array(config.z_bounds, config.z_spacing)

            # Create the 3D Cartesian product grid
            X, Y, Z = np.meshgrid(x_arr, y_arr, z_arr, indexing='ij')
            
            # Flatten grid into an N x 3 coordinate list
            coordinates = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
            
            # Save to strict CharLES format (No headers, space-delimited)
            filepath = os.path.join(self.output_dir, f"SpaceProbe_{region_name}.csv")
            df = pd.DataFrame(coordinates, columns=['x', 'y', 'z'])
            df.to_csv(filepath, index=False, header=False, sep=' ')
            print(f"[INFO] Generated Space Probe: {region_name} -> {len(coordinates)} points")

    def generate_time_probes(self):
        """Generates explicit Cartesian coordinate files for high-frequency temporal probes."""
        for region_name, config in self.state.io_and_probes.time_probes.items():
            
            x_arr = self._generate_1d_array(config.x_bounds, config.x_spacing)
            y_arr = self._generate_1d_array(config.y_bounds, config.y_spacing)
            z_arr = np.array(config.z_planes)

            # Create the 3D Cartesian product grid mapping across defined Z-planes
            X, Y, Z = np.meshgrid(x_arr, y_arr, z_arr, indexing='ij')
            coordinates = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
            
            # Save to strict CharLES format (No headers, space-delimited)
            filepath = os.path.join(self.output_dir, f"TimeProbe_{region_name}.csv")
            df = pd.DataFrame(coordinates, columns=['x', 'y', 'z'])
            df.to_csv(filepath, index=False, header=False, sep=' ')
            print(f"[INFO] Generated Time Probe: {region_name} -> {len(coordinates)} points")

    def execute(self):
        """Main execution runner for the generation module."""
        print("--- Generating Probes ---")
        self.generate_space_probes()
        self.generate_time_probes()
        print("--- Probe Generation Complete ---")