"""
Module: probe_generator
Description: Generates highly customized 3D spatial and temporal point clouds 
             using vectorized NumPy operations. Outputs strict, CharLES-compatible 
             .csv coordinate files (no headers, space-separated). Supports standard
             diagnostics alongside advanced Mean Flow (Cf, BL) and Large Data (SPOD) blocks.
"""

import os
import numpy as np
import pandas as pd
from typing import Any

class ProbeGenerator:
    """
    Translates user-defined bounds and spacing laws into raw X, Y, Z coordinate matrices
    for CharLES POINTCLOUD_PROBE arrays.
    """
    def __init__(self, state: Any, physics_state: Any):
        """
        Args:
            state: SimulationState (parsed securely from YAML)
            physics_state: PhysicsState (computed fluid properties, e.g., delta99)
        """
        self.state = state
        self.physics = physics_state
        self.output_dir = os.path.join(f"output_simulations/{state.identity.run_name}/probe_coordinates")
        os.makedirs(self.output_dir, exist_ok=True)

    def _generate_1d_array(self, bounds: list, config: Any) -> np.ndarray:
        """
        Core mathematical engine generating 1D point distributions based on user spacing rules.
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
            normalized_space = np.cumsum(np.linspace(0, 1, pts))
            normalized_space = (normalized_space - normalized_space.min()) / (normalized_space.max() - normalized_space.min())
            return start + normalized_space * (end - start)
            
        elif config.type == "custom":
            # Injects the explicit coordinate array provided by the user in the YAML
            return np.array(config.custom_vector)
            
        elif config.type in ["mesh_like", "mesh_multiple"]:
            # Leverages the physics_engine's estimated boundary layer thickness (delta99) 
            # to dynamically approximate the required physical clustering near solid surfaces.
            multiplier = config.value if config.value else 1.0
            estimated_pts = int(20 * multiplier * self.physics.estimated_delta99)
            estimated_pts = max(5, estimated_pts) # Ensure at least 5 resolution points exist
            space = np.geomspace(1e-4, 1.0, estimated_pts)
            normalized_space = (space - space.min()) / (space.max() - space.min())
            return start + normalized_space * (end - start)
            
        else:
            raise ValueError(f"Unknown spacing type: {config.type}")

    def _save_pointcloud(self, x_arr: np.ndarray, y_arr: np.ndarray, z_arr: np.ndarray, filename: str, probe_type: str):
        """
        Standardized 3D Cartesian product meshgrid generator and CSV exporter.
        Enforces strict formatting (space-separated, no headers) required by CharLES.
        """
        # Vectorized coordinate expansion natively in C via NumPy
        X, Y, Z = np.meshgrid(x_arr, y_arr, z_arr, indexing='ij')
        coordinates = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
        
        filepath = os.path.join(self.output_dir, filename)
        
        # Pandas export strictly forbidding headers and commas to prevent CharLES C++ parser crashes
        df = pd.DataFrame(coordinates, columns=['x', 'y', 'z'])
        df.to_csv(filepath, index=False, header=False, sep=' ')
        
        print(f"[INFO] Generated {probe_type}: {filename} -> {len(coordinates):,} nodes")

    def generate_space_probes(self):
        """Generates low-frequency statistical bounding coordinate files (Volumetric Snapshots)."""
        for region_name, config in self.state.io_and_probes.space_probes.items():
            x_arr = self._generate_1d_array(config.x_bounds, config.x_spacing)
            y_arr = self._generate_1d_array(config.y_bounds, config.y_spacing)
            z_arr = self._generate_1d_array(config.z_bounds, config.z_spacing)
            self._save_pointcloud(x_arr, y_arr, z_arr, f"SpaceProbe_{region_name}.csv", "Space Probe")

    def generate_time_probes(self):
        """Generates high-frequency spectral arrays at explicit Z-Planes."""
        for region_name, config in self.state.io_and_probes.time_probes.items():
            x_arr = self._generate_1d_array(config.x_bounds, config.x_spacing)
            y_arr = self._generate_1d_array(config.y_bounds, config.y_spacing)
            z_arr = np.array(config.z_planes)
            self._save_pointcloud(x_arr, y_arr, z_arr, f"TimeProbe_{region_name}.csv", "Time Probe")

    def generate_cf_probes(self):
        """Generates heavily compressed wall-normal matrices strictly for calculating du/dy viscous stresses."""
        for region_name, config in self.state.io_and_probes.cf_probes.items():
            x_arr = self._generate_1d_array(config.x_bounds, config.x_spacing)
            y_arr = self._generate_1d_array(config.y_bounds, config.y_spacing)
            z_arr = np.array(config.z_planes)
            self._save_pointcloud(x_arr, y_arr, z_arr, f"CfProbe_{region_name}.csv", "Mean Flow (Cf) Probe")

    def generate_boundary_layer_probes(self):
        """Generates vertical probe lines tailored for continuous integration of boundary layer thicknesses."""
        for region_name, config in self.state.io_and_probes.boundary_layer_probes.items():
            x_arr = self._generate_1d_array(config.x_bounds, config.x_spacing)
            y_arr = self._generate_1d_array(config.y_bounds, config.y_spacing)
            z_arr = np.array(config.z_planes)
            self._save_pointcloud(x_arr, y_arr, z_arr, f"BLProbe_{region_name}.csv", "Mean Flow (BL) Probe")

    def generate_large_data_probes(self):
        """Generates massive 3D volumetric matrices deployed strictly for SPOD/POD and Resolvent Analysis."""
        for region_name, config in self.state.io_and_probes.large_data_probes.items():
            x_arr = self._generate_1d_array(config.x_bounds, config.x_spacing)
            y_arr = self._generate_1d_array(config.y_bounds, config.y_spacing)
            z_arr = self._generate_1d_array(config.z_bounds, config.z_spacing)
            self._save_pointcloud(x_arr, y_arr, z_arr, f"LargeData_{region_name}.csv", "Large Data (SPOD) Probe")

    def execute(self):
        """Main execution orchestrator for the mathematical generator."""
        print("--- Generating Probes ---")
        self.generate_space_probes()
        self.generate_time_probes()
        self.generate_cf_probes()
        self.generate_boundary_layer_probes()
        self.generate_large_data_probes()
        print("--- Probe Generation Complete ---")