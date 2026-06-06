"""
Core package for the CharLES Pre-Processing Pipeline.
Exposes the primary classes to establish a clean, professional API.
"""

from .config_parser import SimulationState
from .physics_engine import (
    PhysicsEngine,
    BoundaryLayerStrategy,
    EckertAnalyticalStrategy,
    SimilaritySolutionStrategy,
    FixedLengthStrategy,
    LocalDNSStrategy
)
from .mesh_planner import MeshPlanner
from .probe_generator import ProbeGenerator
from .template_writer import TemplateWriter

# Defines exactly what is exported when someone runs 'from core import *'
__all__ = [
    "SimulationState",
    "PhysicsEngine",
    "BoundaryLayerStrategy",
    "EckertAnalyticalStrategy",
    "SimilaritySolutionStrategy",
    "FixedLengthStrategy",
    "LocalDNSStrategy",
    "MeshPlanner",
    "ProbeGenerator",
    "TemplateWriter"
]