"""
=====================================================
Finite elements based on a variational formulation.
=====================================================

Elements in this group share the following features:
* large or finite deformation kinematics
* based on variational principles
* numeric integration of stiffness and internal force, thus, allowing for the use of inelastic material
"""

__all__ = (
    "Triangle",
    "Quad",
    "Quad8",
    "Quad9",
    "Truss",
    "Frame2D",
)

from .Frame2D import *
from .Quad import *
from .Quad8 import *
from .Quad9 import *
from .Triangle import *
from .Truss import *

