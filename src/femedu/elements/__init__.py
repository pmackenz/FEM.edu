"""
=====================================================
Finite elements based on a variational formulation.
=====================================================

Elements in this group share the following features:
* small deformation kinematics
* based on variational principles
* numeric integration of stiffness and internal force, thus, allowing for the use of inelastic material
"""

__all__ = (
    "Element",
    "LinearElement",
    "DrawElement",
    "Faces",
    "Face2D",
    "Face3D",
)

from .Faces import *
from .Face2D import *
from .Face3D import *
from .DrawElement import *
from .Element import *
from .LinearElement import *

