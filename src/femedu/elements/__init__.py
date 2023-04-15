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
    "DrawElement",
    "Faces",
    "Face2D",
    "Face3D",
)

