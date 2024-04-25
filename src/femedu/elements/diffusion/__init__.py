"""
=====================================================
Diffusion elements
=====================================================

Elements in this group share the following features:
* gradients with respect to an undeformed configuration
* flux proportional to gradient of the potential.

"""

__all__ = (
    "Triangle",
    "Triangle6",
)

from .Triangle import *
from .Triangle6 import *

