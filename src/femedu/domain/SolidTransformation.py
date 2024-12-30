import numpy as np

from .Transformation import *

class SolidTransformation(Transformation):
    r"""
    A Transformation including only displacement dofs (no rotational dofs)
    """

    def __init__(self, dir1=None, dir2=None, axis=None):
        super().__init__(dir1, dir2, axis)
