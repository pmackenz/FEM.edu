import numpy as np

from .SolidTransformation import *

class Solid2dTransformation(SolidTransformation):
    r"""
    A Transformation including only displacement dofs (no rotational dofs)
    """

    def __init__(self, dir1=None, dir2=None, axis=None):
        super().__init__(dir1, dir2, axis)
