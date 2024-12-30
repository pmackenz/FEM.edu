import numpy as np

from .Transformation import *

class CosseratTransformation(Transformation):
    r"""
    A Transformation including displacement dofs and rotational dofs
    """

    def __init__(self, dir1=None, dir2=None, axis=None):
        super().__init__(dir1, dir2, axis)
