import numpy as np

from .FrameTransformation import *

class BeamTransformation(FrameTransformation):
    r"""
    This is an alias for :py:meth:`FrameTransformation`
    """

    def __init__(self, dir1=None, dir2=None, axis=None):
        super().__init__(dir1, dir2, axis)
