import numpy as np

from .Frame2dTransformation import *

class Beam2dTransformation(Frame2dTransformation):
    r"""
    This is an alias for :py:meth:`Frame2dTransformation`.
    See documentation for the :py:meth:`Frame2dTransformation` class for further detail.

    :param dir1: in-plane vector defining the local x-direction
    :param dir2: in-plane vector used to define the local y-direction
    :param angle: rotation angle in degrees
    """

    def __init__(self, dir1=None, dir2=None, angle=None):
        super().__init__(dir1, dir2, angle)
