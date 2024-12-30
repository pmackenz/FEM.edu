import numpy as np

from .FrameTransformation import *

class Frame2dTransformation(FrameTransformation):
    r"""
    A Transformation including displacement dofs and rotational dofs

    If vectors :math:`{\bf d}_1` and :math:`{\bf d}_2` are given, the local coordinate frame :math:`\{{\bf g}_x,{\bf g}_y\}` is defined as follows

    .. math::
        {\bf g}_x = \frac{1}{||{\bf d}_1||} \: {\bf d}_1

    .. math::
        \hat{\bf g}_y = {\bf d}_2 - ({\bf d}_2\cdot{\bf g}_x) {\bf g}_x

    .. math::
        {\bf g}_y = \frac{1}{||\hat{\bf g}_y||} \: \hat{\bf g}_y

    If instead an angle, :math:`\theta`, is given, the local coordinate frame is defined as

    .. math::
        {\bf g}_x = \cos\theta \: {\bf e}_x + \sin\theta \: {\bf e}_y

    .. math::
        {\bf g}_y = \sin\theta \: {\bf e}_x + \cos\theta \: {\bf e}_y

    where  :math:`\{{\bf e}_x,{\bf e}_y\}` is the global cartesian coordinate frame.

    :param dir1: in-plane vector, **dir1** == :math:`{\bf d}_1`, defining the local x-direction
    :param dir2: in-plane vector, **dir2** == :math:`{\bf d}_2`, used to define the local y-direction
    :param angle: rotation angle, :math:`\theta`, in degrees
    """

    def __init__(self, dir1=None, dir2=None, angle=None):

        dir1_is_vector = dir1 and (isinstance(dir1, np.ndarray) or isinstance(dir1, list))
        dir2_is_vector = dir2 and (isinstance(dir2, np.ndarray) or isinstance(dir2, list))

        if dir1_is_vector and len(dir1)==2:
            dir1 = np.append(dir1, [0.0])

        if dir2_is_vector and len(dir2)==2:
            dir2 = np.append(dir2, [0.0])

        if np.any(angle) and isinstance(angle,float):
            axis = np.array([0.,0.,np.radians(angle)])
        else:
            axis = None

        super().__init__(dir1, dir2, axis)
