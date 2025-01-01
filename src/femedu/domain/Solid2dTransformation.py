import numpy as np

from .SolidTransformation import *

class Solid2dTransformation(SolidTransformation):
    r"""
    A Transformation including only displacement dofs (no rotational dofs)

    If vectors :math:`{\bf d}_1`  is given, the local coordinate frame :math:`\{{\bf g}_x,{\bf g}_y\}` is defined as follows.

    Normalized base

    .. math::
        {\bf g}_x = \frac{1}{||{\bf d}_1||} \: {\bf d}_1

    .. math::
        {\bf g}_y = {\bf e}_3 \times {\bf g}_x

    where :math:`{\bf e}_3` is the default out-of-plane unit vector.

    If instead an angle, :math:`\theta`, is given, the local coordinate frame is defined as

    .. math::
        {\bf g}_x = \cos\theta \: {\bf e}_x + \sin\theta \: {\bf e}_y

    .. math::
        {\bf g}_y = \sin\theta \: {\bf e}_x + \cos\theta \: {\bf e}_y

    where  :math:`\{{\bf e}_x,{\bf e}_y\}` is the global cartesian coordinate frame.

    :param dir1: in-plane vector, **dir1** == :math:`{\bf d}_1`, defining the local x-direction
    :param angle: rotation angle, :math:`\theta`, in degrees
    """

    def __init__(self, dir1=None, axis=None):

        dir1_is_vector = (isinstance(dir1, np.ndarray) or isinstance(dir1, list))

        if dir1_is_vector:
            if len(dir1)==2:
                dir2 = np.array([[0., -1.],[1., 0.],[0.0, 0.0]]) @ dir1
                dir1 = np.append(dir1, [0.0])      # vector perpendicular to the member axis
            elif len(dir1)==3:
                dir2 = np.array([[0., -1., 0.], [1., 0., 0.],[0., 0., 0.]]) @ dir1
            else:
                msg = ""
                raise TypeError(msg)

            dir2_is_vector = True

        super().__init__(dir1, dir2, axis)
