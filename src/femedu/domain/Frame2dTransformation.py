import numpy as np

from .FrameTransformation import *

class Frame2dTransformation(FrameTransformation):
    r"""
    A Transformation including displacement dofs and rotational dofs

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
    :param angle: rotation angle, :math:`\theta`, in *degrees*
    """

    def __init__(self, dir1=None, angle=None):

        dir1_is_vector = (isinstance(dir1, np.ndarray) or isinstance(dir1, list))

        if dir1_is_vector:
            if len(dir1)==2:
                dir2 = np.array([[0., -1.],[1., 0.]]) @ dir1 # vector perpendicular to the member axis
            else:
                msg = f"Wrong spatial dimension of dir1: {dir1}"
                raise TypeError(msg)

            dir2_is_vector = True

        if np.any(angle) and isinstance(angle,float):
            axis = np.radians(angle)
        else:
            axis = None

        super().__init__(dir1, dir2, axis)

        # reset vectors for transformation
        self.known_vectors = []
        self.addVector(['ux','uy'])
        # self.addVector(['rz'])    # rotation does not require transformation
