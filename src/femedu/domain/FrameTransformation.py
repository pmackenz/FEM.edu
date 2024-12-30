import numpy as np

from .CosseratTransformation import *

class FrameTransformation(CosseratTransformation):
    r"""
    A Transformation including displacement dofs and rotational dofs


    If vectors :math:`{\bf d}_1` and :math:`{\bf d}_2` are given, the local coordinate frame :math:`\{{\bf g}_x,{\bf g}_y\},{\bf g}_z\}` is defined as follows

    .. math::
        {\bf g}_x = \frac{1}{||{\bf d}_1||} \: {\bf d}_1

    .. math::
        \hat{\bf g}_y = {\bf d}_2 - ({\bf d}_2\cdot{\bf g}_x) {\bf g}_x

    .. math::
        {\bf g}_y = \frac{1}{||\hat{\bf g}_y||} \: \hat{\bf g}_y

    .. math::
        {\bf g}_z = {\bf g}_y \times {\bf g}_y

    If instead an axial vector, :math:`{\boldsymbol\omega}`, is given, the local coordinate frame is defined as

    .. math::
        {\bf g}_i = {\boldsymbol\Lambda} \: {\bf e}_i
        \qquad i=x,y,z

    where  :math:`\{{\bf e}_x,{\bf e}_y,{\bf e}_z\}` is the global cartesian coordinate frame
    and

    .. math::
        {\boldsymbol\Lambda} := \exp[{\boldsymbol\Omega}]

    is the rotation tensor, generated
    as the matrix exponent of the skew-symmetric tensor :math:`{\boldsymbol\Omega}`. The latter is uniquely
    defined by the following equivalence.

    .. math::
        {\boldsymbol\Omega} \cdot {\bf h} \equiv {\boldsymbol\omega} \times {\bf h}
        \qquad\text{for all}\quad {\bf h}\in\mathbb{R}^3

    :param dir1: in-plane vector, **dir1** == :math:`{\bf d}_1`, defining the local x-direction
    :param dir2: in-plane vector, **dir2** == :math:`{\bf d}_2`, used to define the local y-direction
    :param axis: axial vector, :math:`{\boldsymbol\omega}`, in degrees
    """

    def __init__(self, dir1=None, dir2=None, axis=None):
        super().__init__(dir1, dir2, axis)
