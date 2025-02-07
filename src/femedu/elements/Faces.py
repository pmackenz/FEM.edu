import numpy as np
import sys


class Faces():
    r"""
    abstract class: representing one face of a 2D or 3D element
    """

    def __init__(self, id, *nds, **kwargs):
        self.id = id
        self.nodes = nds
        self.num_nodes = len(nds)
        self.load = [ 0.0, 0.0 ]   # surface traction in local coordinates: [ pn, ps ] are normal and shear loads
        self.flux = 0.0            # surface flux, out-flux is positive

        if self.num_nodes > 0:
            self.dim = nds[0].getPos().size
        else:
            self.dim = 0

        self.initialize()

    def __str__(self):
        s = "{}_{}".format(self.__class__.__name__, self.id)
        return s

    def initialize(self):
        r"""
        This is a virtual method.  Any class derived from :py:class:`Faces` must implement this function.
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def setLoad(self, pn, ps):
        r"""
        Defines surface loading on this Face.
        Coordinates n and s are outward normal and tangent directions, respectively.
        The local tangent is defined from the first to the second (to the third) node.

        :param pn: normal force per unit length
        :param ps: tangential force per unit length
        """
        self.load = [ pn, ps ]

    def setFlux(self, qn, outflux=False):
        r"""
        Defines a scalar surface flux on this Face.
        Coordinates n is the outward normal.
        **qn** is positive if it is an in-flux

        :param qn: normal flux per unit length (or unit area in 3d)
        :param outflux: set to `True` if given value is an out-flux. (same as entering a negative in-flux)
        """
        if outflux:
            self.flux = -qn  # internally, we always consider qn an in-flux
        else:
            self.flux = qn  # internally, we always consider qn an in-flux

    def computeNodalForces(self):
        r"""
        This is a virtual method.  Any class derived from :py:class:`Faces` must implement this function.
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def computeNodalFlux(self):
        r"""
        This is a virtual method.  Any class derived from :py:class:`Faces` must implement this function.
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def isFace(self, X, N):
        r"""
        Abstract interface for test function.

        .. note::

            This method needs to be implemented by every derived class.

        :param X: position vector for a point
        :type X: np.array
        :param N: outward normal vector at **X**
        :type N: np.array
        :return: **True** if **X** and **N** match this face. **False** otherwise.
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)
