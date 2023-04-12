import numpy as np
import sys


class Faces():
    """
    abstract class: representing one face of a 2D or 3D element
    """

    def __init__(self, id, *nds, **kwargs):
        self.id = id
        self.nodes = nds
        self.num_nodes = len(nds)
        self.load = [ 0.0, 0.0 ]   # surface traction in local coordinates: [ pn, ps ] are normal and shear loads

        if self.num_nodes > 0:
            self.dim = nds[0].getPos().size
        else:
            self.dim = 0

        self.initialize()

    def __str__(self):
        s = "{}_{}".format(self.__class__.__name__, self.id)
        return s

    def initialize(self):
        """
        This is a virtual method.  Any class derived from :py:class:`Faces` must implement this function.
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def setLoad(self, pn, ps):
        """
        Defines surface loading on this Face.
        Coordinates n and s are outward normal and tangent directions, respectively.
        The local tangent is defined from the first to th esecond (to the third) node.

        :param pn: normal force per unit length
        :param ps: tangential force per unit length
        """
        self.load = [ pn, ps ]

    def computeNodalForces(self):
        """
        This is a virtual method.  Any class derived from :py:class:`Faces` must implement this function.
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def isFace(self, X, N):
        """
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
