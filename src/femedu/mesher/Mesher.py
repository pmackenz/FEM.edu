import sys
import numpy as np
from copy import deepcopy


class Mesher():
    """

    """

    def __init__(self, model, *pts):
        self.model = model
        self.offset = np.array([0.0, 0.0])
        self.nodes    = []
        self.elements = []

    def shift(self, dx, dy):
        """
        :param dx: move patch to the right
        :param dy: move patch up
        """
        self.offset += np.array([dx, dy])

    def lineMesh(self, NeX, element_type, material, **kwargs):
        """
        nD mesher using line elements

        This is an abstract class - requires implementation in subclass

        :param int NeX: number of elements along the curve
        :param element_type: compatible element type
        :param material: a material object
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def quadMesh(self, NeX, NeY, element_type, material, **kwargs):
        """
        2D mesher using quadrilateral elements

        This is an abstract class - requires implementation in subclass

        :param int NeX: number of elements along the first axis
        :param int NeY: number of elements along the second axis
        :param element_type: compatible element type
        :param material: a material object
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def triangleMesh(self, NeX, NeY, element_type, material, **kwargs):
        """
        2D mesher using triangular elements

        This is an abstract class - requires implementation in subclass

        :param int NeX: number of elements along the first axis
        :param int NeY: number of elements along the second axis
        :param element_type: compatible element type
        :param material: a material object
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def brickMesh(self, NeX, NeY, NeZ, element_type, material, **kwargs):
        """
        Solid mesher using brick elements

        This is an abstract class - requires implementation in subclass

        :param int NeX: number of elements along the first axis
        :param int NeY: number of elements along the second axis
        :param int NeZ: number of elements along the third axis
        :param element_type: compatible element type
        :param material: a material object
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def tetMesh(self, NeX, NeY, NeZ, element_type, material, **kwargs):
        """
        Solid mesher using tetrahedral elements

        This is an abstract class - requires implementation in subclass

        :param int NeX: number of elements along the first axis
        :param int NeY: number of elements along the second axis
        :param int NeZ: number of elements along the third axis
        :param element_type: compatible element type
        :param material: a material object
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def tie(self, other, tol=1.0e-3):
        """
        Check nodes of **other** Mesher object for (near) identical
        location to any node in this Mesher object.  If overlap is identified, make the
        node in **other** a follower of the respective node in this Mesher.

        :param other: Mesher instance
        """
        if isinstance(other, Mesher):
            # tie the patches together
            for nd1 in self.nodes:
                for nd2 in other.nodes:
                    if np.linalg.norm(nd2.getPos() - nd1.getPos()) < tol:
                        nd2.make_follower(nd1)
        else:
            msg = "The tie(other) method requires a Mesher object as argument"
            raise TypeError(msg)

