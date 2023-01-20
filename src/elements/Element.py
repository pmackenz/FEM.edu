import numpy as np
import os
import sys


class Element():
    """
    abstract class: representing a single generic element
    """

    def __init__(self, nodes, material):
        """

        :param nodes:
        :param material:
        """
        self.nodes    = nodes
        self.material = material
        self.dof_list = ('ux','uy')

        self.force    = 0.0
        self.Forces   = None
        self.Kt       = None

    def __str__(self):
        s = \
        """{}: node {} to node {}:
    material properties: {}  strain:{}   stress:{}  
    internal force: {}
    Pe: [ {} {} ]""".format( self.__class__,
                            self.nodes[0].index, self.nodes[1].index,
                            repr(self.material), self.material.getStrain(),
                            self.material.getStress(),
                            self.force, *self.Forces[1] )
        return s

    def __repr__(self):
        return "{}({},{},{})".format(self.__class__,
                                     repr(self.nodes[0]),
                                     repr(self.nodes[1]),
                                     repr(self.material))

    def initalizeDoFs(self):
        for node in nodes:
            pass


    def getForce(self):
        """

        :return:
        """
        self.updateState()
        return self.Forces

    def getStress(self):
        self.updateState()
        return None

    def getStiffness(self):
        """

        :return:
        """
        self.updateState()
        return self.Kt

    def updateState(self):
        """

        """
        msg = "{}(Element): updateState() method has not been implemented".format(self.__class__.__name__)
        raise NotImplementedError(msg)


if __name__ == "__main__":

    sys.path.insert(0, os.path.abspath(".."))

    from Node import *
    from materials import Material

    # testing the Element class
    nd0 = Node(0.0, 0.0)
    nd0.index = 0
    nd1 = Node(3.0, 2.0)
    nd1.index = 1
    params = {'E':100, 'A':1.5, 'fy':1.0e20}
    mat = Material(params)
    elem = Element(nd0, nd1, mat)

    print(nd0)
    print(nd1)

    print("force =", elem.getAxialForce())
    print("nodal forces: ", *elem.getForce())
    print("element stiffness: ", elem.getStiffness())

    # change the nodal displacements
    nd0.setDisp(.1, .05)
    nd1.setDisp(.05, .2)

    print(nd0)
    print(nd1)

    print("force =", elem.getAxialForce())
    print("nodal forces: ", *elem.getForce())
    print("element stiffness: ", elem.getStiffness())


