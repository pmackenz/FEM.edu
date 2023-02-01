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
        self.dof_idx  = {}

        self._requestDofs( tuple() )

        self.force    = 0.0
        self.Forces   = []
        self.Kt       = []


    def __str__(self):
        s = \
        """{}: nodes {}
        material properties: {}  
        strain:{},   stress:{},  internal force: {}
        Pe: {}""".format(self.__class__, self.nodes,
                         repr(self.material), self.material.getStrain(),
                         self.material.getStress(),
                         self.force, self.Forces)
        return s

    def __repr__(self):
        return "{}({},{},{})".format(self.__class__,
                                     repr(self.nodes[0]),
                                     repr(self.nodes[1]),
                                     repr(self.material))


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

    def _requestDofs(self, dof_requests):
        """
        Helper function (internal use) to inform **all** nodes of this element about the needed/used
        degrees of freedom.

        **Remark**: if nodes of different type are to be used by the element, **DO NOT** use this method but
        implement your own overloaded initialization within the constructor of your element.

        :param dof_requests: list of dofs for a typical node in this element
        """
        for node in self.nodes:
            dof_idx = node.request(dof_requests)
            node.linkElement(self)
            self.dof_idx[node] = dof_idx

    def getDofs(self):
        return self.dof_list


if __name__ == "__main__":

    sys.path.insert(0, os.path.abspath(".."))

    from domain import Node
    from materials import Material

    # testing the Element class
    nd0 = Node(0.0, 0.0)
    nd0.index = 0
    nd1 = Node(3.0, 2.0)
    nd1.index = 1
    params = {'E':100, 'A':1.5, 'fy':1.0e20}
    mat = Material(params)
    elem = Element([nd0, nd1], mat)

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


