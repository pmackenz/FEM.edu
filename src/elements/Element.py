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
        self.transforms = [ None for nd in self.nodes ]
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


    def addTransformation(self, T, local_nodes=[]):
        """
        Attach a transformation to any node of the element.

        If no **local_nodes** list is given or an empty list is handed to the function,
        the transformation, **T**, will be applied to all nodes in the element.

        A non-empty **local_nodes** list will apply the transformation to those local nodes listed in that list.
        Local nodes start at 0 and go to N-1, where N is the number of elements in this element.

        A transformation can be removed from a node by assigning :code:`T=None` as the transformation.
        """
        if local_nodes:
            for local_id in local_nodes:
                if local_id >= 0 and local_id < len(self.transforms):
                    self.transforms[local_id] = T
        else:
            self.transforms = [ T for nd in self.nodes ]

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
            dof_idx = node.request(dof_requests, self)
            self.dof_idx[node] = dof_idx

    def getDofs(self):
        """
        returns the dof-codes for this element in a list
        """
        return self.dof_list

    def getCurve(self, factor=0.0):
        """
        Used for plotting this element.

        If a factor other than 0.0 (default: 0.0) is given,
        add factor*displacement to the reference position.

        Represent the element by an arbitrary polygon.
        Collect nodal coordinates of points in lists of X-, Y-, Z- coordinates.
        The Z-coordinate is only needed for 3d plotting.

        :return: plot information in a tuple (X,Y,...)
        """
        return None


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


