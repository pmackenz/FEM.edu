from .Element import *
from domain.Node import *

class Truss(Element):
    """
    class: representing a single 2D or 3D truss element

        self.nodes    = nodes i and j (tuple) \n
        self.material = material parameters  \n
        self.force    = internal force (float) \n
        self.Forces   = internal force vectors (list of np.arrays) \n
        self.Kt       = tangent stiffness (list of np.arrays) \n
    """

    def __init__(self, nodei, nodej, material):
        super().__init__((nodei, nodej), material)

        dim = nodei.getPos().size

        if dim == 3:
            dof_list = ('ux', 'uy', 'uz')
        elif dim == 2:
            dof_list = ('ux', 'uy')
        else:
            raise TypeError("spatial dimension of nodes must be 2 or 3")

        self._requestDofs(dof_list)

        self.L0       = np.linalg.norm(self.nodes[1].getPos() - self.nodes[0].getPos())
        self.force    = 0.0
        self.Forces   = [np.zeros(dim), np.zeros(dim)]
        self.Kt       = [[np.zeros((dim, dim)), np.zeros((dim, dim))],
                         [np.zeros((dim, dim)), np.zeros((dim, dim))]]


    def __str__(self):
        s = \
"""Truss: node {} to node {}:
   material properties: {}  strain:{}   stress:{}  
   internal force: {}
   Pe: [ {} {} ]""".format( self.nodes[0].index, self.nodes[1].index,
                            repr(self.material), self.material.getStrain(),
                            self.material.getStress(),
                            self.force, *self.Forces[1] )
        return s

    def __repr__(self):
        return "Truss({},{},{})".format( repr(self.nodes[0]),
                                         repr(self.nodes[1]),
                                         repr(self.material))

    def getAxialForce(self):
        """

        :return: truss axial force (float)
        """
        self.updateState()
        return self.force

    def getDisp(self,node):
        """
        get displacement vector of a desired node of the element

        **NEEDS TO BE IMPLEMENTED, ONCE TRANSFORMATION CLASS IS COMPLETE**

        :param node: pointer to node
        """


    def updateState(self):
        """
        Update truss axial force, end forces and tangent stiffness

        """
        U0 = self.nodes[0].getDisp()
        X0 = self.nodes[0].getPos()
        U1 = self.nodes[1].getDisp()
        X1 = self.nodes[1].getPos()


        Lvec = (X1 + U1) - (X0 + U0)
        ell = np.linalg.norm(Lvec)
        Nvec = Lvec / ell

        eps = Nvec @ (U1 - U0) / ell
        self.material.setStrain({'xx':eps})
        stress = self.material.getStress()
        sig = stress['xx']
        area   = self.material.getArea()
        self.force = sig * area

        Pe = self.force * Nvec
        self.Forces = [-Pe, Pe]

        Et = self.material.getStiffness()
        ke = (Et * area / ell) * np.outer(Nvec, Nvec)
        self.Kt = [[ke, -ke], [-ke, ke]]


if __name__ == "__main__":
    # testing the Element class
    nd0 = Node(0.0, 0.0)
    nd0.index = 0
    nd1 = Node(3.0, 2.0)
    nd1.index = 1
    params = {'E':100, 'A':1.5, 'fy':1.0e20}
    elem = Truss(nd0, nd1, Material(params))

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


