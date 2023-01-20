from .Element import *
from domain.Node import *

class Truss(Element):
    """
    class: representing a single truss element

        self.nodes    = nodes i and j (tuple)
        self.material = material parameters (

        self.force    = internal force (float)
        self.Forces   = internal force vectors (list of np.arrays)
        self.Kt       = tangent stiffness (list of np.arrays)
    """

    def __init__(self, nodei, nodej, material, dim=3):
        super().__init__((nodei, nodej), material)

        # self.force    = 0.0
        # self.Forces   = [ np.zeros(2), np.zeros(2) ]
        # self.Kt       = [ [np.zeros((2,2)), np.zeros((2,2))], [np.zeros((2,2)), np.zeros((2,2))] ]


        self.setInitialLengthAndDirection()
        if nodei.getPos().size == 3:
            self.dof_list = ('ux','uy','uz')
        elif nodei.getPos().size == 2:
            self.dof_list = ('ux','uy')
        else:
            raise TypeError("dimension of nodes must be 2 or 3")



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
        self.updateState()
        return self.force

    def setInitialLengthAndDirection(self):
        X0 = self.nodes[0].getPos()
        X1 = self.nodes[1].getPos()

        Lvec = X1 - X0
        self.ell = np.linalg.norm(Lvec)
        self.Nvec = Lvec / self.ell

    def updateState(self):
        U0 = self.nodes[0].getDisp()
        U1 = self.nodes[1].getDisp()

        ell = self.ell
        Nvec = self.Nvec

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
        self.Kt = [[ke,-ke],[-ke,ke]]


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


