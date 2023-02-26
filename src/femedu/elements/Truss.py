from .Element import *
from ..domain.Node import *

class Truss(Element):
    """
    Representing a single truss element between 2 nodes.

    The element is using the following dofs:

    * :math:`u` (:code:`ux`) and :math:`v` (:code:`uy`) in a 2D model.
    * :math:`u` (:code:`ux`), :math:`v` (:code:`uy`) and :math:`w` (:code:`uz`) in a 3D model.

    .. list-table:: Internal variables

        * - self.nodes
          - nodes i and j (tuple)
        * - self.material
          - pointer to Material object
        * - self.force
          - internal force (float)
        * - self.Forces
          - nodal force vectors (list of np.arrays)
        * - self.Kt
          - tangent stiffness (list of np.arrays)

    """

    def __init__(self, nodei, nodej, material):
        super().__init__((nodei, nodej), material)
        self.element_type = DrawElement.LINE

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
        :returns: axial force in the truss. Positive value = tension.
        """
        self.updateState()
        return self.force


    def updateState(self):
        """
        Compute internal state, nodal forces, and tangent stiffness for the current state of deformation.
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


