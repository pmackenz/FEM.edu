import numpy as np

from ..Element import *
from ...materials.Material import Material
from ...domain.Node import *

class Truss(Element):
    """
    Representing a single truss element for finite deformation analysis between 2 nodes.

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

        if not (self.material.materialType() == Material.SECTION1D
                or self.material.materialType() == Material.FIBER):
            msg = "Incompatible material type: need FIBER or SECTION1D"
            raise TypeError(msg)

        dim = nodei.getPos().size

        if dim == 3:
            dof_list = ('ux', 'uy', 'uz')
        elif dim == 2:
            dof_list = ('ux', 'uy')
        else:
            raise TypeError("spatial dimension of nodes must be 2 or 3")

        self._requestDofs(dof_list)

        # local coordinate system
        X0 = self.getPos(0)
        X1 = self.getPos(1)

        Lvec = X1 - X0
        self.L0   = np.linalg.norm(Lvec)
        self.Nvec = Lvec / self.L0

        # initialize forces and stiffness
        self.force    = 0.0
        self.Forces   = [np.zeros(dim), np.zeros(dim)]
        self.Kt       = [[np.zeros((dim, dim)), np.zeros((dim, dim))],
                         [np.zeros((dim, dim)), np.zeros((dim, dim))]]


    def __str__(self):
        s = \
"""Truss: {} to {}:
    material properties: {}  strain:{}   stress:{}  
    internal force: {}""".format( self.nodes[0].getID(), self.nodes[1].getID(),
                            repr(self.material), self.material.getStrain(),
                            self.material.getStress(), self.force)
        return s

    def __repr__(self):
        return "NLTruss({},{},{})".format( repr(self.nodes[0].getID()),
                                         repr(self.nodes[1].getID()),
                                         repr(self.material))

    def getAxialForce(self):
        """
        :returns: axial force in the truss. Positive value = tension.
        """
        self.updateState()
        return self.force

    def getInternalForce(self, variable=''):
        self.updateState()

        if variable.lower() == 'f' or variable.lower() == 'fx':
            # axial force
            s   = np.array([0.,1.])
            val = np.array([self.force, self.force])
            return (s,val)

        elif variable.lower() == 'eps' or variable.lower() == 'strain':
            # axial force
            s   = np.array([0.,1.])
            strain = self.material.getStrain()
            val = np.array([strain,strain])
            return (s,val)

        else:
            return (np.empty([0]),np.empty([0]))

    def updateState(self):
        """
        Compute internal state, nodal forces, and tangent stiffness for the current state of deformation.
        """
        U0 = self.getDisp(0)
        X0 = self.getPos(0)
        U1 = self.getDisp(1)
        X1 = self.getPos(1)

        # local coordinate system
        lvec = (X1 + U1) - (X0 + U0)
        ell = np.linalg.norm(lvec)
        nvec = lvec / ell

        # kinematics: Henky strain
        eps = np.log( ell/self.L0 )

        # constitutive behavior
        self.material.setStrain({'xx':eps})
        stress = self.material.getStress()
        sig    = stress['xx']

        # stress resultant
        area   = self.material.getArea()
        self.force = sig * area

        # nodal forces
        Pe = self.force * nvec
        self.Forces = [-Pe, Pe]

        # nodal and element stiffness matrix
        Et = self.material.getStiffness()
        n_outer_n = np.outer(nvec, nvec)
        ke = (Et * area / ell) * n_outer_n + self.force / ell * (np.eye(len(nvec)) - n_outer_n)
        self.Kt = [[ke, -ke], [-ke, ke]]
