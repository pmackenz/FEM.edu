import numpy as np

from ..Element import *
from ...domain.Node import *

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

        self.L0       = np.linalg.norm(self.nodes[1].getPos() - self.nodes[0].getPos())
        self.force    = 0.0
        self.Forces   = [np.zeros(dim), np.zeros(dim)]
        self.Kt       = [[np.zeros((dim, dim)), np.zeros((dim, dim))],
                         [np.zeros((dim, dim)), np.zeros((dim, dim))]]

        """
        Compute internal state, nodal forces, and tangent stiffness for the current state of deformation.
        """
        X0 = self.getPos(0)
        X1 = self.getPos(1)

        # local coordinate system
        Lvec = X1 - X0
        L = np.linalg.norm(Lvec)   # element length (undeformed)
        Nvec = Lvec / L            # axial vector (unity)

        # get initial material stiffness
        self.material.setStrain({'xx':0.0})     # set strain to 0.0 for the initial stiffness
        EA   = self.material.getStiffness() * self.material.getArea()

        # nodal stiffness and element stiffness
        ke = (EA / L) * np.outer(Nvec, Nvec)
        self.Kt = [[ke, -ke], [-ke, ke]]

        # remember for force recovery
        self.L    = L
        self.Nvec = Nvec
        self.EA   = EA


    def __str__(self):
        area = self.material.getArea()
        stress = self.force / area
        strain = self.force / self.EA
        s = \
"""Truss: {} to {}:
    material properties: EA:{}  strain:{}   stress:{}  
    internal force: {}""".format( self.nodes[0].getID(), self.nodes[1].getID(),
                            self.EA, strain, stress, self.force)
        return s

    def __repr__(self):
        return "Truss({},{},{})".format( repr(self.nodes[0].getID()),
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
            strain = self.force / self.EA
            val = np.array([strain,strain])
            return (s,val)

        else:
            return (np.empty([0]),np.empty([0]))

    def updateState(self):
        """
        Stress/force recovery
        """
        L    = self.L
        Nvec = self.Nvec
        EA   = self.EA

        U0 = self.getDisp(0)
        U1 = self.getDisp(1)

        # compute strain
        eps = Nvec @ (U1 - U0) / L

        # compute internal force
        self.force = EA * eps

        # nodal force vector
        Pe = self.force * Nvec
        self.Forces = [-Pe, Pe]

