import numpy as np

from ..Element import *
from ...materials.Material import *
from ...domain.Node import *

class Beam2D(Element):
    """
    Representing a 2D beam element.

    **Assumptions**

    * The beam is defined along the X-axis.  (Y-axis values are ignored)
    * Small displacements and small rotations
    * Navier's and Bernoulli-Euler assumptions (plane sections and shear rigidity)
    * Linear elastic material (may change future releases)

    The element is using dofs :math:`v` (:code:`uy`) and :math:`\\theta=v'` (:code:`rz`)

    .. list-table:: Internal variables

        * - self.nodes
          - nodes i and j (tuple)
        * - self.material
          - pointer to Material object
        * - self.force
          - internal force (list of arrays for shear, _V_, and moment, _M_)
        * - self.Forces
          - nodal force vectors (list of np.arrays)
        * - self.Kt
          - tangent stiffness (list of np.arrays)

    """

    def __init__(self, nodei, nodej, material):
        """
        :param nodei: (pointer to) start Node object
        :param nodej: (pointer to) end Node object
        :param material: a material for axial and flexural integrated behavior
        """
        super().__init__((nodei, nodej), material)
        self.element_type = DrawElement.CURVE

        if not material.materialType() == Material.SECTION1D:
            raise TypeError("received material type {} but require SECTION1D type".format(material.materialType()))

        dim = nodei.getPos().size

        if dim == 2:
            dof_list = ('uy', 'rz')
        else:
            raise TypeError("spatial dimension of nodes must be 2")

        self._requestDofs(dof_list)

        ndof = len(dof_list)

        # initialize element load to zero
        self.distributed_load = 0.0

        self.L0       = np.linalg.norm(self.nodes[1].getPos() - self.nodes[0].getPos())
        self.force    = 0.0
        self.Forces   = [np.zeros(ndof), np.zeros(ndof)]
        self.Kt       = [[np.zeros((ndof, ndof)), np.zeros((ndof, ndof))],
                         [np.zeros((ndof, ndof)), np.zeros((ndof, ndof))]]
        self.internal_forces = {'fi':0.0, 'Vi':0.0, 'Mi':0.0, 'fj':0.0, 'Vj':0.0, 'Mj':0.0}

    def __str__(self):
        s  = "Beam2D: {} to {}:\n".format( self.nodes[0].getID(), self.nodes[1].getID())
        s += "   material {} properties: {}  strain:{}   stress:{}\n".format(self.material.__class__.__name__,
                                                                           self.material.parameters,
                                                                           self.material.getStrain(),
                                                                           self.material.getStress())
        s += "   nodal forces: Vi:{} Mi:{} Vj:{} Mj:{}".format(*self.Forces[0], *self.Forces[1])
        return s

    def __repr__(self):
        return "Beam2D({},{},{})".format( repr(self.nodes[0]),
                                          repr(self.nodes[1]),
                                          repr(self.material))

    def setDistLoad(self, w):
        """
        :param w: uniform load. positive if pointing in the local y-direction.
        """
        self.distributed_load = w

    def resetLoads(self):
        """

        """
        self.setDistLoad(0.0)
        super(Beam2D, self).resetLoads()

    def getInternalForce(self, variable=''):
        """
        computes vectors of normalized locations (**s**: ndarray) for which
        values (**val**: ndarray) are provided.
        Values of **s** are normalized to the interval :math:`[0,1]`.

        :returns: tuple (s, val)
        """
        self.updateState()

        Xi = self.nodes[0].getPos()
        Xj = self.nodes[1].getPos()
        Nvec = Xj - Xi
        L = np.linalg.norm(Nvec)

        if variable.lower() == 'm' or variable.lower() == 'mz':
            # bending moment (in plane)
            if self.distributed_load:
                s   = np.linspace(0,1,10)
            else:
                s   = np.array([0.,1.])

            Ml  = self.internal_forces['Mi']
            Ml += self.internal_forces['Mw']
            Mr  = self.internal_forces['Mj']
            Mr += self.internal_forces['Mw']

            #val = np.array([Ml, Mr])
            val = Ml * (1. - s) + Mr * s - 0.5 * (self.distributed_load * self.loadfactor) * L*L * s * (1 - s)

        elif variable.lower() == 'v' or variable.lower() == 'vy':
            # transverse shear (in-plane)
            if self.distributed_load:
                s   = np.linspace(0,1,10)
            else:
                s   = np.array([0.,1.])

            Vl  = self.internal_forces['Vi']
            Vl -= self.internal_forces['Pw']
            #Vr  = self.internal_forces['Vj']
            #Vr += self.internal_forces['Pw']

            val = Vl + s * (self.distributed_load * self.loadfactor) * L

        else:
            s   = np.array([0.,1.])
            val = np.zeros_like(s)

        return (s,val)

    def updateState(self):
        """
        Compute internal state, nodal forces, and tangent stiffness for the current state of deformation.
        """

        Xi = self.getPos(0)
        Xj = self.getPos(1)

        # 0 ... uy
        # 1 ... theta
        dispi = self.getDisp(0)
        dispj = self.getDisp(1)

        # compute local coordinate system
        L = Xj[0] - Xi[0]

        # - compute curvature
        """
        Treat flexure are linear elastic for now.
        This is achieved by feeding the material curvature==0
        and, thus, receiving the initial flexural tangent stiffness.
        The moment needs to be computed by alternative means.
        """
        vi     = dispi[0]    # the transformation was performed when we received Ui
        thetai = dispi[1]
        vj     = dispj[0]    # the transformation was performed when we received Uj
        thetaj = dispj[1]

        curvature = 0.0

        # - update material set strain
        self.material.setStrain({'axial':0.0, 'flexure':curvature})

        # - force and moment from the material
        stress = self.material.getStress()
        if 'axial' in stress:
            self.force = stress['axial']
        else:
            self.force = 0.0
        if 'flexure' in stress:
            self.moment = stress['flexure']
        else:
            self.moment = 0.0

        # get tangent moduli from the material
        Et = self.material.getStiffness()
        if 'ax-ax' in Et:
            EA = Et['ax-ax']
        else:
            EA = 0.0
        if 'flx-flx' in Et:
            EI = Et['flx-flx']
        else:
            EI = 0.0
        if ('ax-flx' in Et and Et['ax-flx']) or ('flx-ax' in Et and Et['flx-ax']):
            print('coupling of axial and flexural behavior is ignored by this element')

        #  use the linear stiffness

        kfu = 12. * EI / L**3
        kft =  6. * EI / L**2
        kmu =  kft
        kmt =  4. * EI / L
        kmth = kmt / 2.

        KtII = np.array( [[ kfu,  kft],
                          [ kmu, kmt ]] )
        KtIJ = np.array( [[-kfu,  kft],
                          [-kmu, kmth]] )
        KtJI = np.array( [[-kfu, -kft],
                          [ kmu, kmth]] )
        KtJJ = np.array( [[ kfu, -kft],
                          [-kmu, kmt ]] )

        # nodal forces
        Vi =  kfu * (vi - vj) + kft * (thetaj + thetai)
        Vj = -kfu * (vi - vj) - kft * (thetaj + thetai)

        Mi = kmu * (vi - vj) + kmt * thetai + kmth * thetaj
        Mj = kmu * (vi - vj) + kmth * thetai + kmt * thetaj


        # build element load vector and element tangent stiffness

        # .. internal force
        self.Fi = np.array([Vi, Mi])
        self.Fj = np.array([Vj, Mj])
        self.Forces = (self.Fi, self.Fj)

        # .. tangent stiffness
        self.Kt = np.array( [[KtII, KtIJ],[KtJI, KtJJ]] )

        # internal forces at nodes
        if 'Pw' in self.internal_forces:
            Pw = self.internal_forces['Pw']
        else:
            Pw = 0.0

        if 'Mw' in self.internal_forces:
            Mw = self.internal_forces['Mw']
        else:
            Mw = 0.0
        self.internal_forces = {'fi':self.force, 'Vi': Vi, 'Mi':-Mi,
                                'fj':self.force, 'Vj':-Vj, 'Mj': Mj,
                                'Pw':Pw, 'Mw':Mw}

    def computeSurfaceLoads(self):

        Xi = self.getPos(0)
        Xj = self.getPos(1)

        # compute local coordinate system
        Nvec = Xj - Xi
        L = np.linalg.norm(Nvec)
        Nvec /= L
        Svec = np.array([[0,-1],[1,0]]) @ Nvec

        # .. applied element load (reference load)
        if self.distributed_load:
            Pw = self.distributed_load * L / 2.
            Mw = Pw * L / 6.

            #Pi = np.r_[Pw * Svec, np.array([ Mw])]
            #Pj = np.r_[Pw * Svec, np.array([-Mw])]
            Pi = np.array([Pw, Mw])
            Pj = np.array([Pw, -Mw])
            self.Loads = (Pi, Pj)

            self.internal_forces['Pw'] = Pw * self.loadfactor
            self.internal_forces['Mw'] = Mw * self.loadfactor

    # prepare for removal
    def getAxialForce(self):
        raise DeprecationWarning

