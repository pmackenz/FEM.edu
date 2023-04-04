import numpy as np

from ..Element import *
from ...materials.Material import *
from ...domain.Node import *

class Frame2D(Element):
    """
    class: representing a 2D frame element

    **Assumptions**

    * The frame is in the X-Y-plane.
    * Small displacements and moderate rotations (:math:`P-\Delta` model)
    * Navier's and Bernoulli-Euler assumptions (plane sections and shear rigidity)
    * Linear elastic material (may change future releases)

    The element is using dofs :math:`u` (:code:`ux`), :math:`v` (:code:`uy`) and :math:`\\theta` (:code:`rz`)

    .. list-table:: Internal variables

        * - self.nodes
          - nodes i and j (tuple)
        * - self.material
          - pointer to Material object
        * - self.force
          - internal force (list of arrays for axial force, _f_, shear, _V_, and moment, _M_)
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
            dof_list = ('ux', 'uy', 'rz')
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
        s = super(Frame2D, self).__str__()
        if not 'Pw' in self.internal_forces:
            self.internal_forces['Pw'] = 0.0
        if not 'Mw' in self.internal_forces:
            self.internal_forces['Mw'] = 0.0
        s += "\n    internal forces: f0={fi:.2f} V0={Vi:.2f} M0={Mi:.2f} fl={fj:.2f} Vl={Vj:.2f} Ml={Mj:.2f} Pw={Pw:.2f} Mw={Mw:.2f}".format(**self.internal_forces)
        return s

    def setDistLoad(self, w):
        self.distributed_load = w

    def resetLoads(self):
        self.setDistLoad(0.0)
        super(Frame2D, self).resetLoads()

    def getAxialForce(self):
        self.updateState()
        return self.force

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

        elif variable.lower() == 'f' or variable.lower() == 'fx':
            # transverse shear (in-plane)
            Fl = self.internal_forces['fi']
            Fr = self.internal_forces['fj']
            s   = np.array([0.,1.])
            val = np.array([Fl, Fr])

        else:
            s   = np.array([0.,1.])
            val = np.zeros_like(s)

        return (s,val)

    def updateState(self):

        Xi = self.getPos(0)
        Xj = self.getPos(1)

        # 0 ... ux
        # 1 ... uy
        # 2 ... theta
        dispi = self.getDisp(0)
        Ui = dispi[:2]
        dispj = self.getDisp(1)
        Uj = dispj[:2]

        # compute local coordinate system
        Nvec = Xj - Xi
        L = np.linalg.norm(Nvec)
        Nvec /= L
        Svec = np.array([[0,-1],[1,0]]) @ Nvec

        nvec = (Xj + Uj) - (Xi + Ui)
        ell = np.linalg.norm(nvec)
        nvec /= ell
        svec = np.array([[0,-1],[1,0]]) @ nvec

        # - compute axial strain
        #strain = 0.5 * ((ell/L)**2 - 1.)                    # finite deformation
        #psi = Svec @ (Uj - Ui) / L                          # von Karman
        #strain = Nvec @ (Uj - Ui) / L + 0.5 * (psi * psi)   # von Karman
        strain = Nvec @ (Uj - Ui) / L                        # P-delta

        # - compute curvature
        """
        Treat flexure are linear elastic for now.
        This is achieved by feeding the material curvature==0
        and, thus, receiving the initial flexural tangent stiffness.
        The moment needs to be computed by alternative means.
        """
        #vi     = dispi[1]    # the transformation was performed when we received Ui
        #vj     = dispj[1]    # the transformation was performed when we received Uj

        vi = Ui @ Svec  # transformation not yet implemented
        vj = Uj @ Svec  # transformation not yet implemented
        thetai = dispi[2]
        thetaj = dispj[2]

        curvature = 0.0

        # - update material set strain
        self.material.setStrain({'axial':strain, 'flexure':curvature})

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

        # ** axial portion

        # - compute nodal tangent stiffness for axial
        n_tensor_n = np.outer(Nvec, Nvec)                                     # P-delta
        #n_tensor_n = np.outer((Nvec + psi/L * Svec), (Nvec + psi/L * Svec))  # von Karman
        #n_tensor_n = np.outer(nvec, nvec)                                    # finite deformation
        ke  = (EA / L) * n_tensor_n                                           # P-delta && von Karman
        #ke += self.force / L * (np.identity(len(nvec)) - n_tensor_n)         # von Karman
        #ke += self.force / L * (np.identity(len(nvec)) - n_tensor_n)         # finite deformation

        # ** bending portion (P-delta)

        if self.force < -EI / L**2 / 1000000:
            """
            use exact solution for compressive members
            """

            kappa  = np.sqrt(-self.force / EI) * L
            kappa2 = kappa / 2

            cot  = 1./np.tan(kappa)
            csc  = 1./np.sin(kappa)
            cot2 = 1./np.tan(kappa2)
            sec2 = 1./np.cos(kappa2)

            EIfact = EI / (L*L*L * (2. - kappa*cot2))
            kfu  = EIfact * kappa**3 * cot2
            kft  = EIfact * kappa**2 * L
            kmu  = kft

            EIfact = EI * kappa * csc / (L*( kappa*csc - sec2**2))
            kmt  = EIfact * (kappa*cot - 1.)
            kmth = EIfact * (1. - kappa*csc)

        elif self.force > EI / L**2 / 1000000:
            """
            use exact solution for tensile members 
            """
            kappa  = np.sqrt(self.force / EI) * L

            coshk = np.cosh(kappa)
            sinhk = np.sinh(kappa)

            EIfact = EI / (L*L*L*(2. - 2.*coshk + kappa*sinhk))

            kfu = EIfact * kappa**3 * sinhk
            kft = EIfact * kappa**2 * L * (coshk - 1.)
            kmu =  kft
            kmt =  EIfact * kappa * L*L * (kappa*coshk - sinhk)
            kmth = EIfact * kappa * L*L * (sinhk - kappa)

        else:
            """
            use the linear stiffness to avoid random numeric errors
            """

            kfu = 12 * EI / L**3
            kft =  6 * EI / L**2
            kmu =  kft
            kmt =  4 * EI / L
            kmth = kmt / 2

        #s_tensor_s = np.outer(svec, svec)

        n_tensor_n = np.outer(Nvec, Nvec)
        s_tensor_s = np.outer(Svec, Svec)
        KtII = np.block( [[ kfu * s_tensor_s,  kft * Svec[:,np.newaxis]], [ kmu * Svec, kmt ]] )
        KtIJ = np.block( [[-kfu * s_tensor_s,  kft * Svec[:,np.newaxis]], [-kmu * Svec, kmth]] )
        KtJI = np.block( [[-kfu * s_tensor_s, -kft * Svec[:,np.newaxis]], [ kmu * Svec, kmth]] )
        KtJJ = np.block( [[ kfu * s_tensor_s, -kft * Svec[:,np.newaxis]], [-kmu * Svec, kmt ]] )

        Vi =  kfu * (vi - vj) + kft * (thetaj + thetai)
        Vj = -kfu * (vi - vj) - kft * (thetaj + thetai)

        Mi = kmu * (vi - vj) + kmt * thetai + kmth * thetaj
        Mj = kmu * (vi - vj) + kmth * thetai + kmt * thetaj

        #Fi = Vi * svec - self.force * nvec                   # finite deformation
        #Fj = Vj * svec + self.force * nvec                   # finite deformation
        #Fi = Vi * Svec - self.force * (Nvec + psi/L * Svec)  # von Karman
        #Fj = Vj * Svec + self.force * (Nvec + psi/L * Svec)  # von Karman
        Fi = Vi * Svec - self.force * (Nvec)                  # P-delta
        Fj = Vj * Svec + self.force * (Nvec)                  # P-delta

        # ** combine the effects

        # - add nodal forces for axial
        KtII[0:2,0:2] +=  ke
        KtIJ[0:2,0:2] += -ke
        KtJI[0:2,0:2] += -ke
        KtJJ[0:2,0:2] +=  ke

        # ** build element load vector and element tangent stiffness

        # .. internal force
        self.Fi = np.r_[Fi, Mi]
        self.Fj = np.r_[Fj, Mj]
        self.Forces = [self.Fi, self.Fj]

        # .. tangent stiffness
        self.Kt = [[KtII, KtIJ],[KtJI, KtJJ]]

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

            Pi = np.r_[Pw * Svec, np.array([ Mw])]
            Pj = np.r_[Pw * Svec, np.array([-Mw])]
            self.Loads = (Pi, Pj)

            self.internal_forces['Pw'] = Pw * self.loadfactor
            self.internal_forces['Mw'] = Mw * self.loadfactor
