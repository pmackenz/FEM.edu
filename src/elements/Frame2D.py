from .Element import *
from materials.Material import *
from domain.Node import *

class Frame2D(Element):
    """
    class: representing a 2D frame element

        self.nodes    = nodes i and j (tuple)
        self.material = material parameters (

        self.force    = internal force (float)
        self.Forces   = internal force vectors (list of np.arrays)
        self.Kt       = tangent stiffness (list of np.arrays)
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

        self.L0       = np.linalg.norm(self.nodes[1].getPos() - self.nodes[0].getPos())
        self.force    = 0.0
        self.Forces   = [np.zeros(ndof), np.zeros(ndof)]
        self.Kt       = [[np.zeros((ndof, ndof)), np.zeros((ndof, ndof))],
                         [np.zeros((ndof, ndof)), np.zeros((ndof, ndof))]]


    def __str__(self):
        s = \
"""Frame2D: node {} to node {}:
   material properties: {}  strain:{}   stress:{}  
   internal force: {}
   Pi:[ {} {} {} ]  Pj:[ {} {} {} ]""".format( self.nodes[0].index, self.nodes[1].index,
                            repr(self.material), self.material.getStrain(),
                            self.material.getStress(),
                            self.force, *self.Forces[0], *self.Forces[1] )
        return s

    def __repr__(self):
        return "Frame2D({},{},{})".format( repr(self.nodes[0]),
                                         repr(self.nodes[1]),
                                         repr(self.material))

    def getAxialForce(self):
        self.updateState()
        return self.force


    def updateState(self):

        Xi = self.nodes[0].getPos()
        Xj = self.nodes[1].getPos()

        # 0 ... ux
        # 1 ... uy
        # 2 ... theta
        dispi = self.nodes[0].getDisp()
        Ui = dispi[:2]
        dispj = self.nodes[1].getDisp()
        Uj = dispj[:2]

        # compute local coordinate system
        Nvec = Xj - Xi
        L = np.linalg.norm(Nvec)
        Nvec /= L

        nvec = (Xj + Uj) - (Xi + Ui)
        ell = np.linalg.norm(nvec)
        nvec /= ell
        svec = np.array([[0,-1],[1,0]]) @ nvec

        # - compute axial strain
        strain = 0.5 * ((ell/L)**2 - 1.)

        # - compute curvature
        """
        Treat flexure are linear elastic for now.
        This is achieved by feeding the material curvature==0
        and, thus, receiving the initial flexural tangent stiffness.
        The moment needs to be computed by alternative means.
        """
        vi     = dispi[1]    # the transformation was performed when we received Ui
        thetai = dispi[2]
        vj     = dispj[1]    # the transformation was performed when we received Uj
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
        n_tensor_n = np.outer(nvec, nvec)
        ke  = (EA / L) * n_tensor_n
        ke += self.force / L * (np.identity(len(nvec)) - n_tensor_n)

        # ** bending portion

        if self.force < -EI / L**2 / 1000:
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

        elif self.force > EI / L**2 / 1000:
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

        s_tensor_s = np.outer(svec, svec)
        KtII = np.block( [[ kfu * s_tensor_s,  kft * svec[:,np.newaxis]], [ kmu * svec, kmt ]] )
        KtIJ = np.block( [[-kfu * s_tensor_s,  kft * svec[:,np.newaxis]], [-kmu * svec, kmth]] )
        KtJI = np.block( [[-kfu * s_tensor_s, -kft * svec[:,np.newaxis]], [ kmu * svec, kmth]] )
        KtJJ = np.block( [[ kfu * s_tensor_s, -kft * svec[:,np.newaxis]], [-kmu * svec, kmt ]] )

        Vi =  kfu * (vi - vj) + kft * (thetaj + thetai)
        Vj = -kfu * (vi - vj) - kft * (thetaj + thetai)

        Mi = kmu * (vi - vj) + kmt * thetai + kmth * thetaj
        Mj = kmu * (vi - vj) + kmth * thetai + kmt * thetaj

        Fi = Vi * svec - self.force * nvec
        Fj = Vj * svec + self.force * nvec

        # ** combine the effects

        # - add nodal forces for axial
        KtII[0:2,0:2] +=  ke
        KtIJ[0:2,0:2] += -ke
        KtJI[0:2,0:2] += -ke
        KtJJ[0:2,0:2] +=  ke

        # ** build element load vector and element tangent stiffness
        self.Fi = np.r_[Fi, Mi]
        self.Fj = np.r_[Fj, Mj]

        self.Kt = [[KtII, KtIJ],[KtJI, KtJJ]]


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


