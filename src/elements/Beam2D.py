from .Element import *
from materials.Material import *
from domain.Node import *

class Beam2D(Element):
    """
    class: representing a 2D beam element

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
            dof_list = ('uy', 'rz')
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
        s  = "Beam2D: node {} to node {}:\n".format( self.nodes[0].index, self.nodes[1].index)
        s += "   material {} properties: {}  strain:{}   stress:{}\n".format(self.material.__class__.__name__,
                                                                           self.material.parameters,
                                                                           self.material.getStrain(),
                                                                           self.material.getStress())
        s += "   nodal forces: Vi:{} Mi:{} Vj:{} Mi:{}".format(*self.Forces[0], *self.Forces[1])
        return s

    def __repr__(self):
        return "Beam2D({},{},{})".format( repr(self.nodes[0]),
                                         repr(self.nodes[1]),
                                         repr(self.material))

    def getAxialForce(self):
        return 0.0

    def updateState(self):

        Xi = self.nodes[0].getPos()
        Xj = self.nodes[1].getPos()

        # 0 ... uy
        # 1 ... theta
        dispi = self.nodes[0].getDisp()
        dispj = self.nodes[1].getDisp()

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
        self.Fi = np.array([Vi, Mi])
        self.Fj = np.array([Vj, Mj])
        self.Forces = (self.Fi, self.Fj)

        self.Kt = np.array( [[KtII, KtIJ],[KtJI, KtJJ]] )


if __name__ == "__main__":


    # testing the Element class
    nd0 = Node(0.0, 0.0)
    nd1 = Node(3.0, 0.0)
    params = {'E':100, 'A':1.5, 'I':1.0}
    elem = Beam2D(nd0, nd1, ElasticSection(params))

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


