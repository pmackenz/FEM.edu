import numpy as np
from .Element import *
from domain.Node import *

class LinearTriangle(Element):
    """
    class: representing a single truss element
    """

    def __init__(self, node0, node1, node2, material):
        super().__init__((node0, node1, node2), material)
        self.element_type = DrawElement.TRIANGLE

        if node0.getPos().size == 3:
            dof_list = ('ux','uy','uz')
        elif node0.getPos().size == 2:
            dof_list = ('ux','uy')
        else:
            raise TypeError("dimension of nodes must be 2 or 3")

        self._requestDofs(dof_list)

        self.force    = 0.0
        self.Forces   = [ np.zeros(2), np.zeros(2) , np.zeros(2) ]
        self.Kt       = [ [np.zeros((2,2)), np.zeros((2,2)), np.zeros((2,2))],
                          [np.zeros((2,2)), np.zeros((2,2)), np.zeros((2,2))],
                          [np.zeros((2,2)), np.zeros((2,2)), np.zeros((2,2))] ]

        # covariant base vectors (reference system)
        base1 = node1.getPos() - node0.getPos()
        base2 = node2.getPos() - node0.getPos()
        self.gcov = np.vstack((base1, base2))

        # metric (reference system)
        self.GIJ = self.gcov @ self.gcov.T

        # dual base vectors (reference system)
        self.gcont = np.linalg.inv(self.GIJ) @ self.gcov

        self.area = np.sqrt(np.linalg.det(self.GIJ)) / 2.0

    def __str__(self):
        s = \
"""LinearTriangle: nodes {} {} {}:
   material properties: {}  strain:{}   stress:{}  
   internal force: {}
   Pe: [ {} {} ]""".format( self.nodes[0].index, self.nodes[1].index, self.nodes[2].index,
                            repr(self.material), self.material.getStrain(),
                            self.material.getStress(),
                            self.force, *self.Forces[1] )
        return s

    def __repr__(self):
        return "LinearTriangle({},{},{},{})".format( repr(self.nodes[0]),
                                                     repr(self.nodes[1]),
                                                     repr(self.nodes[2]),
                                                     repr(self.material))

    def updateState(self):

        node0 = self.nodes[0]
        node1 = self.nodes[1]
        node2 = self.nodes[2]

        Gs = self.gcont[0]
        Gt = self.gcont[1]
        Gu = -Gs - Gt

        # covariant base vectors (current system)

        gs = node1.getDeformedPos(self.dof_list) - node0.getDeformedPos(self.dof_list)
        gt = node2.getDeformedPos(self.dof_list) - node0.getDeformedPos(self.dof_list)
        gu = -gs - gt

        # metric (current system)
        gcov = np.vstack((gs, gt))
        gIJ = gcov @ gcov.T

        # deformation gradient
        F = np.outer(gs, Gs) + np.outer(gt, Gt)

        convective_strain = 0.5 * (gIJ - self.GIJ)
        almansi_strain = 0.5 * ( np.tensordot(F,F,((0,), (0,))) - np.eye(np.size(gs)) )

        # Dphi0_inv = np.linalg.inv(self.gcov.T)
        # val = Dphi0_inv.T @ convective_strain @ Dphi0_inv  # same as almansi_strain

        eps = almansi_strain

        # update the material state

        strain = {'xx':eps[0,0], 'yy':eps[1,1], 'xy':eps[0,1]+eps[1,0]}

        self.material.setStrain(strain)

        # 2nd Piola-Kirchhoff stress
        stress = self.material.getStress()

        S = np.array( [[stress['xx'],stress['xy']],[stress['xy'],stress['yy']]] )

        # 1st Piola-Kirchhoff stress
        P = F @ S

        # store stress for reporting
        self.stress = {'xx':P[0,0], 'xy':P[0,1], 'yx':P[1,0], 'yy':P[1,1]}

        # tractions
        ts = P @ Gs
        tt = P @ Gt
        tu = P @ Gu

        # initialize arrays
        gx = Gs[0] * gs + Gt[0] * gt
        gy = Gs[1] * gs + Gt[1] * gt

        # compute the kinematic matrices
        GI = (Gu, Gs, Gt)

        Bu = [Gu[0]*gx, Gu[1]*gy, Gu[1]*gx + Gu[0]*gy]
        Bs = [Gs[0]*gx, Gs[1]*gy, Gs[1]*gx + Gs[0]*gy]
        Bt = [Gt[0]*gx, Gt[1]*gy, Gt[1]*gx + Gt[0]*gy]

        BI = ( np.array(Bu), np.array(Bs), np.array(Bt) )

        # internal force
        self.Forces = [
            tu * self.area * self.material.getThickness(),
            ts * self.area * self.material.getThickness(),
            tt * self.area * self.material.getThickness()
            ]

        # tangent stiffness
        Ct = self.material.getStiffness()

        Kt = []
        One = np.eye(2, dtype=np.float64)
        for Gi, Bi in zip(GI, BI):
            Krow = []
            Ti = Gi @ S
            for Gj, Bj in  zip(GI, BI):
                GIJ = Ti @ Gj
                Krow.append( Bi.T @ Ct @ Bj + GIJ * One)
            Kt.append( Krow )

        self.Kt = Kt

    def getStress(self):
        return self.Stress


if __name__ == "__main__":
    # testing the Element class
    nd0 = Node(0.0, 0.0)
    nd0.index = 0
    nd1 = Node(3.0, 2.0)
    nd1.index = 1
    nd2 = Node(2.0, 4.0)
    nd2.index = 2
    params = {'E':100, 'A':1.5, 'fy':1.0e20}
    elem = LinearTriangle(nd0, nd1, nd2, Material(params))

    print(nd0)
    print(nd1)
    print(nd2)

    print("stress =", elem.getStress())
    print("nodal forces: ", *elem.getForce())
    print("element stiffness: ", elem.getStiffness())

    # change the nodal displacements
    nd0.setDisp(.1, .05)
    nd1.setDisp(.05, .2)
    nd2.setDisp(-.05, .1)

    print(nd0)
    print(nd1)
    print(nd2)

    print("force =", elem.getAxialForce())
    print("nodal forces: ", *elem.getForce())
    print("element stiffness: ", elem.getStiffness())


