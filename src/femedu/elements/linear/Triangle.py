import numpy as np
from ..Element import *
from ...domain.Node import *

class Triangle(Element):
    """
    class: representing a 3-noded plane triangle
    """

    def __init__(self, node0, node1, node2, material, label=None):
        super().__init__((node0, node1, node2), material, label=label)
        self.element_type = DrawElement.TRIANGLE
        self.createFaces()

        if node0.getPos().size == 3:
            dof_list = ('ux','uy','uz')
            ndof = 3
        elif node0.getPos().size == 2:
            dof_list = ('ux','uy')
            ndof = 2
        else:
            raise TypeError("dimension of nodes must be 2 or 3")

        self._requestDofs(dof_list)

        self.distributed_load = [0.0, 0.0, 0.0]
        self.force    = 0.0
        self.Forces   = [ np.zeros(ndof) for k in range(len(self.nodes)) ]
        self.Kt       = [ [ np.zeros(ndof) for k in range(len(self.nodes)) ] for m in range(len(self.nodes)) ]

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
        s = super(Triangle, self).__str__()
        s += "\n    strain: xx={xx:.3e} yy={yy:.3e} xy={xy:.3e} zz={zz:.3e}".format(**self.material.getStrain())
        s += "\n    stress: xx={xx:.3e} yy={yy:.3e} xy={xy:.3e} zz={zz:.3e}".format(**self.material.getStress())
        if np.array(self.distributed_load).any():
            s += "\n    element forces added to node:"
            for i, P in enumerate(self.Loads):
                Pi = np.array(P) * self.loadfactor
                s += "\n        {}: {}".format(self.nodes[i].getID(), Pi)
        return s

    def setSurfaceLoad(self, face, pn, ps=0):
        """
        .. list-table::
            :header-rows: 1

            * - face ID
              - nodes defining that face
            * - 0
              - :py:obj:`node 0` to :py:obj:`node 1`
            * - 1
              - :py:obj:`node 1` to :py:obj:`node 2`
            * - 2
              - :py:obj:`node 2` to :py:obj:`node 0`


        :param face: face ID for the laoded face
        :param pn: magnitude of distributed normal load per area. Tension on a surface is positive.
        :param ps: magnitude of distributed shear load per area. Positive direction is defined as shown in the above table.
        """
        if face >= 0 and face <= 2:
            self.faces[face].setLoad(pn, ps)

    def resetLoads(self):
        super(Triangle, self).resetLoads()

    def updateState(self):

        node0 = self.nodes[0]
        node1 = self.nodes[1]
        node2 = self.nodes[2]

        Gs = self.gcont[0]
        Gt = self.gcont[1]
        Gu = -Gs - Gt

        # covariant base vectors (reference system)
        gso = self.gcov[0]
        gto = self.gcov[1]
        guo = -gso - gto

        # covariant base vectors (current system)
        gs = node1.getDeformedPos() - node0.getDeformedPos()
        gt = node2.getDeformedPos() - node0.getDeformedPos()
        gu = -gs - gt

        # metric (current system)
        # gcov = np.vstack((gs, gt))
        # gIJ = gcov @ gcov.T

        # deformation gradient
        F = np.outer(gs, Gs) + np.outer(gt, Gt)

        # strain
        eps = 0.5 * ( F + F.T ) - np.eye(np.size(gs))

        # update the material state
        strain = {'xx':eps[0,0], 'yy':eps[1,1], 'xy':eps[0,1]+eps[1,0]}
        self.material.setStrain(strain)

        # stress
        self.stress = self.material.getStress()

        S = np.array( [[self.stress['xx'],self.stress['xy']],[self.stress['xy'],self.stress['yy']]] )

        # tractions
        ts = S @ Gs
        tt = S @ Gt
        tu = S @ Gu

        # initialize arrays
        gxo = Gs[0] * gso + Gt[0] * gto
        gyo = Gs[1] * gso + Gt[1] * gto

        # compute the kinematic matrices
        GI = (Gu, Gs, Gt)

        Bu = [Gu[0]*gxo, Gu[1]*gyo, Gu[1]*gxo + Gu[0]*gyo]
        Bs = [Gs[0]*gxo, Gs[1]*gyo, Gs[1]*gxo + Gs[0]*gyo]
        Bt = [Gt[0]*gxo, Gt[1]*gyo, Gt[1]*gxo + Gt[0]*gyo]

        BI = ( np.array(Bu), np.array(Bs), np.array(Bt) )

        # internal force
        self.Forces = [
            tu * self.area,
            ts * self.area,
            tt * self.area
            ]

        # tangent stiffness
        Ct = self.material.getStiffness() * self.area

        Kt = []
        for Gi, Bi in zip(GI, BI):
            Krow = []
            for Gj, Bj in  zip(GI, BI):
                Krow.append( Bi.T @ Ct @ Bj )
            Kt.append( Krow )

        self.Kt = Kt

        # .. applied element load (reference load)
        self.computeSurfaceLoads()

    def computeSurfaceLoads(self):
        """
        compute surface loads using faces

        This method should be called during :py:meth:`updateState()` by every
        element supporting surface loads

        """
        self.Loads = [ np.zeros_like(self.Forces[I]) for I in range(len(self.nodes)) ]

        for I, face in enumerate(self.faces):
            loads = face.computeNodalForces()

            # indexing
            J = I+1
            if J>2:
                J -= 3

            # add to element load vectors
            self.Loads[I] += loads[0]
            self.Loads[J] += loads[1]
            if loads.shape[0]>2:
                numNodes = len(self.nodes)
                if numNodes == 6:
                    K = I + 3
                elif numNodes == 8 or numNodes == 9:
                    K = I + 4
                else:
                    msg = "Force data provided from Face2D inconsistent with element data"
                    raise TypeError(msg)

                self.Loads[K] += loads[2]

    def getStress(self):
        return self.Stress

    def mapGaussPoints(self, var, target_node=None):
        r"""
        Initiate mapping of Gauss-point values to nodes.
        This method is an internal method and should not be called by the user.
        Calling that method explicitly will cause faulty nodal values.

        :param var: variable code for a variable to be mapped from Gauss-points to nodes
        :param target_node: pointer to a node.  If given, the element will map only to that node.  Default is map to all nodes.
        """
        # this element has a single gauss-point at s = t = u = 1/3

        stresses = ('sxx','syy','szz','sxy','syz','szx')
        membrane = ('nxx','nyy','nxy')
        strains  = ('epsxx','epsyy','epszz','epsxy','epsyz','epszx')

        value = 0.0

        if var.lower() in stresses:
            key = var[1:3].lower()
            tensor = self.material.getStress()   # 2nd-Piola-Kirchhoff stress
            if key in tensor:
                thickness = self.material.getThickness()
                value = tensor[key]/thickness

        if var.lower() in membrane:
            key = var[1:3].lower()
            tensor = self.material.getStress()   # 2nd-Piola-Kirchhoff stress
            if key in tensor:
                value = tensor[key]

        elif var.lower() in strains:
            key = var[1:3].lower()
            tensor = self.material.getStrain()
            if key in tensor:
                value = tensor[key]

        wi     = self.area / 3.
        val_wi = value * wi

        for node in self.nodes:
            if target_node == None or node == target_node:
                node._addToMap(wi, val_wi)
