import numpy as np
from copy import deepcopy

from ..Element import *
from ...domain.Node import *
from ...utilities import TriangleShapes, TriangleIntegration, GPdataType

class Triangle6(Element):
    r"""
    class: representing a 6-noded plane triangle
    """

    def __init__(self, node0, node1, node2, node3, node4, node5, material, label=None):
        super().__init__((node0, node1, node2, node3, node4, node5), material, label=label)

        ndim = node0.getPos().size

        if ndim == 3:
            dof_list = ('ux','uy','uz')
            ndof = 3
        elif ndim == 2:
            dof_list = ('ux','uy')
            ndof = 2
        else:
            raise TypeError("dimension of nodes must be 2 or 3")

        self.initialize(
            type=DrawElement.TRIANGLE,
            dofs=dof_list
                        )

        # select shape functions and integrator for the quadratic triangle
        self.interpolation = TriangleShapes()
        integrator    = TriangleIntegration(order=2)

        (self.xis, self.wis) = integrator.parameters()

        self.gpData = [ GPdataType() for i  in range(len(self.xis)) ]

        gpt = 0
        gp2nd_map = []

        for xi, wi, gpData in zip(self.xis,self.wis,self.gpData):

            gpData.material = deepcopy(material)

            shape   = self.interpolation.shape(  # requesting shape function array
                order=2,          # polynomial order per direction: quadratic
                s=xi[0], t=xi[1], # local coordinates for current position
                n=(0,0))          # no derivative with respect to (s,t)

            dshape1 = self.interpolation.shape(  # d shape function / d xi1
                order=2,          # polynomial order per direction: quadratic
                s=xi[0], t=xi[1], # local coordinates for current position
                n=(1,0))          # first derivative with respect to s

            dshape2 = self.interpolation.shape(  # d shape function / d xi2
                order=2,          # polynomial order per direction: quadratic
                s=xi[0], t=xi[1], # local coordinates for current position
                n=(0,1))          # first derivative with respect to t

            # covariant base vectors (reference system)
            X = np.zeros(ndim)
            base1 = np.zeros(ndim)
            base2 = np.zeros(ndim)
            for N, dN1, dN2, node in zip(shape, dshape1, dshape2, self.nodes):
                X     +=  N  * node.getPos()
                base1 += dN1 * node.getPos()
                base2 += dN2 * node.getPos()

            gpData.X = X
            self.gcov = np.vstack((base1, base2))
            gpData.base = self.gcov

            # metric (reference system)
            GIJ = self.gcov @ self.gcov.T
            gpData.state['GIJ'] = GIJ

            # dual base vectors (reference system)
            gpData.dual_base = np.linalg.inv(GIJ) @ self.gcov

            gpData.J = np.sqrt(np.linalg.det(GIJ))

            # populate gauss-point to nodes map
            raw_map = self.interpolation.shape(  # requesting shape function array
                order=1,  # polynomial order per direction: quadratic
                s=xi[0], t=xi[1],  # local coordinates for current position
                n=(0, 0))  # no derivative with respect to (s,t)
            map = raw_map * wi * gpData.J
            gp2nd_map.append(map)
            gpt += 1

        self.ngpts      = gpt                    # number of gauss points
        self._gp2nd_map = np.array(gp2nd_map).T  # gauss-point to nodes map array


    def __str__(self):
        s = super(Triangle6, self).__str__()
        for gp, gpdata in enumerate(self.gpData):
            s += "\n    strain {}: xx={xx:.3e} yy={yy:.3e} xy={xy:.3e} zz={zz:.3e}".format(gp, **gpdata.material.getStrain())
            s += "\n    stress {}: xx={xx:.3e} yy={yy:.3e} xy={xy:.3e} zz={zz:.3e}".format(gp, **gpdata.material.getStress())
        if np.array(self.distributed_load).any():
            s += "\n    element forces added to node:"
            for i, P in enumerate(self.Loads):
                Pi = np.array(P) * self.loadfactor
                s += "\n        {}: {}".format(self.nodes[i].getID(), Pi)
        return s

    def setSurfaceLoad(self, face, pn, ps=0):
        r"""
        .. list-table::
            :header-rows: 1

            * - face ID
              - nodes defining that face
            * - 0
              - :py:obj:`node 0` to :py:obj:`node 4` to :py:obj:`node 1`
            * - 1
              - :py:obj:`node 1` to :py:obj:`node 5` to :py:obj:`node 2`
            * - 2
              - :py:obj:`node 2` to :py:obj:`node 6` to :py:obj:`node 0`


        :param face: face ID for the laoded face
        :param pn: magnitude of distributed normal load per area. Tension on a surface is positive.
        :param ps: magnitude of distributed shear load per area. Positive direction is defined as shown in the above table.
        """
        if face >= 0 and face <= 2:
            self.faces[face].setLoad(pn, ps)

    def resetLoads(self):
        super(Triangle6, self).resetLoads()

    def updateState(self):

        # initializes internal force and tangent stiffness to zero arrays of the appropriate size.
        self.reset_matrices()

        for xi, wi, gpData in zip(self.xis,self.wis,self.gpData):

            # grab pre-computed quantities

            Gs = gpData.dual_base[0]
            Gt = gpData.dual_base[1]

            # shape functions

            dshape1 = self.interpolation.shape(  # d shape function / d xi1
                order=2,  # polynomial order per direction: quadratic
                s=xi[0], t=xi[1],  # local coordinates for current position
                n=(1, 0))  # first derivative with respect to s

            dshape2 = self.interpolation.shape(  # d shape function / d xi2
                order=2,  # polynomial order per direction: quadratic
                s=xi[0], t=xi[1],  # local coordinates for current position
                n=(0, 1))  # first derivative with respect to t

            # grab covariant base vectors (reference system) ...
            gso = self.gcov[0]
            gto = self.gcov[1]
            gu0 = -gso - gto

            # initialize covariant base vectors (current system) ...
            gs = np.zeros_like(Gs)  #
            gt = np.zeros_like(Gt)

            # ... compute
            for dN1, dN2, node in zip(dshape1, dshape2, self.nodes):
                gs += dN1 * node.getDeformedPos(self)
                gt += dN2 * node.getDeformedPos(self)
            # gu = -gs - gt

            # deformation gradient
            F = np.outer(gs, Gs) + np.outer(gt, Gt)

            # strain
            eps = 0.5 * ( F + F.T ) - np.eye(np.size(gso))
            gpData.state['strain'] = eps

            # update the material state
            strain = {'xx':eps[0,0], 'yy':eps[1,1], 'xy':eps[0,1]+eps[1,0]}
            gpData.material.setStrain(strain)

            # stress
            stress = gpData.material.getStress()
            gpData.state['stress'] = stress

            S = np.array( [[stress['xx'],stress['xy']],[stress['xy'],stress['yy']]] )

            # initialize components of the B-matrix ...
            dshapeX = dshape1 * Gs[0] + dshape2 * Gt[0]
            dshapeY = dshape1 * Gs[1] + dshape2 * Gt[1]

            # initialize arrays
            gxo = Gs[0] * gso + Gt[0] * gto
            gyo = Gs[1] * gso + Gt[1] * gto

            # compute the kinematic matrices
            GI = np.vstack((dshapeX,dshapeY)).T

            BI = np.array([ [dNx*gxo, dNy*gyo, dNx*gyo + dNy*gxo] for (dNx,dNy) in zip(dshapeX,dshapeY) ])

            # internal force
            area = gpData.J * wi

            # stress * area
            stress_vec = np.array([stress['xx'], stress['yy'], stress['xy']]) * area

            # tangent stiffness
            Ct = gpData.material.getStiffness() * area

            for Krow, Fi, Gi, Bi in zip(self.Kt, self.Forces, GI, BI):
                Fi += Bi.T @ stress_vec
                KGi = Gi @ S * area
                for KIJ, Gj, Bj in  zip(Krow, GI, BI):
                    KGij = KGi @ Gj
                    KIJ += Bi.T @ Ct @ Bj + KGij * np.eye(np.size(gs))

        # .. applied element load (reference load)
        self.computeSurfaceLoads()

    def computeSurfaceLoads(self):
        r"""
        compute surface loads using faces

        This method should be called during :py:meth:`updateState()` by every
        element supporting surface loads

        """
        self.Loads = [ np.zeros_like(self.Forces[I]) for I in range(len(self.nodes)) ]

        for I, face in enumerate(self.faces):
            loads = face.computeNodalForces()

            # indexing: face-ID matches the start node
            M = I+3  # the middle node
            J = I+1  # the end node
            if J>2:
                J -= 3

            # add to element load vectors
            self.Loads[I] += loads[0]
            self.Loads[M] += loads[1]
            self.Loads[J] += loads[2]


    def getStress(self):
        stress = [ data.state['stress'] for data in self.gpData ]
        return stress

    def mapGaussPoints(self, var, target_node=None):
        r"""
        Initiate mapping of Gauss-point values to nodes.
        This method is an internal method and should not be called by the user.
        Calling that method explicitly will cause faulty nodal values.

        :param var: variable code for a variable to be mapped from Gauss-points to nodes
        :param target_node: pointer to a node.  If given, the element will map only to that node.  Default is map to all nodes.
        """
        stresses = ('sxx','syy','szz','sxy','syz','szx')
        membrane = ('nxx','nyy','nxy')
        strains  = ('epsxx','epsyy','epszz','epsxy','epsyz','epszx')

        values = np.zeros( self.ngpts )

        if var.lower() in stresses:
            key = var[1:3].lower()
            values = []
            for gpdata in self.gpData:   # gauss-point loop
                if key in gpdata.state['stress']:
                    thickness = gpdata.material.getThickness()
                    values.append(gpdata.state['stress'][key]/thickness)
                else:
                    values.append(0.0)

        if var.lower() in membrane:
            key = var[1:3].lower()
            values = []
            for gpdata in self.gpData:   # gauss-point loop
                if key in gpdata.state['stress']:
                    values.append(gpdata.state['stress'][key])
                else:
                    values.append(0.0)

        if var.lower() in strains:
            key = var[1:3].lower()
            values = []
            for gpdata in self.gpData:   # gauss-point loop
                if key in gpdata.state['strain']:
                    values.append(gpdata.state['strain'][key])
                else:
                    values.append(0.0)

        values = np.array(values)

        for node, map in zip(self.nodes, self._gp2nd_map):
            if target_node == None or node == target_node:
                wndi   = np.sum(map)
                val_wi = map @ values
                node._addToMap(wndi, val_wi)