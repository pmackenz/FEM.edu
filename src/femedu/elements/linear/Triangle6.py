import numpy as np
from copy import deepcopy

from ..Element import *
from ...domain.Node import *
from ...utilities import TriangleShapes, TriangleIntegration, GPdataType

class Triangle6(Element):
    """
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
        self.integrator    = TriangleIntegration()

        (xis, wis) = self.integrator.parameters()

        self.gpData = [ GPdataType() for i  in range(len(xis)) ]

        for xi, wi, gpData in zip(xis,wis,self.gpData):

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
            gcov = np.vstack((base1, base2))
            gpData.base = gcov

            # metric (reference system)
            GIJ = gcov @ gcov.T
            gpData.state['GIJ'] = GIJ

            # dual base vectors (reference system)
            gpData.dual_base = np.linalg.inv(GIJ) @ gcov

            gpData.J = np.sqrt(np.linalg.det(GIJ))


    def __str__(self):
        s = super(Triangle6, self).__str__()
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

        (xis, wis) = self.integrator.parameters()

        for xi, wi, gpData in zip(xis,wis,self.gpData):

            # grab pre-computed quantities

            Gs = gpData.dual_base[0]
            Gt = gpData.dual_base[1]
            Gu = -Gs - Gt

            # shape functions

            dshape1 = self.interpolation.shape(  # d shape function / d xi1
                order=2,  # polynomial order per direction: quadratic
                s=xi[0], t=xi[1],  # local coordinates for current position
                n=(1, 0))  # first derivative with respect to s

            dshape2 = self.interpolation.shape(  # d shape function / d xi2
                order=2,  # polynomial order per direction: quadratic
                s=xi[0], t=xi[1],  # local coordinates for current position
                n=(0, 1))  # first derivative with respect to t

            # initialize covariant base vectors (current system) ...
            gs = np.zeros_like(Gs) #
            gt = np.zeros_like(Gt)

            # ... compute
            for dN1, dN2, node in zip(dshape1, dshape2, self.nodes):
                gs += dN1 * node.getPos()
                gt += dN2 * node.getPos()

            gu = -gs - gt

            # deformation gradient
            F = np.outer(gs, Gs) + np.outer(gt, Gt)

            # strain
            eps = 0.5 * ( F + F.T ) - np.eye(np.size(gs))
            gpData.state['strain'] = eps

            # update the material state
            strain = {'xx':eps[0,0], 'yy':eps[1,1], 'xy':eps[0,1]+eps[1,0]}
            gpData.material.setStrain(strain)

            # stress
            stress = gpData.material.getStress()
            gpData.state['stress'] = stress

            S = np.array( [[stress['xx'],stress['xy']],[stress['xy'],stress['yy']]] )

            # tractions
            ts = S @ Gs
            tt = S @ Gt
            tu = S @ Gu

            # initialize components of the B-matrix ...
            Bs = dshape1
            Bt = dshape2

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
            area = gpData.J * wi

            # self.Forces = [
            #     tu * area,
            #     ts * area,
            #     tt * area
            #     ]

            # tangent stiffness
            Ct = gpData.material.getStiffness() * area

            for Krow, Gi, Bi in zip(self.Kt, GI, BI):
                for KIJ, Gj, Bj in  zip(Krow, GI, BI):
                    KIJ += Bi.T @ Ct @ Bj

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



