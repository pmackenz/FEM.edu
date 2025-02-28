r"""
===================================================================
Patch test for the finite deformation 9-node quadrilateral
===================================================================
Full bi-quadratic interpolation

.. code::

             1
          x     y
      x^2   x*y   y^2
        x^2*y x*y^2
          x^2*y^2


"""

import numpy as np
from copy import deepcopy

from ..Element import *
from ...domain.Node import *
from ...utilities import QuadIntegration, QuadShapes, GPdataType

class Quad9(Element):
    r"""
    class: representing a plane 9-node quadrilateral

    This element works as 2D plate, loaded in-plane, and as a 3D membrane element.

    * For 2D plate behavior, define nodes as two-dimensional nodes
    * For 3D membrane behavior, define nodes as three-dimensional nodes
    """

    def __init__(self, node0, node1, node2, node3, node4, node5, node6, node7, node8, material, label=None):
        super(Quad9, self).__init__((node0, node1, node2, node3, node4, node5, node6, node7, node8), material, label=label)
        self.element_type = DrawElement.QUAD
        self.createFaces()

        ndim = node0.getPos().size

        if ndim == 3:
            dof_list = ('ux','uy','uz')
            ndof = 3
        elif ndim == 2:
            dof_list = ('ux','uy')
            ndof = 2
        else:
            raise TypeError("dimension of nodes must be 2 or 3")

        self._requestDofs(dof_list)

        self.distributed_load = [0.0, 0.0, 0.0, 0.0]   # face loads along the perimeter of the element
        self.force    = 0.0
        self.Forces   = [ np.zeros(ndof) for k in range(len(self.nodes)) ]
        self.Kt       = [ [ np.zeros((ndof,ndof)) for i in range(len(self.nodes)) ] for j in range(len(self.nodes)) ]
        self.ndof = ndof

        #self.material  = []   # one material object per integration point
        #self.stress    = []   # hold gauss-point stress (1st Piola-Kirchhoff stress)
        #self.strain    = []   # hold gauss-point strain
        #self.J         = []   # jacobian at integration point

        #self.Grad      = []   # derivative of shape functions with respect to global coords

        X  = np.array([ node.getPos() for node in self.nodes ])

        # initialization step
        self.integrator = QuadIntegration(order=4)
        self.xis, self.wis = self.integrator.parameters()

        self.gpData = [ GPdataType() for i  in range(len(self.xis)) ]

        gpt = 0
        gp2nd_map = []

        self.interpolation = QuadShapes()

        for xi, wi, gpData in zip(self.xis, self.wis, self.gpData):

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

            Grad = np.vstack((dshape1, dshape2))

            # covariant base vectors (reference system)
            Xgp   = np.zeros(ndim)
            base1 = np.zeros(ndim)
            base2 = np.zeros(ndim)
            for N, dN1, dN2, node in zip(shape, dshape1, dshape2, self.nodes):
                Xgp   +=  N  * node.getPos()
                base1 += dN1 * node.getPos()
                base2 += dN2 * node.getPos()

            gpData.X = Xgp
            gcov = np.vstack((base1, base2))
            gpData.base = gcov

            # metric (reference system)
            GIJ = gcov @ gcov.T
            gpData.state['GIJ'] = GIJ

            # dual base vectors (reference system)
            gpData.dual_base = np.linalg.inv(GIJ) @ gcov

            gpData.J = np.sqrt(np.linalg.det(GIJ))



            # reference configuration
            # -------------------------

            #
            # $ D\Phi_0 $
            #
            DPhi0 = (Grad @ X).T
            #self.J.append( np.linalg.det(DPhi0) )  # material model already includes thickness

            # dual base (contra-variant)
            gpData.Grad = np.linalg.inv(DPhi0).T @ Grad

            # populate gauss-point to nodes map
            raw_map= self.interpolation.shape(  # requesting shape function array
                order=2,                        # polynomial order per direction: quadratic
                s=xi[0], t=xi[1],               # local coordinates for current position
                n=(0, 0))                       # no derivative with respect to (s,t)
            map = raw_map * wi * gpData.J
            gp2nd_map.append(map)
            gpt += 1

        self.ngpts      = gpt                    # number of gauss points
        self._gp2nd_map = np.array(gp2nd_map).T  # gauss-point to nodes map array


    def __str__(self):
        s = super(Quad9, self).__str__()
        for igpt, gpData in enumerate(self.gpData):
            material = gpData.material
            s += "\n    strain ({}): xx={xx:.3e} yy={yy:.3e} xy={xy:.3e} zz={zz:.3e}".format(igpt,**material.getStrain())
            s += "\n    stress ({}): xx={xx:.3e} yy={yy:.3e} xy={xy:.3e} zz={zz:.3e}".format(igpt,**material.getStress())
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
              - :py:obj:`node 2` to :py:obj:`node 6` to :py:obj:`node 3`
            * - 3
              - :py:obj:`node 3` to :py:obj:`node 7` to :py:obj:`node 0`


        :param face: face ID for the laoded face
        :param pn: magnitude of distributed normal load per area. Tension on a surface is positive.
        :param ps: magnitude of distributed shear load per area. Positive direction is defined as shown in the above table.
        """
        if face >= 0 and face <= 3:
            self.faces[face].setLoad(pn, ps)

    def resetLoads(self):
        super(Quad9, self).resetLoads()

    def updateState(self):

        # initialization step
        ndim = self.nodes[0].getPos().size

        nnds = len(self.nodes)
        ndof = self.ndof       # mechanical element

        self.Forces = [ np.zeros(ndof) for k in range(nnds) ]
        Kt = [ [ np.zeros((ndof,ndof)) for i in range(nnds) ] for j in range(nnds) ]

        # create array of deformed nodal coordinates
        xt = np.array([ node.getDeformedPos() for node in self.nodes ])

        for xi, wi, gpData in zip(self.xis, self.wis, self.gpData):

            # grab pre-computed quantities
            Grad = gpData.Grad
            wi  *= gpData.J

            # spatial configuration
            # -------------------------

            # deformation gradient
            F = (Grad @ xt).T

            # strain
            eps = 0.5 * ( F.T @ F - np.eye(ndim) )

            # update the material state
            strain = {'xx':eps[0,0], 'yy':eps[1,1], 'xy':eps[0,1]+eps[1,0]}
            gpData.state['strain'] = strain
            gpData.material.setStrain(strain)

            # 2nd Piola-Kirchhoff stress
            stress =gpData.material.getStress()

            S = np.array( [[stress['xx'],stress['xy']],[stress['xy'],stress['yy']]] )

            # 1st Piola-Kirchhoff stress
            P = F @ S

            # store stress for reporting
            stress = {'xx':P[0,0], 'xy':P[0,1], 'yx':P[1,0], 'yy':P[1,1]}
            gpData.state['stress'] = stress

            # compute kinematic matrices
            BI = [ np.array([ Grad[0,K]*F[:,0],                       # the XX component
                              Grad[1,K]*F[:,1],                       # the YY component
                              Grad[0,K]*F[:,1] + Grad[1,K]*F[:,0] ])  # the XY component is XY + YX ("gammaXY = 2 epsXY")
                   for K in range(nnds) ]

            # internal forces
            for i, force in enumerate(self.Forces):
                # self.Forces[i] += P @ Grad[:,i] * wi  # wi already includes J
                force += P @ Grad[:,i] * wi  # wi already includes J

            # tangent stiffness
            Ct = gpData.material.getStiffness() * wi

            One = np.eye(2, dtype=np.float64)
            for I, Bi in enumerate(BI):
                Ti   = Grad[:,I] @ S
                GCti = Bi.T @ Ct
                for J, Bj in enumerate(BI):
                    GIJ = Ti @ Grad[:,J] * wi
                    Kt[I][J] += GCti @ Bj + GIJ * One

            # on to the next integration point

        self.Kt = Kt

    def computeSurfaceLoads(self):
        r"""
        compute surface loads using faces

        This method should be called during :py:meth:`updateState()` by every
        element supporting surface loads

        Nodes are mapped between the element and each Face2D as follows

        .. list-table:: Mapping element node IDs to face node IDs
            :header-rows: 1

            * - Face ID
              - Face nodes (0,1,2)
            * - 0
              - (0,4,1)
            * - 0
              - (1,5,2)
            * - 0
              - (2,6,3)
            * - 0
              - (3,7,0)


        """
        self.Loads = [ np.zeros_like(self.Forces[I]) for I in range(len(self.nodes)) ]

        for I, face in enumerate(self.faces):
            loads = face.computeNodalForces()

            # indexing
            J = I+1
            if J>3:
                J -= 4

            # add to element load vectors
            self.Loads[I] += loads[0]
            self.Loads[J] += loads[2]
            if loads.shape[0]>2:
                numNodes = len(self.nodes)
                if numNodes == 8 or numNodes == 9:
                    K = I + 4
                else:
                    msg = "Force data provided from Face2D inconsistent with element data"
                    raise TypeError(msg)

                self.Loads[K] += loads[1]

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