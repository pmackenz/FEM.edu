"""
======================================
8-node serendipity quadrilateral
======================================
incomplete bi-quadratic interpolation

.. code::

             1
          x     y
      x^2   x*y   y^2
        x^2*y x*y^2


"""
import numpy as np
from copy import deepcopy

from ..Element import *
from ...domain.Node import *
from ...utilities import QuadIntegration, QuadShapes

class Quad8(Element):
    """
    class: representing a plane 8-node quadrilateral

    This element works as 2D plate, loaded in-plane, and as a 3D membrane element.

    * For 2D plate behavior, define nodes as two-dimensional nodes
    * For 3D membrane behavior, define nodes as three-dimensional nodes
    """

    def __init__(self, node0, node1, node2, node3, node4, node5, node6, node7, material):
        super(Quad8, self).__init__((node0, node1, node2, node3, node4, node5, node6, node7), material)
        self.element_type = DrawElement.QUAD
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

        self.distributed_load = [0.0, 0.0, 0.0, 0.0]   # face loads along the perimeter of the element
        self.force    = 0.0
        self.Forces   = [ np.zeros(ndof) for k in range(len(self.nodes)) ]
        self.Kt       = [ [ np.zeros(ndof) for k in range(len(self.nodes)) ] for m in range(len(self.nodes)) ]
        self.ndof = ndof

        self.material  = []   # one material object per integration point
        self.stress    = []   # hold gauss-point stress (1st Piola-Kirchhoff stress)
        self.J         = []   # jacobian at integration point

        self.Grad      = []   # derivative of shape functions with respect to global coords

        X  = np.array([ node.getPos() for node in self.nodes ])

        # initialization step
        integrator = QuadIntegration(order=3)
        xis, wis = integrator.parameters()

        gpt = 0

        interpolation = QuadShapes()

        for xi, wi in zip(xis, wis):

            dphi_ds = interpolation.shape(1, *xi, n=(1,0), serendipity=True)
            dphi_dt = interpolation.shape(1, *xi, n=(0,1), serendipity=True)
            Grad = np.vstack((dphi_ds, dphi_dt))

            # reference configuration
            # -------------------------

            #
            # $ D\Phi_0 $
            #
            DPhi0 = (Grad @ X).T
            self.J.append( np.linalg.det(DPhi0) )  # material model already includes thickness

            # dual base (contra-variant)
            self.Grad.append( np.linalg.inv(DPhi0).T @ Grad)

            self.material.append(deepcopy(material))
            self.stress.append({})
            gpt += 1

        pass


    def __str__(self):
        s = super(Quad8, self).__str__()
        for igpt, material in enumerate(self.material):
            s += "\n    strain ({}): xx={xx:.3e} yy={yy:.3e} xy={xy:.3e} zz={zz:.3e}".format(igpt,**material.getStrain())
            s += "\n    stress ({}): xx={xx:.3e} yy={yy:.3e} xy={xy:.3e} zz={zz:.3e}".format(igpt,**material.getStress())
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
        super(Quad8, self).resetLoads()

    def updateState(self):

        # initialization step
        integrator = QuadIntegration(order=3)

        nnds = len(self.nodes)
        ndof = self.ndof       # mechanical element

        self.Forces = [ np.zeros(ndof) for k in range(nnds) ]
        R  = [ np.zeros((ndof,ndof)) for i in range(nnds) ]
        Kt = [ [ np.zeros((ndof,ndof)) for i in range(nnds) ] for j in range(nnds) ]

        # create array of deformed nodal coordinates
        xt = np.array([ node.getDeformedPos() for node in self.nodes ])

        gpt = 0

        # interpolation = QuadShapes()   # we are doing that and the isoparametric transformation in the constructor

        xis, wis = integrator.parameters()

        for xi, wi in zip(xis, wis):

            # dphi_ds = interpolation.shape(1, *xi, n=(1,0))
            # dphi_dt = interpolation.shape(1, *xi, n=(0,1))
            # Grad = np.vstack((dphi_ds, dphi_dt))

            # reference configuration
            # -------------------------

            Grad = self.Grad[gpt]  # pre-computed in __init__
            wi *= self.J[gpt]      # J includes the thickness of the plate

            # spatial configuration
            # -------------------------

            # deformation gradient
            F = (Grad @ xt).T

            # compute Green-Lagrange strain tensor
            eps = 0.5 * ( np.tensordot(F,F,((0,), (0,))) - np.eye(self.ndof) )

            # update the material state
            strain = {'xx':eps[0,0], 'yy':eps[1,1], 'xy':eps[0,1]+eps[1,0]}

            self.material[gpt].setStrain(strain)

            # 2nd Piola-Kirchhoff stress
            stress = self.material[gpt].getStress()

            S = np.array( [[stress['xx'],stress['xy']],[stress['xy'],stress['yy']]] )

            # 1st Piola-Kirchhoff stress
            P = F @ S

            # store stress for reporting
            self.stress[gpt] = {'xx':P[0,0], 'xy':P[0,1], 'yx':P[1,0], 'yy':P[1,1]}

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
            Ct = self.material[gpt].getStiffness() * wi

            One = np.eye(2, dtype=np.float64)
            for I, Bi in enumerate(BI):
                Ti   = Grad[:,I] @ S
                GCti = Bi.T @ Ct
                for J, Bj in enumerate(BI):
                    GIJ = Ti @ Grad[:,J] * wi
                    Kt[I][J] += GCti @ Bj + GIJ * One

            # on to the next integration point
            gpt += 1

        self.Kt = Kt

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
            if J>3:
                J -= 4

            # add to element load vectors
            self.Loads[I] += loads[0]
            self.Loads[J] += loads[1]
            if loads.shape[0]>2:
                numNodes = len(self.nodes)
                if numNodes == 8 or numNodes == 9:
                    K = I + 4
                else:
                    msg = "Force data provided from Face2D inconsistent with element data"
                    raise TypeError(msg)

                self.Loads[K] += loads[2]

    def getStress(self):
        return self.Stress



