"""
==============================================
4-node quadrilateral - small displacement
==============================================
incomplete bi-linear interpolation

.. code::

         1
      x     y
        x*y


"""
import numpy as np
from copy import deepcopy

from ..Element import *
from ...domain.Node import *
from ...utilities import QuadIntegration, QuadShapes

class Quad(Element):
    """
    class: representing a plane 4-node quadrilateral

    This element works as 2D plate, loaded in-plane, and as a 3D membrane element.

    * For 2D plate behavior, define nodes as two-dimensional nodes
    * For 3D membrane behavior, define nodes as three-dimensional nodes
    """

    def __init__(self, node0, node1, node2, node3, material):
        super(Quad, self).__init__((node0, node1, node2, node3), material)
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
        integrator = QuadIntegration(order=2)
        xis, wis = integrator.parameters()

        gpt = 0

        interpolation = QuadShapes()

        for xi, wi in zip(xis, wis):

            dphi_ds = interpolation.shape(1, *xi, n=(1,0))
            dphi_dt = interpolation.shape(1, *xi, n=(0,1))
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
        s = super(Quad, self).__str__()
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
              - :py:obj:`node 0` to :py:obj:`node 1`
            * - 1
              - :py:obj:`node 1` to :py:obj:`node 2`
            * - 2
              - :py:obj:`node 2` to :py:obj:`node 3`
            * - 3
              - :py:obj:`node 3` to :py:obj:`node 0`


        :param face: face ID for the laoded face
        :param pn: magnitude of distributed normal load per area. Tension on a surface is positive.
        :param ps: magnitude of distributed shear load per area. Positive direction is defined as shown in the above table.
        """
        if face >= 0 and face <= 3:
            self.faces[face].setLoad(pn, ps)

    def resetLoads(self):
        super(Quad, self).resetLoads()

    def updateState(self):

        # initialization step
        integrator = QuadIntegration(order=2)

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

            # compute small strain tensor
            eps = 0.5 * ( F + F.T ) - np.eye(self.ndof)

            # update the material state
            strain = {'xx':eps[0,0], 'yy':eps[1,1], 'xy':eps[0,1]+eps[1,0]}

            self.material[gpt].setStrain(strain)

            # 2nd Piola-Kirchhoff stress
            stress = self.material[gpt].getStress()

            S = np.array( [[stress['xx'],stress['xy']],[stress['xy'],stress['yy']]] )

            # store stress for reporting
            self.stress[gpt] = {'xx':S[0,0], 'xy':S[0,1], 'yx':S[1,0], 'yy':S[1,1]}

            Ones = np.eye(self.ndof)
            # compute kinematic matrices
            BI = [ np.array([ Grad[0,K]*Ones[:,0],                         # the XX component
                              Grad[1,K]*Ones[:,1],                         # the YY component
                              Grad[0,K]*Ones[:,1] + Grad[1,K]*Ones[:,0] ]) # the XY component is XY + YX ("gammaXY = 2 epsXY")
                   for K in range(nnds) ]

            # internal forces
            for i, force in enumerate(self.Forces):
                force += S @ Grad[:,i] * wi  # wi already includes J

            # tangent stiffness
            Ct = self.material[gpt].getStiffness() * wi

            for I, Bi in enumerate(BI):
                GCti = Bi.T @ Ct
                for J, Bj in enumerate(BI):
                    Kt[I][J] += GCti @ Bj

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



