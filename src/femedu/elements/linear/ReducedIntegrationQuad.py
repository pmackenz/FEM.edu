"""
==============================================
4-node quadrilateral - small displacement
==============================================
complete bi-linear interpolation
with selective reduced integration for shear

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

class ReducedIntegrationQuad(Element):
    """
    class: representing a plane 4-node quadrilateral

    This element works as 2D plate, loaded in-plane, and as a 3D membrane element.

    * For 2D plate behavior, define nodes as two-dimensional nodes
    * For 3D membrane behavior, define nodes as three-dimensional nodes
    """

    def __init__(self, node0, node1, node2, node3, material, label=None):
        super(ReducedIntegrationQuad, self).__init__((node0, node1, node2, node3), material, label=label)
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

        ## initialization step

        interpolation = QuadShapes()

        # reduced integration
        # -------------------------
        xi = (0.0, 0.0)   # gauss point at the center
        wi = 4.0          # weight at the center

        dphi_ds = interpolation.shape(1, *xi, n=(1, 0))
        dphi_dt = interpolation.shape(1, *xi, n=(0, 1))
        Grad = np.vstack((dphi_ds, dphi_dt))

        # reference configuration

        DPhi0 = (Grad @ X).T
        self.J0 = np.linalg.det(DPhi0)   # material model already includes thickness

        # dual base (contra-variant)
        self.Grad0 = np.linalg.inv(DPhi0).T @ Grad

        self.material0 = deepcopy(material)
        self.stress0   = {}

        # full integration
        # -------------------------
        integrator = QuadIntegration(order=2)
        self.xis, self.wis = integrator.parameters()

        gpt = 0
        gp2nd_map = []

        for xi, wi in zip(self.xis, self.wis):

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

            # populate gauss-point to nodes map
            map = interpolation.shape(1, *xi, n=(0,0)) * wi
            gp2nd_map.append(map)

            self.material.append(deepcopy(material))
            self.stress.append({})
            gpt += 1

        self.ngpts      = gpt                    # number of gauss points
        self._gp2nd_map = np.array(gp2nd_map).T  # gauss-point to nodes map array

    def __str__(self):
        s = super(ReducedIntegrationQuad, self).__str__()
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
        super(ReducedIntegrationQuad, self).resetLoads()

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

        # interpolation = QuadShapes()   # we are doing that and the isoparametric transformation in the constructor

        # reduced integration
        # -------------------------
        xi = (0.,0.)
        wi = 4.0

        Grad = self.Grad0  # pre-computed in __init__
        wi *= self.J0  # J includes the thickness of the plate

        # spatial configuration
        # -------------------------

        # deformation gradient
        F = (Grad @ xt).T

        # compute small strain tensor
        eps = 0.5 * (F + F.T) - np.eye(self.ndof)

        # update the material state
        strain0 = {'xx': eps[0, 0], 'yy': eps[1, 1], 'xy': eps[0, 1] + eps[1, 0]}

        self.material0.setStrain(strain0)

        # 2nd Piola-Kirchhoff stress
        stress = self.material0.getStress()

        S = np.array([[stress['xx'], stress['xy']], [stress['xy'], stress['yy']]])

        # store stress for reporting
        self.stress0 = {'xx': S[0, 0], 'xy': S[0, 1], 'yx': S[1, 0], 'yy': S[1, 1]}

        Ones = np.eye(self.ndof)
        # compute kinematic matrices
        BI = [np.array([Grad[0, K] * Ones[:, 1] + Grad[1, K] * Ones[:, 0]])   # the XY component is XY + YX ("gammaXY = 2 epsXY")
              for K in range(nnds)]

        # internal forces
        ##for i, force in enumerate(self.Forces):
        ##    force += S @ Grad[:, i] * wi  # wi already includes J

        # tangent stiffness
        Ct = self.material0.getStiffness() * wi

        for I, Bi in enumerate(BI):
            GCti = Bi.T * Ct[2,2]
            for J, Bj in enumerate(BI):
                Kt[I][J] += GCti @ Bj


        # full integration
        # -------------------------
        gpt = 0

        for xi, wi in zip(self.xis, self.wis):

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
            # ... replacing the shear strain at the gauss-points by the center shear strain
            #strain = {'xx':eps[0,0], 'yy':eps[1,1], 'xy':eps[0,1]+eps[1,0]}
            strain = {'xx':eps[0,0], 'yy':eps[1,1], 'xy':strain0['xy']}

            self.material[gpt].setStrain(strain)

            # 2nd Piola-Kirchhoff stress
            stress = self.material[gpt].getStress()

            S = np.array( [[stress['xx'],stress['xy']],[stress['xy'],stress['yy']]] )

            # store stress for reporting
            self.stress[gpt] = {'xx':S[0,0], 'xy':S[0,1], 'yx':S[1,0], 'yy':S[1,1]}

            Ones = np.eye(self.ndof)
            # compute kinematic matrices
            BI = [ np.array([ Grad[0,K]*Ones[:,0],                         # the XX component
                              Grad[1,K]*Ones[:,1]])                        # the YY component
                   for K in range(nnds) ]

            # internal forces
            for i, force in enumerate(self.Forces):
                force += S @ Grad[:,i] * wi  # wi already includes J

            # tangent stiffness
            Ct = self.material[gpt].getStiffness() * wi

            for I, Bi in enumerate(BI):
                GCti = Bi.T @ Ct[0:2,0:2]
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

        if var.lower() in membrane:
            key = var[1:3].lower()
            values = []
            for gpdata in self.stress:   # gauss-point loop
                if key in gpdata:
                    values.append(gpdata[key])
                else:
                    values.append(0.0)

        if var.lower() in stresses:
            key = var[1:3].lower()
            values = []
            for gpdata, material in zip(self.stress, self.material):   # gauss-point loop
                if key in gpdata:
                    thickness = material.getThickness()
                    values.append(gpdata[key]/thickness)
                else:
                    values.append(0.0)

        if var.lower() in strains:
            key = var[1:3].lower()
            values = []
            for gpdata in self.strain:   # gauss-point loop
                if key in gpdata:
                    values.append(gpdata[key])
                else:
                    values.append(0.0)

        values = np.array(values)

        for node, Ji, wi, map in zip(self.nodes, self.J, self.wis, self._gp2nd_map):
            if target_node == None or node == target_node:
                val_wi = map @ values
                node._addToMap(wi * Ji, val_wi * Ji)