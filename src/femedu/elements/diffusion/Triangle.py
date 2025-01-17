import numpy as np
from ..LinearElement import *
from ...domain.Node import *
from ...materials.Material import Material

class Triangle(LinearElement):
    r"""
    class: representing a 3-node triangle for diffusion problems
    """

    def __init__(self, node0, node1, node2, material):

        if not material.isMaterialType(Material.DIFFUSION):
            msg = f"Incompatible material type: {self.__class__.__name__} requires type Material.DIFFUSION"
            raise TypeError(msg)

        super().__init__((node0, node1, node2), material)

        self.initialize(
            type=DrawElement.TRIANGLE,
            dofs=['T']
                        )

        # covariant base vectors (reference system)
        base1 = node1.getPos() - node0.getPos()
        base2 = node2.getPos() - node0.getPos()
        self.gcov = np.vstack((base1, base2))

        # metric (reference system)
        self.GIJ = self.gcov @ self.gcov.T

        self.area = np.sqrt(np.linalg.det(self.GIJ)) / 2.0

        # dual base vectors (reference system)
        if np.linalg.det(self.GIJ)<1.0e-10:
            print("Peter, we are having a problem")

        self.gcont = np.linalg.inv(self.GIJ) @ self.gcov

        # dual basis in triangle coordinates
        Gs = self.gcont[0]
        Gt = self.gcont[1]
        Gu = -Gs - Gt
        self.gcont = np.vstack((Gu, self.gcont))

        multiplier  = self.area
        multiplier *= self.material.getThickness()
        multiplier *= self.material.getDiffusivity()

        self.Kt = np.array( [[ Gu @ Gu, Gu @ Gs, Gu @ Gt ],
                             [ Gs @ Gu, Gs @ Gs, Gs @ Gt ],
                             [ Gt @ Gu, Gt @ Gs, Gt @ Gt ]] ) * multiplier


    def __str__(self):
        s = super(Triangle, self).__str__()
        gradPhi = self.material.getGrad()
        s += "\n    grad phi: x={:.3e} y={:.3e}".format(gradPhi[0], gradPhi[1])
        flux = self.material.getFlux()
        s += "\n    flux:     x={:.3e} y={:.3e}".format(flux[0], flux[1])
        if np.array(self.distributed_load).any():
            s += "\n    element flux added to node:"
            for i, P in enumerate(self.Loads):
                Pi = np.array(P) * self.loadfactor
                s += "\n        {}: {}".format(self.nodes[i].getID(), Pi)
        return s

    def setSurfaceLoad(self, face, qn):
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


        :param face: face ID for the loaded face
        :param qn: magnitude of distributed heat flux :math:`q_n`. Influx on a surface is positive.
        """
        if face >= 0 and face <= 2:
            self.faces[face].setLoad(qn, 0.0)

    def resetLoads(self):
        super(Triangle, self).resetLoads()

    def updateState(self):

        # nodal potential/temperature
        T0 = self.getDisp(0,dofs='T')
        T1 = self.getDisp(1,dofs='T')
        T2 = self.getDisp(2,dofs='T')

        # dual base vector in triangle coordinates
        Gu = self.gcont[0]
        Gs = self.gcont[1]
        Gt = self.gcont[2]

        # gradient
        gradT = T0 * Gu + T1 * Gs + T2 * Gt

        # update the material state
        self.material.setGrad(gradT)

        # stress: this is actually the element flux
        flux = self.material.getFlux()
        self.stress = {'Tx':gradT[0], 'Ty':gradT[1], 'qx':flux[0], 'qy':flux[1]}

        multiplier  = self.area
        multiplier *= self.material.getThickness()

        # internal force
        self.Forces = [
            -Gu @ flux * multiplier,
            -Gs @ flux * multiplier,
            -Gt @ flux * multiplier
            ]

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
            loads = face.computeNodalFlux()

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

    def mapGaussPoints(self, var):
        r"""
        Initiate mapping of Gauss-point values to nodes.
        This method is an internal method and should not be called by the user.
        Calling that method explicitly will cause faulty nodal values.

        :param var: variable code for a variable to be mapped from Gauss-points to nodes
        """
        # this element has a single gauss-point at s = t = u = 1/3

        flux     = ('qx','qy','Tx','Ty')

        value = 0.0

        if var in flux:
            flow = self.stress
            if var in flow:
                value = flow[var]

        wi     = self.area / 3.
        val_wi = value * wi

        for node in self.nodes:
            node._addToMap(wi, val_wi)






