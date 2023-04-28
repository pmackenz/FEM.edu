import numpy as np
from ..Element import *
from ...domain.Node import *
from ...materials.Material import Material

class Triangle(Element):
    """
    class: representing a 3-node trangle for diffusion problems
    """

    def __init__(self, node0, node1, node2, material):
        super().__init__((node0, node1, node2), material)
        self.element_type = DrawElement.TRIANGLE
        self.createFaces()

        if not self.material.isMaterialType(Material.DIFFUSION):
            msg = f"Incompatible material type: {self.__class__.__name__} requires type Material.DIFFUSION"
            raise TypeError(msg)

        self._requestDofs('T')

        self.distributed_load = [0.0, 0.0, 0.0]
        self.force    = 0.0
        self.Forces   = [ np.zeros(1) for k in range(len(self.nodes)) ]
        self.Kt       = [ [ np.zeros(1) for k in range(len(self.nodes)) ] for m in range(len(self.nodes)) ]

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
        T0 = self.getDisp(0,'T')
        T1 = self.getDisp(1,'T')
        T2 = self.getDisp(2,'T')

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
        self.stress = {'delTx':flux[0], 'delTy':flux[1]}

        multiplier  = self.area
        multiplier *= self.material.getThickness()

        # internal force
        self.Forces = [
            Gu @ flux * multiplier,
            Gs @ flux * multiplier,
            Gt @ flux * multiplier
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



