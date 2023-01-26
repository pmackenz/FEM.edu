import numpy as np

class Node():
    """
    class: representing a single Node
    """


    def __init__(self, x0, y0, z0=None):
        """

        :param x0: Initial position (List)

        """
        if isinstance(z0, (int, float)):
            self.pos = np.array([x0, y0, z0])
        else:
            self.pos = np.array([x0, y0])

        self.index    = -1
        self.disp     = np.array([])
        self.dofs     = {}
        self.elements = []
        self.fixity   = []
        self.force    = np.array([])
        self._hasLoad = False

    def __str__(self):
        s = \
        """Node {}: {}
        x:{}, fix:{}, 
        P:{}, u:{}""".format(self.index, self.dofs,
                             self.pos, self.fixity, self.force, self.disp)
        return s

    def __repr__(self):
        return "Node{}({}, x={}, u={})".format(self.index, self.dofs, self.pos, self.disp)

    def request(self, dof_list):
        """
        send list or individual dof code. Common codes:

        .. list-table::
            :header-rows: 1

            * - code
              - description
            * - **ux**
              - displacement in x-direction
            * - **uy**
              - displacement in y-direction
            * - **uz**
              - displacement in z-direction
            * - **rx**
              - rotation about x-axis
            * - **ry**
              - rotation about y-axis
            * - **rz**
              - rotation about z-axis


        :param: dof_list ... list of dof-codes required by calling element
        """
        dof_idx = []
        for dof in dof_list:
            if dof not in self.dofs:
                self.dofs[dof] = len(self.dofs)
                self.fixity.append(False)
                self.force = np.append(self.force, 0)
                self.disp = np.append(self.disp, 0)
            dof_idx.append(self.dofs[dof])

        return tuple(dof_idx)


    def fixDOF(self, dofs):
        """

        :param idx:
        """
        if isinstance(dofs, str):
            self.fixity[self.dofs[dofs]] = True
        else:
            for dof in dofs:
                if isinstance(dof, str):
                    self.fixity[self.dofs[dof]] = True
                else:
                    raise TypeError("fix DOF using name or index")


    def __floordiv__(self, other):
        """

        :param other:
        :return: self
        """
        self.fixDOF(other)
        return self

    def isFixed(self, idx):
        """

        :param idx:
        """
        return self.fixity[idx]

    def setDisp(self, u, v, dof_list=None):
        """

        :param u:
        :param v:
        """

        self.disp = np.array([u, v])

    def getDisp(self, dof_list=None):
        """

        :return: nodal displacement vector
        """
        return self.disp

    def getPos(self, dof_list=None):
        """

        :return: initial position vector
        """

        return self.pos

    def getDeformedPos(self, dof_list=None, factor=1.0):
        """
        Return deformed position :math:`{\\bf x} = {\\bf X} + f \\: {\\bf u}`

        :param factor: deformation magnification factor, :math:`f`.
        :return: deformed position vector, :math:`{\\bf x}`.
        """
        return self.pos + factor * self.disp

    def addLoad(self, P, dof_list=None):
        self.force   += np.array([Px, Py])
        self._hasLoad = True

    def setLoad(self, load_list, dof_list=None):
        if isinstance(dof_list, str):
            idx = self.dofs[dof_list]
            self.force[idx] = load_list
        else:
            for (load, dof) in zip(load_list, dof_list):
                idx = self.dofs[dof]
                self.force[idx] = load
        self._hasLoad = True

    def resetLoad(self):
        self.disp = np.zeros(len(self.dofs))
        self._hasLoad = False

    def getLoad(self, dof_list=None):
        return self.force

    def hasLoad(self):
        return self._hasLoad

    def resetDisp(self):
        self.disp = np.zeros(len(self.dofs))

    def resetAll(self):
        self.resetDisp()
        self.resetLoad()


if __name__ == "__main__":
    # testing the Node class
    node = Node(2.0, 3.5)
    node.index = 42
    node.setLoad(1.2, 3.4)
    node.addLoad(5.6, 7.8)
    node.setDisp(0.1234, -4.321)
    node.fixDOF(1)   # fixes y-direction

    print(repr(node))
    print(node)
