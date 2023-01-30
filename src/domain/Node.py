import numpy as np

class Node():
    """
    class: representing a single Node
    """

    def __init__(self, x ,y, u=0, v=0):
        """

        :param x:
        :param y:
        :param u:
        :param v:
        """
        self.pos      = np.array([x,y])
        self.index    = -1
        self.disp     = np.array([u,v])
        self.fixity   = [False, False]
        self.force    = np.zeros(2)
        self._hasLoad = False

    def __str__(self):
        s = \
"""Node {}:
   x:{}   y:{}
   fix:{} fix:{}
   Px:{}  Py:{}
   u:{}   v:{}""".format( self.index,
                           self.pos[0],   self.pos[1],
                           *self.fixity,
                           self.force[0], self.force[1],
                           self.disp[0],  self.disp[1])
        return s

    def __repr__(self):
        return "Node({},{},u={},v={})".format(self.pos[0],self.pos[1],self.disp[0],self.disp[1])

    def request(self, dof_list):
        """
        send list of dof codes. Common codes:

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
        pass

    def fixDOF(self, idx):
        """

        :param idx:
        """
        self.fixity[idx] = True

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
        self.disp = np.array([u,v])

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

    def addLoad(self, Px, Py, dof_list=None):
        self.force   += np.array([Px, Py])
        self._hasLoad = True

    def setLoad(self, Px, Py, dof_list=None):
        self.force    = np.array([Px, Py])
        self._hasLoad = True

    def getLoad(self):
        return self.force

    def hasLoad(self):
        return self._hasLoad

    def resetDisp(self):
        self.disp = np.zeros( 2 )

    def resetLoad(self):
        self.disp = np.zeros( 2 )
        self._hasLoad = False

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
