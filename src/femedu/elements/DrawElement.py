
import numpy as np

class PlotCurve():
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []

    def __add__(self, Xt):
        if not isinstance(Xt, np.ndarray):
            raise TypeError

        if len(Xt)>=2:
            self.x.append(Xt[0])
            self.y.append(Xt[1])
        if len(Xt)>2:
            self.z.append(Xt[2])

        return self

    def asTuple(self):
        return (np.array(self.x), np.array(self.y), np.array(self.z))


class DrawElement():
    """
    class for drawing any element

    element types are:

    .. list-table::

       * - UNKNOWN
         -
       * - LINE
         - truss, simple frame
       * - CURVE
         - beam, nice frame
       * - TRIANGLE
         - plates, shells
       * - TETRAHEDRON
         - plates, shells
       * - QUAD
         - plates, shells
       * - BRICK
         - continuum

    """

    """
    defining element class parameters
    
    Add a new one for every general type of elements.
    """
    UNKNOWN     = 0x000000
    LINE        = 0x000001  # truss, simple frame
    CURVE       = 0x000002  # beam, nice frame
    TRIANGLE    = 0x000004  # plates, shells
    TETRAHEDRON = 0x000008  # plates, shells
    QUAD        = 0x000010  # plates, shells
    BRICK       = 0x000020  # continuum

    def __init__(self):
        self.element_type = self.UNKNOWN

    def isType(self, test_type):
        return (test_type & self.element_type)

    def draw(self, factor=0.0, modeshape=False, **kwargs):
        """
        Returns a series of coordinate vectors representing the _x_, _y_, and _z_ values
        of points for plotting the deformed shape the current element (:code:`self`).
        For the undeformed element, set :code:`factor=0.0`.

        :param factor:  magnification factor for deformation
        :returns:  a tuple (list) of coordinate vectors
        """
        kwargs['modeshape'] = modeshape

        if self.element_type == self.LINE:
            return self.drawLine(factor, **kwargs)
        elif self.element_type == self.CURVE :
            return self.drawCurve(factor, **kwargs)
        elif self.element_type == self.TRIANGLE :
            return self.drawTriangle(factor, **kwargs)
        elif self.element_type == self.TETRAHEDRON :
            return self.drawTetrahedron(factor, **kwargs)
        elif self.element_type == self.QUAD :
            return self.drawQuad(factor, **kwargs)
        elif self.element_type == self.BRICK :
            return self.drawBrick(factor, **kwargs)
        else:
            raise NotImplementedError


    def drawLine(self, factor, **kwargs):
        """
        implementation of a generic :code:`LINE` type
        """
        c = PlotCurve()
        if len(self.nodes) >= 2:
            for node in self.nodes[:2]:
                c += node.getDeformedPos(self, factor, **kwargs)
        return c.asTuple()

    def drawCurve(self, factor, **kwargs):
        """
        implementation of a generic :code:`CURVE` type
        """
        Xi = self.nodes[0].getPos()
        Xj = self.nodes[1].getPos()

        #
        Lvec = Xj - Xi
        L = np.linalg.norm(Lvec)
        n = Lvec / L

        if n.shape[0] == 2:
            # assuming 2D
            s = np.array([-n[1], n[0]])
        elif n.shape[0] == 3:
            # assuming 2D
            k = np.linalg.cross(n, np.array([0.1, np.pi, 0.0]))
            kk = np.linalg.norm(k)
            if kk <= 0.1:
                k = np.linalg.cross(n, np.array([-0.1, np.pi, 0.0]))
                kk = np.linalg.norm(k)
            k /= kk
            s = np.linalg.cross(kk, n)
        else:
            return tuple()

        # local displacements
        # 0 ... ux
        # 1 ... uy
        # 2 ... theta
        dispi = self.nodes[0].getDisp(dofs=('ux','uy','rz'), **kwargs)
        dispj = self.nodes[1].getDisp(dofs=('ux','uy','rz'), **kwargs)

        # ui     = dispi[0]    # the transformation was performed when we received dispi
        # vi     = dispi[1]    # the transformation was performed when we received dispi
        ui     = dispi[0:2] @ n
        vi     = dispi[0:2] @ s
        thetai = dispi[2] * L

        # uj     = dispj[0]    # the transformation was performed when we received dispi
        # vj     = dispj[1]    # the transformation was performed when we received dispj
        uj     = dispj[0:2] @ n
        vj     = dispj[0:2] @ s
        thetaj = dispj[2] * L

        qu = np.array([ui, uj]) * factor
        qv = np.array([vi, thetai, vj, thetaj]) * factor

        # shape functions
        xsi = np.linspace(0, 1, 11)
        phiU = [1. - xsi, xsi]
        phiV = [
            [1., 0.972, 0.896, 0.784, 0.648, 0.5, 0.352, 0.216, 0.104, 0.028, 0.],
            [0., 0.081, 0.128, 0.147, 0.144, 0.125, 0.096, 0.063, 0.032, 0.009, 0.],
            [0., 0.028, 0.104, 0.216, 0.352, 0.5, 0.648, 0.784, 0.896, 0.972, 1.],
            [0., -0.009, -0.032, -0.063, -0.096, -0.125, -0.144, -0.147, -0.128, -0.081, 0.]
        ]

        # deformed shape local
        xl = xsi * L
        xl += qu @ phiU
        yl = qv @ phiV

        # deformed shape global
        xvecs = Xi + np.outer(xl, n) + np.outer(yl, s)

        return (xvecs[:, 0], xvecs[:, 1])


    def drawTriangle(self, factor, **kwargs):
        """
        implementation of a generic :code:`TRIANGLE` type
        """
        c = PlotCurve()
        if len(self.nodes) >= 3:
            for node in self.nodes[:3]:
                c += node.getDeformedPos(self, factor, **kwargs)
            c+= self.nodes[0].getDeformedPos(self, factor, **kwargs)
        return c.asTuple()

    def drawTetrahedron(self, factor, **kwargs):
        """
        implementation of a generic :code:`TETRAHEDRON` type
        """
        c = PlotCurve()
        if len(self.nodes) >= 4:
            for node in self.nodes:
                c += node.getDeformedPos(self, factor, **kwargs)
            c+= self.nodes[0].getDeformedPos(self, factor, **kwargs)
        return c.asTuple()

    def drawQuad(self, factor, **kwargs):
        """
        implementation of a generic :code:`QUAD` type
        """
        c = PlotCurve()
        if len(self.nodes) >= 4:
            for node in self.nodes[:4]:
                c += node.getDeformedPos(self, factor, **kwargs)
            c+= self.nodes[0].getDeformedPos(self, factor, **kwargs)
        return c.asTuple()

    def drawBrick(self, factor, **kwargs):
        """
        implementation of a generic :code:`BRICK` type
        """
        c = PlotCurve()
        if len(self.nodes) >= 8:
            for node in self.nodes[:4]:
                c += node.getDeformedPos(self, factor, **kwargs)
            c+= self.nodes[0].getDeformedPos(self, factor, **kwargs)
        return c.asTuple()



