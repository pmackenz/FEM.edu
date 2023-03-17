import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splev, splrep

from .Mesher import *
from ..domain import Node, System

class CurveMesher(Mesher):
    """
    Mesher for beams and frame components.
    """

    def __init__(self, model, *pts):
        super(CurveMesher, self).__init__(model, *pts)

        if len(pts) < 2:
            msg = "{} requires at least 2 points but only {} given".format(self.__class__.__name__, len(pts))
            raise TypeError(msg)

        self.npts    = len(pts)
        self.ndim    = len(pts[0])
        self.support = np.array(pts)

    def mesh(self, Ne, element_type, material, **kwargs):
        """
        mesh a curve in 2d/3d using any kind of line element.

        The given set of points are interpolated using B-splines

        :param int Ne: number of elements along the curve
        :param element_type: compatible element type (BEam, Frame2D, ...)
        :param material: a 1D material object or section model
        :returns (nodes, elements):
        :var nodes: list of created :py:class:`Node` objects
        :var elements: list of created (subclass of) :py:class:`Element` objects
        """
        elements = []
        nodes    = []

        ndeg = np.min( [self.npts-1, 3] )
        ndim = self.ndim

        # source points local coordinates
        s = np.linspace(0, self.npts-1, self.npts)
        # target points local coordinates
        t = np.linspace(0, self.npts-1, Ne+1)

        interpolated = []

        for j in range(ndim):
            xsource = self.support[:,j]
            spl = splrep(s, xsource, k = ndeg)
            xtarget = splev(t, spl)
            interpolated.append(xtarget)

        fem_points = np.array(interpolated).T

        lastnode = None
        for coords in fem_points:
            node = Node(*coords)
            nodes.append(node)
            if lastnode:
                elements.append(element_type(lastnode, node, deepcopy(material)))
            lastnode = node

        if self.model and isinstance(self.model, System):
            self.model.addNode(*nodes)
            self.model.addElement(*elements)

        return (nodes, elements)

    def _create_doc(self):
        """
        create text and figure for the sphinx documentation
        """
        doc_string = """
        .. code:: python
        
            model = System()
            mesher = CurveMesher(model, (0,0),(1.5,.25),(2,1),(3.,1.5))
            mesher.mesh(10, Frame2D, ElasticSection())
        
        .. figure:: CurveMesher01.png
            :alt: CurveMesher
            :align: center
            
            Meshing a curve with 2-noded elements.
            
        """
        print(doc_string)

        points = [(0,0),(1.5,.25),(2,1),(3.,1.5)]
        support = np.array( points )
        npts = len(points)
        Ne = 10

        deg = np.min( [npts-1, 3] )

        s = np.linspace(0, npts-1, npts)
        x = support[:,0]
        y = support[:,1]
        splx = splrep(s, x, k = deg)
        sply = splrep(s, y, k = deg)

        s2 = np.linspace(0, self.npts-1, Ne+1)
        x2 = splev(s2, splx)
        y2 = splev(s2, sply)
        plt.plot(x, y, 'bo', ms=10, label="Input points")
        plt.plot(x2, y2, '-ro', label="Created mesh")
        plt.gca().set_aspect('equal')
        plt.gca().set_aspect('equal')
        plt.legend()
        plt.savefig("CurveMesher01.png", bbox_inches='tight')
        #plt.show()



