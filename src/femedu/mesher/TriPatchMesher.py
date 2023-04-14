from .Mesher import *
from .PatchMesher import *
from ..domain import Node, System
from ..elements.Element import Element


class TriPatchMesher(Mesher):
    """
    Mesher for plate and shell elements on a triangular domain,
    including in-plane and out-of plane loaded plates.
    """

    def __init__(self, model, *pts):
        super(TriPatchMesher, self).__init__(model, *pts)
        if len(pts)<3:
            raise TypeError("insufficient number of points given")

        X = list( pts[:3] )

        if len(pts) > 3 and pts[3]:
            X.append( pts[3] )
        else:
            Xm = 0.5*(np.array(X[0]) + np.array(X[1]))
            X.append( Xm )

        if len(pts) > 4 and pts[4]:
            X.append( pts[4] )
        else:
            Xm = 0.5*(np.array(X[1]) + np.array(X[2]))
            X.append( Xm )

        if len(pts) > 5 and pts[5]:
            X.append( pts[5] )
        else:
            Xm = 0.5*(np.array(X[2]) + np.array(X[0]))
            X.append( Xm )

        self.X = np.array(X)


    def quadMesh(self, Ne, element_type, material, **kwargs):
        """
        2D mesher using quadrilateral elements

        This is an abstract class - requires implementation in subclass

        :param int Ne: number of elements along each side
        :param element_type: compatible element type
        :param material: a material object
        :returns (nodes, elements):
        :var nodes: list of created :py:class:`Node` objects
        :var elements: list of created (subclass of) :py:class:`Element` objects
        """
        # if not element_type.isType(Element.QUAD):
        #     msg="Incompatible element type: {} in {}".format(element_type.__class__.__name__, sys._getframe().f_code.co_name)
        #     raise TypeError(msg)

        if Ne % 2:
            msg = "{} requires an even number elements per side. Ne={} was given.".format(self.__class__.__name__, Ne)
            raise TypeError(msg)

        # getting reference points

        # corners
        pt0 = self.map(0.0, 0.0)
        pt1 = self.map(1.0, 0.0)
        pt2 = self.map(0.0, 1.0)

        # mid points
        pt3 = self.map(0.5, 0.0)
        pt4 = self.map(0.5, 0.5)
        pt5 = self.map(0.0, 0.5)

        # center point
        pt6 = self.map(1./3., 1./3.)

        # quarter points on t=0
        pt7 = self.map(0.25, 0.0)
        pt8 = self.map(0.75, 0.0)

        # quarter points on u=0
        pt9 = self.map(0.75, 0.25)
        pt10 = self.map(0.25, 0.75)

        # quarter points on s=0
        pt11 = self.map(0.0, 0.75)
        pt12 = self.map(0.0, 0.25)

        # interior helper points
        pt13 = self.map(5./12., 2./12.)
        pt14 = self.map(5./12., 5./12.)
        pt15 = self.map(2./12., 5./12.)

        if self.model and isinstance(self.model, System):

            nodes    = []
            elements = []

            patch0 = PatchMesher(self.model, pt0, pt3, pt6, pt5, pt7, pt13, pt15, pt12)
            nds, elms = patch0.triangleMesh(Ne//2, Ne//2, element_type, material)
            nodes += nds
            elements += elms

            patch0 = PatchMesher(self.model, pt1, pt4, pt6, pt3, pt9, pt14, pt13, pt8)
            nds, elms = patch0.triangleMesh(Ne//2, Ne//2, element_type, material)
            nodes += nds
            elements += elms

            patch0 = PatchMesher(self.model, pt2, pt5, pt6, pt4, pt11, pt15, pt14, pt10)
            nds, elms = patch0.triangleMesh(Ne//2, Ne//2, element_type, material)
            nodes += nds
            elements += elms

            self.model.addNode(*nodes)
            self.model.addElement(*elements)

        self.nodes    = nodes
        self.elements = elements

        return (nodes, elements)

    def triangleMesh(self, Ne, element_type, material, **kwargs):
        """
        2D mesher using triangular elements

        This is an abstract class - requires implementation in subclass

        :param int Ne: number of elements along each side
        :param element_type: compatible element type
        :param material: a material object
        :returns (nodes, elements):
        :var nodes: list of created :py:class:`Node` objects
        :var elements: list of created (subclass of) :py:class:`Element` objects
        """
        # if not element_type.isType(Element.TRIANGLE):
        #     msg="Incompatible element type: {} in {}".format(element_type.__class__.__name__, sys._getframe().f_code.co_name)
        #     raise TypeError(msg)

        # local coordinates
        s = np.linspace(0.,1.,Ne+1)
        t = np.linspace(0.,1.,Ne+1)

        nodes    = []
        elements = []

        for j in range(Ne+1):
            row = []
            for i in range(Ne+1-j):
                # mapping from local to global coordinates
                X, Y = self.map(s[i],t[j])
                nd = Node(X,Y)
                row.append(nd)
            nodes.append(row)

        for j in range(Ne):
            for i in range(Ne-j):
                node_list = [ nodes[j][i], nodes[j][i+1], nodes[j+1][i] ]
                elem = element_type(*node_list, deepcopy(material))
                elements.append(elem)
                if j>0:
                    node_list = [ nodes[j][i+1], nodes[j][i], nodes[j-1][i+1] ]
                    elem = element_type(*node_list, deepcopy(material))
                    elements.append(elem)

        nodes = [ nd for row in nodes for nd in row ]

        if self.model and isinstance(self.model, System):
            self.model.addNode(*nodes)
            self.model.addElement(*elements)

        self.nodes    = nodes
        self.elements = elements

        return (nodes, elements)

    def map(self,s,t):
        """
        maps the local coordinates (s,t) from a bi-unit-square
        to the actual position (x,y) in the global model space
        """
        u = 1. - s - t
        shp = [ u*(2*u-1), s*(2*s-1), t*(2*t-1), 4*u*s, 4*s*t, 4*t*u ]
        x = shp @ self.X + self.offset

        return (x[0], x[1])

