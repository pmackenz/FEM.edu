from dis import JUMP_BACKWARD

from fontTools.misc.textTools import caselessSort

from .Mesher import *
from ..domain import Node, System

class PatchMesher(Mesher):
    """
    Mesher for plate and shell elements on a quadrilateral domain,
    including in-plane and out-of plane loaded plates.
    """

    def __init__(self, model, *pts):
        super(PatchMesher, self).__init__(model, *pts)
        if len(pts)<4:
            raise TypeError("insufficient number of points given")

        X = list( pts[:4] )

        if len(pts) > 4 and pts[4]:
            X.append( pts[4] )
        else:
            Xm = 0.5*(np.array(X[0]) + np.array(X[1]))
            X.append( Xm )

        if len(pts) > 5 and pts[5]:
            X.append( pts[5] )
        else:
            Xm = 0.5*(np.array(X[1]) + np.array(X[2]))
            X.append( Xm )

        if len(pts) > 6 and pts[6]:
            X.append( pts[6] )
        else:
            Xm = 0.5*(np.array(X[2]) + np.array(X[3]))
            X.append( Xm )

        if len(pts) > 7 and pts[7]:
            X.append( pts[7] )
        else:
            Xm = 0.5*(np.array(X[3]) + np.array(X[0]))
            X.append( Xm )

        if len(pts) > 8 and pts[8]:
            X.append( pts[8] )
        else:
            Xm = 0.25*(np.array(X[4]) + np.array(X[5]) + np.array(X[6]) + np.array(X[7]))
            X.append( Xm )

        self.X = np.array(X)


    def quadMesh(self, NeX, NeY, element_type, material, **kwargs):
        """
        2D mesher using quadrilateral elements

        :param int NeX: number of elements along the first axis
        :param int NeY: number of elements along the second axis
        :param element_type: compatible element type
        :param material: a material object
        :returns (nodes, elements):
        :var nodes: list of created :py:class:`Node` objects
        :var elements: list of created (subclass of) :py:class:`Element` objects
        """
        # if not element_type.isType(Element.QUAD):
        #     msg="Incompatible element type: {} in {}".format(element_type.__class__.__name__, sys._getframe().f_code.co_name)
        #     raise TypeError(msg)

        match element_type.__name__:
            case 'Quad':
                return self._quad4_mesh(NeX, NeY, element_type, material, **kwargs)
            case 'Quad8':
                return self._quad8_mesh(NeX, NeY, element_type, material, **kwargs)
            case 'Quad9':
                return self._quad9_mesh(NeX, NeY, element_type, material, **kwargs)
            case _:
                msg = "{} does not support element type {}".format(self.__class__.__name__, element_type.__name__)
                raise TypeError(msg)


    def _quad4_mesh(self, NeX, NeY, element_type, material, **kwargs):

        # local coordinates
        ss = np.linspace(-1.,1.,NeX+1)
        tt = np.linspace(-1.,1.,NeY+1)

        nodes    = []
        elements = []

        for t in tt:
            row = []
            for s in ss:
                # mapping from local to global coordinates
                X, Y = self.map(s,t)
                nd = Node(X,Y)
                row.append(nd)
            nodes.append(row)

        for j in range(NeY):
            for i in range(NeX):
                node_list = [ nodes[j][i], nodes[j][i+1], nodes[j+1][i+1], nodes[j+1][i] ]
                elem = element_type(*node_list, deepcopy(material))
                elements.append(elem)

        nodes = [ nd for row in nodes for nd in row ]

        if self.model and isinstance(self.model, System):
            self.model.addNode(*nodes)
            self.model.addElement(*elements)

        self.nodes    = nodes
        self.elements = elements

        return (nodes, elements)

    def _quad8_mesh(self, NeX, NeY, element_type, material, **kwargs):
        msg = "{} does not yet support element type {}".format(self.__class__.__name__, element_type.__name__)
        raise TypeError(msg)

    def _quad9_mesh(self, NeX, NeY, element_type, material, **kwargs):

        # local coordinates
        ss = np.linspace(-1.,1.,2*NeX+1)
        tt = np.linspace(-1.,1.,2*NeY+1)

        nodes    = []
        elements = []

        for t in tt:
            row = []
            for s in ss:
                # mapping from local to global coordinates
                X, Y = self.map(s,t)
                nd = Node(X,Y)
                row.append(nd)
            nodes.append(row)

        for J in range(0, 2*NeY, 2):
            for I in range(0, 2*NeX, 2):
                node_list = [ nodes[J][I], nodes[J][I+2], nodes[J+2][I+2], nodes[J+2][I],      # corner nodes
                              nodes[J][I+1], nodes[J+1][I+2], nodes[J+2][I+1], nodes[J+1][I],  # midside nodes
                              nodes[J+1][I+1] ]                                                # center node
                elem = element_type(*node_list, deepcopy(material))
                elements.append(elem)

        nodes = [ nd for row in nodes for nd in row ]

        if self.model and isinstance(self.model, System):
            self.model.addNode(*nodes)
            self.model.addElement(*elements)

        self.nodes    = nodes
        self.elements = elements

        return (nodes, elements)

    def triangleMesh(self, NeX, NeY, element_type, material, **kwargs):
        r"""
        2D mesher using triangular elements

        This is an abstract class - requires implementation in subclass

        :param int NeX: number of elements along the first axis
        :param int NeY: number of elements along the second axis
        :param element_type: compatible element type
        :param material: a material object
        :returns (nodes, elements):
        :var nodes: list of created :py:class:`Node` objects
        :var elements: list of created (subclass of) :py:class:`Element` objects
        """

        match element_type.__name__:
            case 'Triangle':
                return self._triangle3_mesh(NeX, NeY, element_type, material, **kwargs)
            case 'Triangle6':
                return self._triangle6_mesh(NeX, NeY, element_type, material, **kwargs)
            case _:
                msg = "{} does not support element type {}".format(self.__class__.__name__, element_type.__name__)
                raise TypeError(msg)

    def _triangle3_mesh(self, NeX, NeY, element_type, material, **kwargs):

        # local coordinates
        ss = np.linspace(-1.,1.,NeX+1)
        tt = np.linspace(-1.,1.,NeY+1)

        nodes    = []
        elements = []

        for t in tt:
            row = []
            for s in ss:
                # mapping from local to global coordinates
                X, Y = self.map(s,t)
                nd = Node(X,Y)
                row.append(nd)
            nodes.append(row)

        for j in range(NeY):
            for i in range(NeX):
                node_list = [ nodes[j][i], nodes[j][i+1], nodes[j+1][i] ]
                elem = element_type(*node_list, deepcopy(material))
                elements.append(elem)
                node_list = [ nodes[j+1][i+1], nodes[j+1][i], nodes[j][i+1] ]
                elem = element_type(*node_list, deepcopy(material))
                elements.append(elem)

        nodes = [ nd for row in nodes for nd in row ]

        if self.model and isinstance(self.model, System):
            self.model.addNode(*nodes)
            self.model.addElement(*elements)

        self.nodes    = nodes
        self.elements = elements

        return (nodes, elements)

    def _triangle6_mesh(self, NeX, NeY, element_type, material, **kwargs):

        # local coordinates
        ss = np.linspace(-1., 1., 2*NeX + 1)
        tt = np.linspace(-1., 1., 2*NeY + 1)

        nodes = []
        elements = []

        for t in tt:
            row = []
            for s in ss:
                # mapping from local to global coordinates
                X, Y = self.map(s, t)
                nd = Node(X, Y)
                row.append(nd)
            nodes.append(row)

        for j in range(0, 2*NeY, 2):
            for i in range(0, 2*NeX, 2):
                node_list = [nodes[j][i], nodes[j][i+2], nodes[j+2][i],           # corner nodes
                             nodes[j][i+1], nodes[j+1][i+1], nodes[j+1][i]]       # midside nodes
                elem = element_type(*node_list, deepcopy(material))
                elements.append(elem)
                node_list = [nodes[j+2][i+2], nodes[j+2][i], nodes[j][i+2],       # corner nodes
                             nodes[j+2][i+1], nodes[j+1][i+1], nodes[j+1][i+2]]   # midside nodes
                elem = element_type(*node_list, deepcopy(material))
                elements.append(elem)

        nodes = [nd for row in nodes for nd in row]

        if self.model and isinstance(self.model, System):
            self.model.addNode(*nodes)
            self.model.addElement(*elements)

        self.nodes = nodes
        self.elements = elements

        return (nodes, elements)

    def map(self,s,t):
        """
        maps the local coordinates (s,t) from a bi-unit-square
        to the actual position (x,y) in the global model space
        """
        ss = [0.5*s*(s-1), (1-s)*(s+1), 0.5*s*(s+1)]
        tt = [0.5*t*(t-1), (1-t)*(t+1), 0.5*t*(t+1)]
        shp = [ss[0]*tt[0],
               ss[2]*tt[0],
               ss[2]*tt[2],
               ss[0]*tt[2],
               ss[1]*tt[0],
               ss[2]*tt[1],
               ss[1]*tt[2],
               ss[0]*tt[1],
               ss[1]*tt[1] ]
        x = shp @ self.X + self.offset

        return (x[0], x[1])


