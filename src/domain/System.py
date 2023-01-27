import numpy as np

from .Node     import *
from elements.Element import *
from plotter.Plotter import *


class System():
    """
    class: representing a System model
    """

    def __init__(self):
        """

        """
        self.nodes    = []
        self.elements = []
        self.plotter  = Plotter()
        self.disp     = np.array([])
        self.loads    = np.zeros_like(self.disp)

    def __str__(self):
        s = "System object"
        for node in self.nodes:
            s += "\n" + repr(node)
        for elem in self.elements:
            s += "\n" + repr(elem)
        return s

    def __repr__(self):
        return "System()"

    def addNode(self, *nodes):
        """

        :param newNode: a :ref:`Node` object
        """
        for newNode in nodes:
            if newNode not in self.nodes:
                newNode.index = len(self.nodes)
                self.nodes.append(newNode)
            else:
                print('addNode: node {} already exists in system and was not added again'.format(newNode.index))

    def __add__(self, other):
        if isinstance(other, Node):
            self.addNode(other)
        elif isinstance(other, Element):
            self.addElement(other)
        else:
            raise TypeError
        return self

    def addElement(self, newElement):
        """

        :param newElement: a :ref:`Element` object
        """
        self.elements.append(newElement)

    def solve(self):
        """

        """

        # compute size parameters
        ndof = 0
        for node in self.nodes:
            node.setStart(ndof)
            ndof += node.ndofs
        Rsys = np.zeros(ndof)
        Ksys = np.zeros((ndof, ndof))

        # assemble loads
        for node in self.nodes:
            if node.hasLoad():
                K = node.start + np.arange(node.ndofs)
                Rsys[K] += node.getLoad()

        # Element Loop: assemble element forces and stiffness
        for element in self.elements:
            Fe = element.getForce()     # Element State Update occurs here
            for (i,ndI) in enumerate(element.nodes):
                K = ndI.start + np.arange(ndI.ndofs)
                Rsys[K] -= Fe[i]
                for (j,ndJ) in enumerate(element.nodes):
                    M = ndJ.start + np.arange(ndJ.ndofs)
                    Ksys[K[:, np.newaxis], M] += element.Kt[i][j]

        # apply boundary conditions
        for node in self.nodes:
            for dof in node.dofs:
                if node.isFixed(dof):
                    idx = node.start + node.dofs[dof]
                    Rsys[idx]      = 0.0
                    Ksys[:, idx]   = np.zeros(ndof)
                    Ksys[idx, :]   = np.zeros(ndof)
                    Ksys[idx, idx] = 1.0

        # stability check for system matrix
        (vals, vecs) = np.linalg.eig(Ksys)
        for (lam, v) in zip(vals, vecs.T):
            if np.abs(lam) < 1.0e-2:
                print(f"lambda = {lam:16.12e}")
                print(v)

        # solve for displacements
        U = np.linalg.solve(Ksys, Rsys)

        # update nodal displacements
        for node in self.nodes:
            K = node.start + np.arange(node.ndofs)
            node._updateDisp(U[K])

        # recompute residual force
        Rsys = np.zeros(ndof)
        for element in self.elements:
            Fe = element.getForce()     # Element State Update occurs here
            for (i,ndI) in enumerate(element.nodes):
                K = ndI.start + np.arange(ndI.ndofs)
                Rsys[K] -= Fe[i]

        #
        self.Rsys = Rsys
        self.disp = U

    def plot(self, factor=1.0):
        """

        :param factor: deformation magnification factor
        """

        vertices = [ node.getPos() for node in self.nodes ]
        lines    = [ [ elem.nodes[k].index for k in [0,1] ] for elem in self.elements ]
        self.plotter.setMesh(vertices, lines)

        ndof = len(self.Rsys)
        R = self.Rsys.copy().reshape((ndof//2, 2))
        self.plotter.setReactions(R)

        disp = [ factor * node.getDisp() for node in self.nodes ]
        self.plotter.setDisplacements(disp)
        self.plotter.displacementPlot()

        values = [ elem.getAxialForce() for elem in self.elements ]
        self.plotter.setValues(values)
        self.plotter.valuePlot()

    def report(self):
        """
        print a text-based summary report

        """
        s  = "\nSystem Analysis Report\n"
        s += "=====================\n"
        s += "\nNodes:\n"
        s += "---------------------\n"
        for node in self.nodes:
            for ln in str(node).split('\n'):
                s += "  " + ln + "\n"
        s += "\nElements:\n"
        s += "---------------------\n"
        for elem in self.elements:
            for ln in str(elem).split('\n'):
                s += "  " + ln + "\n"
        print(s)

    def resetDisp(self):
        """
        Resets the displacement vector.
        """
        for node in self.nodes:
            node.resetDisp()

    def resetLoad(self):
        """
        Resets the load vector.
        """
        for node in self.nodes:
            node.resetLoad()

    def resetAll(self):
        """
        Resets load and displacement vectors.
        """
        self.resetDisp()
        self.resetLoad()


if __name__ == "__main__":

    from elements  import Element
    from materials import Material

    # testing the System class
    B = 6.0*12
    H = 8.0*12
    params = {'E':1000000., 'A':10., 'nu':0.0, 'fy':1.e30}

    model = System()

    nd0 = Node(0.0, 0.0)
    nd1 = Node(  B, 0.0)
    nd2 = Node(2*B, 0.0)
    nd3 = Node(3*B, 0.0)
    nd4 = Node(4*B, 0.0)
    nd5 = Node(0.5*B, H)
    nd6 = Node(1.5*B, H)
    nd7 = Node(2.5*B, H)
    nd8 = Node(3.5*B, H)

    model.addNode(nd0)
    model.addNode(nd1)
    model.addNode(nd2)
    model.addNode(nd3)
    model.addNode(nd4)
    model.addNode(nd5)
    model.addNode(nd6)
    model.addNode(nd7)
    model.addNode(nd8)

    model.addElement(Element(nd0, nd1, Material(params)))  # bottom 1
    model.addElement(Element(nd1, nd2, Material(params)))  # bottom 2
    model.addElement(Element(nd2, nd3, Material(params)))  # bottom 3
    model.addElement(Element(nd3, nd4, Material(params)))  # bottom 4

    model.addElement(Element(nd5, nd6, Material(params)))  # upper 1
    model.addElement(Element(nd6, nd7, Material(params)))  # upper 2
    model.addElement(Element(nd7, nd8, Material(params)))  # upper 3

    model.addElement(Element(nd0, nd5, Material(params)))  # up right diag 1
    model.addElement(Element(nd1, nd6, Material(params)))  # up right diag 2
    model.addElement(Element(nd2, nd7, Material(params)))  # up right diag 3
    model.addElement(Element(nd3, nd8, Material(params)))  # up right diag 4

    model.addElement(Element(nd1, nd5, Material(params)))  # up left diag 1
    model.addElement(Element(nd2, nd6, Material(params)))  # up left diag 2
    model.addElement(Element(nd3, nd7, Material(params)))  # up left diag 3
    model.addElement(Element(nd4, nd8, Material(params)))  # up left diag 4

    # boundary conditions
    nd0.fixDOF(0)
    nd0.fixDOF(1)
    nd4.fixDOF(1)

    # load upper nodes
    nd5.setLoad(0.0, -1.0)
    nd6.setLoad(0.0, -1.0)
    nd7.setLoad(0.0, -1.0)

    # solve the system
    model.solve()

    print(model)

    # write out report
    model.report()

