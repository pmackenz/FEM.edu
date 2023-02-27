import numpy as np

from .Node     import *
from ..elements.Element import *
from ..solver.LinearSolver import LinearSolver
from ..plotter.ElementPlotter import ElementPlotter as Plotter


class System():
    """
    class: representing a System model
    """

    def __init__(self):
        self.nodes    = []
        self.elements = []
        self.plotter  = Plotter()
        self.disp     = np.array([])
        self.loads    = np.zeros_like(self.disp)

        # global analysis settings
        self.loadfactor = 1.0

        self.initRecorder()
        self.trackStability(False)

        self.solver = LinearSolver()
        self.solver.connect(self.nodes, self.elements)

    def __str__(self):
        s = "System object"
        for node in self.nodes:
            s += "\n" + repr(node)
        for elem in self.elements:
            s += "\n" + repr(elem)
        return s

    def __repr__(self):
        return "System()"

    def setSolver(self, solver):
        """
        This method will change the current solver to the provided solver.

        Upon successful update of the solver, the old solver state will be
        pushed to the new solver.

        :param solver: a pointer to a solver object.
        """
        state = {}
        if self.solver:
            state = self.solver.fetchState()
        if solver:
            self.solver = solver
            self.solver.pushState(state)

    def addNode(self, *nodes):
        """
        Add one or more nodes to the model.

        :param newNode: a :code:`Node` object
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

        :param newElement: an :code:`Element` object
        """
        self.elements.append(newElement)

    def setLoadFactor(self, lam):
        self.solver.setLoadFactor(lam)

    def resetDisplacements(self):
        self.solver.resetDisplacements()

    def setConstraint(self, nodeIdx, nvec, Ubar):

        self.en = np.zeros(self.sdof)
        idx = np.arange(nodeIdx*self.ndof, (nodeIdx+1)*self.ndof)
        self.en[idx] = nvec
        self.targetU = Ubar

        self.hasConstraint = True
        self.useArcLength = False

    def getSolver(self):
        """
        Provides a pointer to the current solver instance.

        This function is used to get access to the solver and give instructions directly to that solver.
        """
        return self.solver

    def initArcLength(self, load=1., alpha=0.0, tolerance=1.0e-12):
        # store analysis parameter(s)
        self.alpha = alpha
        self.TOL = tolerance

        # make sure we start at an equilibrium point
        self.NewtonSolver()

        # store current configuration as last converged
        self.lastConverged = {'U':self.sysU.copy(), 'lambda':0.0, 'ds':0}

        # use load control to solve for the new equilibrium state
        self.hasConstraint = False   # this forces load control
        self.loadfactor += load      # add reference load level
        self.NewtonSolver()          # find equilibrium configuration for given load level

        # compute the arc-length for that step and store as target arc length
        delU  = self.sysU - self.lastConverged['U']
        dload = self.loadfactor - self.lastConverged['lambda']
        self.arclength2 = delU@delU + self.alpha * dload*dload * self.P@self.P

        # store the current point as last converged point
        self.previousConverged = self.lastConverged
        self.lastConverged = {'U':self.sysU.copy(), 'lambda':0.0, 'ds':np.sqrt(self.arclength2)}

        # set solver parameters
        self.hasConstraint = True
        self.useArcLength = True

    def stepArcLength(self, verbose=False):

        if not self.hasConstraint:
            # this method makes no sense for load control
            return

        # set solver parameters
        self.useArcLength = True

        # store current state as local variable
        Un = {'U':self.sysU.copy(), 'lambda':0.0, 'ds':np.sqrt(self.arclength2)}

        # set suitable trial state
        self.sysU *= 2.0
        self.sysU -= self.previousConverged['U']
        self.loadfactor *= 2.0
        self.loadfactor -= self.previousConverged['lambda']

        # store current state as "last converged solution"
        self.previousConverged = self.lastConverged
        self.lastConverged = Un

        if verbose:
            print(self.lastConverged)

        # solve for next point on the equilibrium path
        self.NewtonSolver(verbose)



# -------- recorder methods -------------

    def initRecorder(self):
        self.record = False
        self.recorder = {'U':[], 'lambda':[], 'stability':None}

    def startRecorder(self):
        """
        starts the recorder
        """
        self.record = True

    def stopRecorder(self):
        """
        stops the recorder
        """
        self.record = False

    def recordThisStep(self):
        """
        record current state of the system
        """
        self.recorder['U'].append(self.sysU.copy())
        self.recorder['lambda'].append(self.loadfactor)
        if self.track_stability:
            self.recorder['stability'].append(self.checkStability())
    def trackStability(self, on=True):
        if on:
            self.track_stability = True
            if not self.recorder['stability']:
                self.recorder['stability'] = [ 0.0 for x in self.recorder['lambda'] ]
        else:
            self.track_stability = False

    def fetchRecord(self):
        """
        :return: a tuple of ndarrays: (loadfactors, list of system deformations, stability index)
        """
        return (np.array(self.recorder['lambda']), np.array(self.recorder['U']), np.array(self.recorder['stability']))

# ------- solver methods ----------

    def solve(self, **kwargs):
        """
        Solve system of equations and find state of deformation for the given load level.
        """
        if self.solve:
            self.solver.solve()
        else:
            msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
            raise NotImplementedError(msg)

# ------------ plot methods -----------------

    def plot(self, factor=1.0, filename=None):
        """
        Create mesh plot showing the undeformed and the deformed system

        If **filename** is given, store the plot to that file.
        Use proper file extensions to indicate the desired format (.png, .pdf)

        :param factor: deformation magnification factor
        :param filename:  filename (str)
        """

        self.plotter.setMesh(self.nodes, self.elements)

        # ndof = len(self.Rsys)
        # R = self.Rsys.copy().reshape((ndof//2, 2))
        # self.plotter.setReactions(R)

        self.plotter.displacementPlot(factor=factor, file=filename)

    def valuePlot(self, variable, factor=0.0, filename=None):
        """
        Create a false color contour plot for the selected variable.
        A value of zero (0.0) will be assigned for any variable not
        provided by an element or a node.

        If **filename** is given, store the plot to that file.
        Use proper file extensions to indicate the desired format (.png, .pdf)

        :param variable: string code for variable to show
        :param factor: deformation magnification factor (default is undeformed)
        :param filename: filename (str)
        """

        self.plotter.setMesh(self.nodes, self.elements)
        self.plotter.valuePlot(variable_name=variable, factor=factor, file=filename)

    def beamValuePlot(self, variable, factor=0.0, filename=None):
        """
        Create a traditional beam value plot, i.e., moment and shear diagrams.

        If **filename** is given, store the plot to that file.
        Use proper file extensions to indicate the desired format (.png, .pdf)

        :param variable: string code for variable
        :param deformed: True | **False**
        :param file: filename (str)
        """

        self.plotter.setMesh(self.nodes, self.elements)
        self.plotter.beamValuePlot(variable_name=variable, factor=factor, file=filename)


    def report(self):
        """
        print a text-based summary report

        """
        s  = "\nSystem Analysis Report\n"
        s += "=======================\n"
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

# ------------ operational support methods --------------

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

    from ..elements  import Element
    from ..materials import Material

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

