from multiprocessing.connection import answer_challenge
from xmlrpc.client import INVALID_ENCODING_CHAR

import numpy as np

###
### matplotlib is needed for deprecated methods only
###
### remove this section onve ALL plotting has been moved into suitable plotter classes
###
import matplotlib.pyplot as plt
from fontTools.misc.textTools import caselessSort

###
### end of deprecated section
###

from .Node import Node
from ..elements.Element import *
from ..solver.LinearSolver import LinearSolver
from ..plotter.ElementPlotter import ElementPlotter as Plotter
from ..recorder import *


class System():
    """
    class: representing a System model
    """

    def __init__(self, verbose=False):
        self.nodes       = []
        self.elements    = []
        self.constraints = []
        self.plotter     = Plotter()

        self.verbose = verbose

        self.disp        = np.array([])
        self.loads       = np.zeros_like(self.disp)

        # global analysis settings
        self.loadfactor = 1.0

        self.solver = LinearSolver()
        self.solver.connect(self, self.nodes, self.elements, self.constraints)

        self.initRecorder()
        self.trackStability(False)

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
                newNode.setLoadFactor(self.loadfactor)
                self.nodes.append(newNode)
            elif self.verbose:
                print('addNode: node {} already exists in system and was not added again'.format(newNode.getID()))

    def __add__(self, other):
        if isinstance(other, Node):
            self.addNode(other)
        elif isinstance(other, Element):
            self.addElement(other)
        else:
            raise TypeError
        return self

    def addElement(self, *newElements):
        """

        :param newElement: an :code:`Element` object
        """
        for elem in newElements:
            elem.setLoadFactor(self.loadfactor)
            self.elements.append(elem)

    def addConstraint(self, *newConstraints):
        """

        :param newConstraints: one or more :py:class:`Constraint` or a subclass objects
        """
        for constraint in newConstraints:
            self.constraints.append(constraint)

# --------- load control functions ----------------------

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
        r"""
        Provides a pointer to the current solver instance.

        This function is used to get access to the solver and give instructions directly to that solver.
        """
        return self.solver

# --------- modeling assist functions ----------------

    def _is_node_on_line(self, node, pos, dir, tol=1.e-3):
        tvec = np.array(dir,dtype=np.float64)
        tvec /= np.linalg.norm(tvec)  # normalized tangent vector

        R = node.getPos() - np.array(pos)
        distX = tvec @ R
        distR = R - distX * tvec
        dist = np.linalg.norm(distR)

        return (dist < tol)


    def findNodesAt(self, pos, tol=1.e-3):
        r"""
        Find any node within distance tol from the given target position (**pos**).

        If a 2d position is given for a 3d-problem, the input target defines a cylinder with a vertical axis (3rd direction) and radius **tol**.

        If a 3d position is given for a 2d-problem, the 3rd coordinate of **pos** will be ignored.

        :param pos: a coordinate tuple (2d or 3d) defining th etarget position
        :param tol: search radius for find. Any node within distance **tol** will be returned.
        :return: a list of tuples (node_pointer, distance)
        """
        ans = []

        for node in self.nodes:
            dist = node.distanceTo(pos)
            if dist<tol:
                ans.append((node,dist))

        return ans

    def findNodesAlongLine(self, pos, dir, tol=1.e-3):
        r"""
        Find any node within distance tol from the given target position (**pos**).

        If a 2d position is given for a 3d-problem, the input target defines a cylinder with a vertical axis (3rd direction) and radius **tol**.

        If a 3d position is given for a 2d-problem, the 3rd coordinate of **pos** will be ignored.

        :param pos: a coordinate tuple (2d or 3d) defining th etarget position
        :param dir: direction vector defining the orientation of the line
        :param tol: search radius for find. Any node within distance **tol** will be returned.
        :return: a list of tuples (node_pointer, distance)
        """
        ans = []

        for node in self.nodes:
            if self._is_node_on_line(node, pos, dir, tol=tol):
                ans.append((node, node.distanceTo(pos)))

        return ans

    def findNodesOnPlane(self, pos, dir, tol=1.e-3):
        r"""
        Find any node within distance tol from the reference plane.

        The plane is defined by

        .. math::
            | {\bf n} \cdot ( {\bf X}_{node} - {\bf X}_0 ) | < TOL

        where
        :math:`{\bf n}` is the unit normal vector to the plane (normalized **dir**),
        :math:`{\bf X}_{node}` is the position vector of any node, and
        :math:`{\bf X}_0 )` is the position vector of a reference point on the plane (**pos**),

        :param pos: a coordinate tuple (2d or 3d) defining th etarget position
        :param dir: normal vector to the target plane
        :param tol: search radius for find. Any node within distance **tol** will be returned.
        :return: a list of tuples (node_pointer, distance)
        """
        X0 = np.array(pos)             # reference point
        nvec = np.array(dir)
        nvec /= np.linalg.norm(nvec)   # normalized normal vector

        ans = []

        for node in self.nodes:
            dist = np.abs( nvec @ (node.getPos() - X0) )
            if dist<tol:
                ans.append((node,dist))

        return ans

    def findFacesAlongLine(self, pos, dir, orientation=0, tol=1.e-2):
        r"""
        Find any element faces within distance tol from the given target position (**pos**).
        If a 2d position is given for a 3d-problem, the input target defines a cylinder with a vertical axis (3rd direction) and radius **tol**.
        If a 3d position is given for a 2d-problem, the 3rd coordinate of **pos** will be ignored.


        .. list-table::
            :header-rows: 1

            * - orientation
              - filter
            * - `0`
              - Any face aligned with the plane will be returned
            * - `+1`
              - Only faces with an outward normal pointing to the right of **dir** will be selected.
            * - `-1`
              - Only faces with an outward normal pointing to the left of **dir** will be selected.



        :param pos: a coordinate tuple (2d or 3d) defining th etarget position
        :param dir: direction vector defining the orientation of the line
        :param orientation: search filter
        :param tol: search radius for find. Any node within distance **tol** will be returned.
        :return: a list of tuples (node_pointer, distance)
        """
        ans = []

        Xo = np.array(pos, dtype=np.float64)
        To = np.array(dir, dtype=np.float64)
        To /= np.linalg.norm(To)

        match len(dir):
            case 2:
                No = np.array([[ 0.0, 1.0],[ -1.0, 0.0]]) @ To  # rotate to the right
            case 3:
                msg = "findFacesAlongLine requires a 2d-analysis.  Use findEdgesAlongLine() for shells and findFacesOnPlane() for 3d solids."
                raise TypeError(msg)

        boundary_nodes = self.findNodesAlongLine(pos, dir, tol=tol)

        for node, dnode in boundary_nodes:
            for elem in node.elements:
                #print('+', elem)
                for face in elem.faces:
                    for x, area in zip(face.pos, face.area):
                        match orientation:
                            case 0:   # both sides
                                dist      = np.abs((x - Xo) @ No)
                                area_test = np.abs( No @ area / np.linalg.norm(area) )
                            case -1:  # left side only
                                dist      = -(x - Xo) @ No
                                area_test = -No @ area / np.linalg.norm(area)
                            case 1:   # right side only
                                dist      = (x - Xo) @ No
                                area_test = No @ area / np.linalg.norm(area)

                        if  dist < tol and area_test > (1. - 10. * tol):
                            # avoid duplication
                            if (elem,face) not in ans:
                                ans.append((elem,face))

        return ans

    def findFacesOnPlane(self, pos, dir, orientation=0, tol=1.e-3):
        r"""
        Find any element faces within distance tol from the given target position (**pos**).

        If a 2d position is given for a 3d-problem, the input target defines a cylinder with a vertical axis (3rd direction) and radius **tol**.

        If a 3d position is given for a 2d-problem, the 3rd coordinate of **pos** will be ignored.


        .. list-table::
            :header-rows: 1

            * - orientation
              - filter
            * - `0`
              - Any face aligned with the plane will be returned
            * - `+1`
              - Only faces with an outward normal pointing in the direction of **dir** will be selected.
            * - `-1`
              - Only faces with an outward normal pointing against the direction of **dir** will be selected.


        :param pos: a coordinate tuple (2d or 3d) defining th etarget position
        :param dir: normal vector to the target plane
        :param orientation: search filter
        :param tol: search radius for find. Any node within distance **tol** will be returned.
        :return: a list of tuples (element_ptr, face_ptr)
        """
        ans = []
        return ans

# --------- displacement control functions ----------------
    def setDisplacementControl(self, node, dof, target):
        r"""
        activate displacement control for the next load step

        :param node:     pointer to the controlling node
        :param dof:      dof code for controlled dof
        :param target:   target displacement value
        """
        if self.solver:
            self.solver.setDisplacementControl(node, dof, target)

# --------- Arc-length control functions: forward to Solver ------------

    def initArcLength(self, load_increment=1., alpha=0.0, tolerance=1.0e-12):
        r"""
        Initializes parameters for the arc-length constraint.

        .. math::

           g({\bf u}, \lambda)) := \alpha ||\bar {\bf P}|| (\lambda-\lambda_n)^2 + ({\bf u} - {\bf u}_n)({\bf u} - {\bf u}_n) - \Delta s^2 = 0

        .. note::

           This feature requires a nonlinear solver. Review the :py:meth:`setSolver` function.

        :param load_increment:   load increment used to calibrate the constraint
        :param alpha:            load contribution factor
        :param tolerance:        convergence tolerance
        """
        if self.solver:
            self.solver.initArcLength(load_increment=load_increment, alpha=alpha, tolerance=tolerance)

    def stepArcLength(self, verbose=False, max_iter=10):
        r"""
        Progresses the model state by one arc-length.

        .. note::

            You need to initialize arc-length control by one call to
            :py:meth:`initArcLength` at least once to set all necessary parameters.

        :return normR: the norm of the generalized residuum from the last iteration step
        """

        if self.solver:
            loadfactor, normR = self.solver.stepArcLength(verbose=verbose, max_iter=max_iter)
            self.setLoadFactor(loadfactor)

    # --------- recorder methods ------------------------------

    def initRecorder(self, **kwargs):
        r"""
        initializes data arrays for gathering of load history data

        :keyword variables: list of variables or d.o.f.-codes to be recorded
        :type variables: list-type
        :keyword nodes:  nodes for which to record
        :type nodes: list of :py:class:`Node`
        :keyword elements:  elements for which to record
        :type elements: list of :py:class:`Element`

        .. code::

            # examples:

            # system recorder tracking the stability index
            initRecorder(variables=['stability',])

            # tracking 'ux' at nodes X1 and X7
            initRecorder(variables=['ux',], nodes=[X1, X7])

            # tracking all components of stress at Elem2 and Elem42
            initRecorder(variables=['stress',], elements=[Elem2, Elem42, ...])

        """
        if 'variables' in kwargs:
            if 'lam' not in kwargs['variables']:
                kwargs['variables'] = list(kwargs['variables']) + ['lam']
            #if 'stability' not in kwargs['variables']:
            #    kwargs['variables'] = list(kwargs['variables']) + ['stability']
        else:
            kwargs['variables'] = ['lam','stability']

        if 'nodes' in kwargs and isinstance(kwargs['nodes'],str):
            if kwargs['nodes'].lower() == 'all':
                kwargs['nodes'] = self.nodes
            else:
                msg = f"Unknown node identifier `nodes={kwargs['nodes']} in Recorder`"
                raise KeyError(msg)

        if 'elements' in kwargs and isinstance(kwargs['elements'],str):
            if kwargs['elements'].lower == 'all':
                kwargs['elements'] = self.elements
            else:
                msg = f"Unknown node identifier `elements={kwargs['elements']} in Recorder`"
                raise KeyError(msg)

        self.recorder = ModelRecorder(**kwargs)

    def startRecorder(self):
        """
        starts the recorder
        """
        if self.recorder:
            self.recorder.enable()
        if self.solver:
            self.solver.startRecorder()
        for node in self.nodes:
            node.startRecorder()
        for elem in self.elements:
            elem.startRecorder()

    def pauseRecorder(self):
        """
        pauses the recorder.  You can restart the recorder by calling
        `startRecorder()`.
        """
        if self.recorder:
            self.recorder.disable()
        if self.solver:
            self.solver.pauseRecorder()
        for node in self.nodes:
            node.pauseRecorder()
        for elem in self.elements:
            elem.pauseRecorder()

    def stopRecorder(self):
        """
        stops the recorder
        """
        if self.recorder:
            self.recorder.disable()
        if self.solver:
            self.solver.stopRecorder()
        for node in self.nodes:
            node.stopRecorder()
        for elem in self.elements:
            elem.startRecorderstopRecorder()

    def recordThisStep(self):
        """
        record current state of the system
        """

        for node in self.nodes:
            node.recordThisStep(self.loadfactor)

        for elem in self.elements:
            elem.recordThisStep(self.loadfactor)

        if self.recorder and self.recorder.isActive():
            data = {'lam':self.loadfactor}
            if self.track_stability:
                data['stability'] = self.solver.checkStability(num_eigen=1)
            else:
                data['stability'] = np.nan

            self.recorder.addData(data)

    def trackStability(self, on=True):
        self.track_stability = on

    def fetchRecord(self, keys=[]):
        """
        Request recorded time history data for the listed keys.
        If a single key is given as a string, a single `np.array()` is returned.
        If a list of keys is given, a `list of np.array()` is returned.

        :returns: time history data for the listed keys.
        """
        if keys:
            return self.recorder.fetchRecord(keys=keys)
        else:
            return []

# ------- solver methods ----------

    def solve(self, **kwargs):
        """
        Solve system of equations and find state of deformation for the given load level.
        """
        if self.solver:
            self.solver.solve(**kwargs)
            if self.solver.hasConstraint:
                # spread the news about the new load level throughout the system
                self.setLoadFactor(self.solver.loadfactor)
        else:
            msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
            raise NotImplementedError(msg)

    def checkStability(self, **kwdargs):
        r"""
        Computes the stability index as

        * :math:`\mathop{det}([{\bf K}_t])` for systems with less than 25 d.o.f.s
        * :math:`\min\lambda_i` where :math:`\lambda_i` are the eigenvalues of :math:`{\bf K}_t`

        **Implemented** by :py:class:`Solver`.

        :param verbose: set to **True** for log info
        :param num_eigen: if set to a value greater than 0, show the **num_eigen** eigenvalues
                        closest to the current load level.
        :returns: stability index
        """
        if self.solver:
            return self.solver.checkStability(**kwdargs)

    def getBucklingMode(self, mode=0, **kwargs):
        if self.solver:
            return self.solver.getBucklingMode(mode=mode, **kwargs)
        else:
            return (np.nan,None)

# ------------ plot methods -----------------

    def plot(self, factor=1.0, show_reactions=True, show_loads=True, force_limit=1.0e-6, **kwargs):
        """
        Create mesh plot showing the undeformed and the deformed system

        If **filename** is given, store the plot to that file.
        Use proper file extensions to indicate the desired format (.png, .pdf)

        :param factor: deformation magnification factor
        :param show_reactions: plot reaction forces if this evaluates to **True**
        :param show_loads:  plot nodal load vectors if this evaluates to **True**
        :param force_limit: don't plot forces smaller than this value
        :param filename:  filename (str)
        """

        self.plotter.setMesh(self.nodes, self.elements)

        if show_loads:

            kwargs['show_loads'] = True
            R = self.solver.getNodalLoads(dofs=('ux','uy','rz'), cut_off=force_limit)
            self.plotter.setNodalLoads(R)

        if show_reactions:

            kwargs['show_reactions'] = True
            R = self.solver.getNodalReactions(dofs=('ux','uy','rz'), cut_off=force_limit)
            self.plotter.setReactions(R)

        self.plotter.displacementPlot(factor=factor, **kwargs)

    def valuePlot(self, variable, factor=0.0, filename=None, **kwargs):
        r"""
        Create a false color contour plot for the selected variable.
        A value of zero (0.0) will be assigned for any variable not
        provided by an element or a node.

        .. list-table:: known variables
            :header-rows: 1

            * - keyword
              - description
            * - 'ux', 'uy', 'uz'
              - displacement component in the `x`, `y` or `z` direction
            * - 'rx', 'ry', 'rz'
              - rotation around the `x`, `y` or `z` axis
            * - 'sxx', 'syy', 'szz'
              - normal stress :math:`\sigma_{xx}`, :math:`\sigma_{yy}`, :math:`\sigma_{zz}`
            * - 'sxy', 'syz', 'szx'
              - shear stress :math:`\sigma_{xy}`, :math:`\sigma_{yz}`, :math:`\sigma_{zx}`
            * - 'epsxx', 'epsyy', 'epszz'
              - normal strain :math:`\varepsilon_{xx}`, :math:`\varepsilon_{yy}`, :math:`\varepsilon_{zz}`
            * - 'epsxy', 'epsyz', 'epszx'
              - engineering shear strain :math:`\gamma_{xy}=2\varepsilon_{xy}`,
                :math:`\gamma_{yz}=2\varepsilon_{yz}`,
                :math:`\gamma_{zx}=2\varepsilon_{zx}`
            * - 'nxx', 'nyy', 'nxy'
              - membrane forces :math:`n_{xx}=h\sigma_{xx}`, :math:`n_{yy}=h\sigma_{yy}`, :math:`n_{xy}=h\sigma_{xy}`
            * - 'T'
              - temperature
            * - 'qx', 'qy', 'qz'
              - `x`, `y` or `z` component of the temperature gradient, :math:`q_i = \partial T/\partial x_i`


        :param variable: string code for variable to show
        :param factor: deformation magnification factor (default is undeformed)
        :param filename: if **filename** is given, store the plot to that file.
                Use proper file extensions to indicate the desired format (.png, .pdf)
        :type filename: str
        """
        varIsAtGaussPoint = True

        for node in self.nodes:
            if node.hasDOF(variable):
                varIsAtGaussPoint = False
                break

        if varIsAtGaussPoint:
            for node in self.nodes:
                node._resetGaussPointMap(variable)
            for elem in self.elements:
                elem.mapGaussPoints(variable)

        self.plotter.setMesh(self.nodes, self.elements)
        self.plotter.valuePlot(variable_name=variable,
                               factor=factor,
                               gaussptvalue=varIsAtGaussPoint,
                               filename=filename, **kwargs)

    def beamValuePlot(self, variable, factor=0.0, filename=None, **kwargs):
        """
        Create a traditional beam value plot, i.e., moment and shear diagrams.


        :param variable: string code for variable
        :param deformed: True | **False**
        :param filename: if **filename** is given, store the plot to that file.
                Use proper file extensions to indicate the desired format (.png, .pdf)
        :type filename: str
        """

        self.plotter.setMesh(self.nodes, self.elements)
        self.plotter.beamValuePlot(variable_name=variable, factor=factor, filename=filename, **kwargs)

    def historyPlot(self, varX='', varY='', filename=None, **kwargs):
        """
        Create a generic X-Y plot using recorder data for load-level (horizontal)
        and varY (vertical).

        If **filename** is given, store the plot to that file.
        Use proper file extensions to indicate the desired format (.png, .pdf)

        :param varX:     X-axis parameter;
                         This must be a variable code previously set by :py:meth:`initRecorder`.

                         If :py:attr:`varX` is a :py:class:`str`, it is assumed to be a global parameter/variable, e.g.,
                         load factor **lam**, or pseudo time **time** (time not yet supported)

                         If :py:attr:`varX` is a :py:class:`tuple` containing a variable code and a :py:class:`Node` pointer,
                         the respective nodal DOF defined by the parameter will be used.

                         If :py:attr:`varX` is a :py:class:`tuple` containing a variable code and a :py:class:`Element` pointer,
                         the respective element DOF defined by the parameter will be used. The parameter must identify a scalar parameter.
                         Hence, you may use :py:data:`'sig:xx'` but **not** :py:data:`'stress'`.

        :param varY:     a variable code previously set by :py:meth:`initRecorder`
        :param filename:
        :type filename: str
        :param node: (optional)
        :param nodes: (optional)
        """
        if not self.recorder:
            return

        if not varX:
            variableX = 'lam'
            src = [None]
        elif isinstance(varX, list) or isinstance(varX, tuple):
            variableX = varX[0]
            src = [varX[1]]
        else:
            variableX = varX
            src = [None]

        if isinstance(varY, list) or isinstance(varY, tuple):
            variableY = varY
        else:
            variableY = [varY]

        if 'nodes' in kwargs:
            data = self.recorder.fetchRecord([variableX, *variableY], source=src + kwargs['nodes'])
        elif 'node' in kwargs:
            data = self.recorder.fetchRecord([variableX, *variableY], source=src + [kwargs['node']])
        elif src:
            data = self.recorder.fetchRecord([variableX, *variableY], source=src + [None]*len(variableY))
        else:
            data = self.recorder.fetchRecord([variableX, *variableY])


        X = Record()
        for key, record in data.items():
            if record.key == variableX:
                X = record

        Y = []
        for key, record in data.items():
            if record.key in variableY:
                Y.append( record )

        if 'title' not in kwargs:
            kwargs['title'] = "Load History Plot"

        self.plotter.xyPlot(X, Y, filename=filename, **kwargs)

    def plotDOF(self, dofs=[]):
        """
        .. warning::

            This method has been marked **DEPRECATED** !

            Use the :code:`System.plot()` method with appropriate :code:`**kwargs` instead.
        """


        if self.track_stability:
            fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,4))
        else:
            fig, ax1 = plt.subplots(1,1, figsize=(5,4))

        if not dofs:
            dofs = range(self.sdof)

        loadfactors   = self.recorder['lambda']
        displacements = np.array(self.recorder['U'])

        for idx in dofs:

            nodeID    = idx // self.ndof + 1
            direction = idx % self.ndof

            if direction==0:
                style = '--x'
                lbl = '$ u_{{{}}} $'.format(nodeID)
            elif direction==1:
                style = '-o'
                lbl = '$ v_{{{}}} $'.format(nodeID)
            elif direction==2 and direction < self.ndim:
                style = ':+'
                lbl = '$ w_{{{}}} $'.format(nodeID)
            else:
                if direction < self.ndim:
                    style = '-+'
                    lbl = '$ u_{{{},{}}} $'.format(nodeID,direction)
                else:
                    style = '-.o'
                    if self.ndim < 3:
                        lbl = '$ \\theta_{{{}}} $'.format(nodeID)
                    else:
                        lbl = '$ \\theta_{{{},{}}} $'.format(nodeID,direction-self.ndim)

            ax1.plot(displacements[:,idx],loadfactors, style,label=lbl)

        ax1.grid(True)
        ax1.set_xlabel('component of displacement, $ u_i$')
        ax1.set_ylabel('load factor, $\\lambda$')
        ax1.legend()

        if self.track_stability:

            detKt = np.array(self.recorder['stability'])

            ax2.plot(detKt, loadfactors,'-r')
            ax2.plot(np.zeros_like(loadfactors), loadfactors,'-k', lw=2)
            ax2.grid(True)
            ax2.set_xlabel('stability index, $ det({\\bf K}_t) $ or $ \\min |\\lambda| $')
            ax2.set_ylabel('load factor, $\\lambda$')

        plt.savefig("history_plots.png", bbox_inches='tight')
        plt.show()

    def plotSystem(self, filename=""):
        r"""
        This is a shortcut for

        .. code::

            System.plot(filename=filename,
                        factor=0.0,
                        title="Undeformed system",
                        show_loads=1,
                        show_reactions=0)


        :param filename: if given, save the plot to that file
        """
        self.plot(filename=filename, factor=0.0, title="Undeformed system", show_loads=1, show_reactions=0)

    def pushU(self):
        """
        Store the current displacement vector for later restore using :code:`popU()`.
        """
        for node in self.nodes:
            node.pushU()

    def popU(self):
        """
        Restore a previously pushed displacement vector (using :code:`pushU()`).
        """
        for node in self.nodes:
            node.popU()

    def plotBucklingMode(self, factor=1.0, mode=0, filename=None, **kwargs):
        """
        .. warning::

            This method has been marked **DEPRECATED** !

            Use the :code:`System.plot()` method with appropriate :code:`**kwargs` instead.

        Select a mode shape by setting `mode` to an integer, where `0` stands for the mode with the smallest
        absolute eigenvalue, `1` the second smallest, and so forth.  If the resulting eigenvalue (shown in the plot)
        is negative, the load-level associated with that mode is below the current load level.

        If a `filename` is given, e.g., `filename="buckling_mode.png"`, the plot will be saved to the respective file.
        The file type will be determined from the given extensions, e.g., `.png` for a PNG format.

        :param factor: Scale the mode shape by this factor.
        :param mode:   which mode shape shall be plotted?
        :param filename: string
        :param kwargs: not specified here.  These are handed over to the :code:`Solver(...)` instance.
        """

        # perform eigenvalue analysis
        lam = self.getBucklingMode(mode=mode, **kwargs)
        title_text = f"Mode Shape for $ \\lambda = {lam:.2f} $"
        self.plot(factor=factor, filename=filename, title=title_text, modeshape=True)

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

    def setLoadFactor(self, lam):
        r"""
        Set the target load factor to **lam**

        The entered load pattern is considered a reference load,
        to be multiplied by a scalar load factor, :math:`\lambda`.

        If no load factor is set explicitly, a factor of 1.0 is assumed, i.e., the
        entire entered load is applied in full.
        """
        self.loadfactor = lam

        # inform elements about the load factor
        for elem in self.elements:
            elem.setLoadFactor(lam)

        # inform nodes about the new load factor
        for node in self.nodes:
            node.setLoadFactor(lam)

        # let solver know about the new load factor
        if self.solver:
            self.solver.setLoadFactor(lam)

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

    msg = """
    This file is not an executable script.
    See the FEM.edu manual for details on how to run the program.
    """
    print(msg)
