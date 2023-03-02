import sys
import numpy as np
import scipy as sp

class Solver():
    """
    Abstract class for any solver implementation.

    This class describes the functions needed by any solver
    """

    def __init__(self):
        """
        Initialize a solver instance with empty elements and nodes lists
        """
        self.loadfactor = 1.0
        self.elements = []
        self.nodes = []
        self.sdof = 0

        # numeric iteration tolerance
        self.TOL = 1.0e-6

        self.hasConstraint = False
        self.targetU       = 0.0

        self.useArcLength = False
        self.lastConverged = {}

        # shall result be recorded?
        self.record = False

    def connect(self, nodes, elems):
        self.nodes = nodes
        self.elements = elems

    def fetchState(self):
        """
        Fetch the current :code:`state` of the solver.

        .. list-table:: **state** is defined as a dictionary with the following contents:

            * - **nodes**
              - list of node pointers (required)
            * - **elements**
              - list of element pointers (required)
            * - **P0**
              - system vector of initial forces
            * - **Pref**
              - system vector of reference forces
            * - **u1**
              - system vector of current (converged) displacements
            * - **un**
              - system vector of previous (converged) displacements
            * - **lam1**
              - load level of current (converged) displacements
            * - **lamn**
              - load level of previous (converged) displacements


        :return: state of the solver
        """
        state = {}
        state['nodes']    = self.nodes
        state['elements'] = self.elements
        state['lam1']     = self.loadfactor

        return state

    def pushState(self, state):
        """
        Pushes :code:`state` to the solver.
        The solver will use that data to update it's internal state.

        .. list-table:: **state** is defined as a dictionary with the following contents:

            * - **nodes**
              - list of node pointers (required)
            * - **elements**
              - list of element pointers (required)
            * - **P0**
              - system vector of initial forces
            * - **Pref**
              - system vector of reference forces
            * - **u1**
              - system vector of current (converged) displacements
            * - **un**
              - system vector of previous (converged) displacements
            * - **lam1**
              - load level of current (converged) displacements
            * - **lamn**
              - load level of previous (converged) displacements

        :param state: state of the solver
        """
        if 'nodes' in state:
            self.nodes = state['nodes']
        else:
            raise TypeError("'nodes' missing from state")

        if 'elements' in state:
            self.elements = state['elements']
        else:
            raise TypeError("'elements' missing from state")

        if 'lam1' in state:
            self.loadfactor = state['lam1']
        else:
            raise TypeError("'lam1' missing from state")

    def setLoadFactor(self, lam):
        """
        Set the target load factor to **lam**

        The entered load pattern is considered a reference load,
        to be multiplied by a scalar load factor, :math:`\lambda`.

        If no load factor is set explicitly, a factor of 1.0 is assumed, i.e., the
        entire entered load is applied in full.
        """
        self.loadfactor = lam


    def assemble(self, force_only=False):
        """
        A general assembler for mixed element types.

        This method will build the out-of-balance force vector (residuum :math:`{\\bf R}`)
        and the tangent stiffness matrix (:math:`{\\bf K}_t`) used by most solvers.

        Specialized solvers may overload this method.

        .. note::

            The solver will apply the global load factor to the reference load returned by nodes
            and elements.  While nodes and elements are aware of that load factor, they shall
            apply it **only if asked explicity** by passing :code:`apply_load_factor=True` to the
            respective access functions.

        :param force_only: set to **True** if only the residual force needs to be assembled
        """

        # compute size parameters
        ndof = 0
        for node in self.nodes:
            node.setStart(ndof)
            ndof += node.ndofs
        Rsys = np.zeros(ndof)
        Ksys = np.zeros((ndof, ndof))

        self.sdof = ndof  # number of system d.o.f.s

        # assemble loads
        for node in self.nodes:
            if node.hasLoad():
                idx = node.start + np.arange(node.ndofs)
                Rsys[idx] += node.getLoad() * self.loadfactor

        # Element Loop: assemble element forces and stiffness
        for element in self.elements:
            Fe = element.getForce()     # Element State Update occurs here
            Pe = element.getLoad()      # Element State Update occurs here
            for (i,ndI) in enumerate(element.nodes):
                idxK = ndI.start + np.arange(ndI.ndofs)
                if isinstance(Pe[i], np.ndarray):
                    Rsys[idxK] -= Fe[i] - self.loadfactor * Pe[i]
                else:
                    Rsys[idxK] -= Fe[i]
                if not force_only:
                    for (j,ndJ) in enumerate(element.nodes):
                        idxM = ndJ.start + np.arange(ndJ.ndofs)
                        Ksys[idxK[:, np.newaxis], idxM] += element.Kt[i][j]

        # apply boundary conditions
        if not force_only:
            for node in self.nodes:
                for dof in node.dofs:
                    if node.isFixed(dof):
                        idx = node.start + node.dofs[dof]
                        Rsys[idx]      = 0.0
                        Ksys[:, idx]   = np.zeros(ndof)
                        Ksys[idx, :]   = np.zeros(ndof)
                        Ksys[idx, idx] = 1.0e3

            self.Kt = Ksys

        self.R  = Rsys


    def solve(self, **kwargs):
        """

        :return:
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def initialize(self):
        """

        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def reset(self):
        """
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def checkStability(self, verbose=True):
        """
        Computes the stability index as

        * :math:`\mathop{det}([{\\bf K}_t])` for systems with less than 25 d.o.f.s
        * :math:`\min\lambda_i` where :math:`\lambda_i` are the eigenvalues of :math:`{\\bf K}_t`

        :returns: stability index
        """
        if self.sdof < 5:
            detKt = np.linalg.det(self.Kt)
            msg = f"\n ** Stability check: det(Kt) = {detKt}\n"
        else:
            evals = sp.linalg.eigvals(self.Kt)
            idx = np.argmin(np.abs(evals))
            detKt = evals[idx]  # we need to get the sign back
            msg = f"\n ** Stability check: (smallest eigenvalue of Kt) = {detKt}\n"

        if verbose:
            print(msg)

        return detKt

    def getBucklingMode(self, mode=0, **kwargs):
        """
        Perform an eigen-analysis on :math:`{\\bf K}_t` for the requested **mode**.
        Default is the mode with the smallest absolute eigenvalue (:math:`\min\{\lambda_i\}`)

        The mode shape will be pushed to the nodes.

        .. note::

            Use **System.plotBucklingMode()** for plotting of the eigenmode (see :doc:`System_class`).
            This plotting method will be replaced by :code:`System.plot(...)` in a future release.

        :return: the eigenvalue, :math:`\lambda_{\mathtt{mode}}`
        """
        if not isinstance(mode,int) or mode < 0 or mode >= self.Kt.shape[0]:
            raise TypeError(f"mode out of range: must be an int between 0 and the number of d.o.f.s")

        w, v = sp.linalg.eigh(self.Kt, subset_by_index=[mode, mode])
        lam = w[0]
        U = v[:,0]

        # update nodal displacements
        for node in self.nodes:
            idxK = node.start + np.arange(node.ndofs)
            node.setDisp(U[idxK], modeshape=True)

        return lam

    def checkResiduum(self, report=False, force_only=True):

        # compute residual force and tangent stiffness
        self.assemble(force_only=force_only)

        normR = np.dot(self.R, self.R)

        # Add constraint violation in case we are using displacement control
        if self.hasConstraint:

            if self.useArcLength:
                # arc-length control
                delU  = self.sysU - self.lastConverged['U']
                dload = self.loadfactor - self.lastConverged['lambda']
                self.g = self.arclength2 - delU@delU - self.alpha * dload*dload * self.P@self.P

            else:
                # displacement control
                self.g = self.targetU - np.dot(self.sysU, self.en)  # this could be done more efficiently
                                                                    # since most values in self.en are zero
            normR += self.g*self.g

        else:
            self.g = 0.0

        normR = np.sqrt(normR)

        if report:
            #print('-')
            if self.hasConstraint:
                print(f"Constraint violation: {self.g:12.4e}")

            # R = self.getResiduum()
            # for i in range(len(R)):
            #     print("R({}): {:12.4e} {:12.4e}".format(i,*R[i]))

        print(f"norm of the out-of-balance force: {normR:12.4e}")

        # convergence check
        return normR

    def resetForces(self):
        """
        Reset force vector to **all zeros**.
        """
        self.P = np.zeros(self.sdof)

    def resetDisplacements(self):
        """
        Reset displacement vector to **all zeros**.
        """
        self.sysU = np.zeros(self.sdof)

