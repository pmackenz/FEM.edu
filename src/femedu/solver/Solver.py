import sys
import numpy as np
import scipy as sc

import matplotlib.pyplot as plt

class Solver():
    """
    Abstract class for any solver implementation.

    This class describes the functions needed by any solver
    """

    def __init__(self):
        """
        Initialize a solver instance with empty elements and nodes lists
        """
        self.loadfactor  = 1.0
        self.elements    = []
        self.nodes       = []
        self.constraints = []
        self.sdof = 0

        # numeric iteration tolerance
        self.TOL = 1.0e-6

        self.hasConstraint = False
        self.targetU       = 0.0

        self.useArcLength  = False
        self.lastConverged = {}

        # shall result be recorded?
        self.record = False

    def connect(self, nodes, elems, constraints):
        self.nodes       = nodes
        self.elements    = elems
        self.constraints = constraints

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
            if node.isLead():
                node.setStart(ndof)
                ndof += node.ndofs

        for constraint in self.constraints:
            constraint.setStart(ndof)
            ndof += constraint.countConditions()

        self.sdof = ndof  # number of system d.o.f.s

        Rsys = np.zeros(ndof)
        Ksys = np.zeros((ndof, ndof))

        # assemble loads
        for node in self.nodes:
            if node.isLead() and node.hasLoad():
                idx = node.start + np.arange(node.ndofs)
                Rsys[idx] += node.getLoad() * self.loadfactor

        # Element Loop: assemble element forces and stiffness
        for element in self.elements:
            Fe = element.getForce()     # Element State Update occurs here
            Pe = element.getLoad()      # Element State Update occurs here
            for (i,ndI) in enumerate(element.nodes):
                idxK = ndI.lead.start + ndI.getIdx4Element(element)
                if isinstance(Pe[i], np.ndarray):
                    Rsys[idxK] -= Fe[i] - self.loadfactor * Pe[i]
                else:
                    Rsys[idxK] -= Fe[i]
                if not force_only:
                    for (j,ndJ) in enumerate(element.nodes):
                        idxM = ndJ.lead.start + ndJ.getIdx4Element(element)
                        Ksys[idxK[:, np.newaxis], idxM] += element.Kt[i][j]

        # apply boundary conditions
        if not force_only:
            for node in self.nodes:
                for dof in node.dofs:
                    if node.isFixed(dof):
                        idx = node.lead.start + node.dofs[dof]
                        Rsys[idx]      = 0.0
                        Ksys[:, idx]   = np.zeros(ndof)   # the range might need adjustment for constraints
                        Ksys[idx, :]   = np.zeros(ndof)   # the range might need adjustment for constraints
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

    def on_converged(self):
        """
        This function needs to be called once a converged state was achieved by the solver.

        It tells all components to update its state to "converged"
        """
        for node in self.nodes:
            node.on_converged()

        for elem in self.elements:
            elem.on_converged()

        for const in self.constraints:
            const.on_converged()

    def revert(self):
        """
        This function needs to be called if the iterative procedure fails to converge.

        It will revert the entire system to the last converged state.
        """
        for node in self.nodes:
            node.revert()

        for elem in self.elements:
            elem.revert()

        for const in self.constraints:
            const.revert()

    def checkStability(self, verbose=True, **kwargs):
        """
        Computes the stability index as

        * :math:`\mathop{det}([{\\bf K}_t])` for systems with less than 25 d.o.f.s
        * :math:`\min\lambda_i` where :math:`\lambda_i` are the eigenvalues of :math:`{\\bf K}_t`

        :param verbose: set to **True** for log info
        :param num_eigen: if set to a value greater than 0, show the **num_eigen** eigenvalues
                        closest to the current load level.
        :returns: stability index
        """
        if 'num_eigen' in kwargs:
            num_eigen = kwargs['num_eigen']
        else:
            num_eigen = 0

        if self.sdof < 10 and not num_eigen:
            detKt = np.linalg.det(self.Kt)
            msg = f"\n ** Stability check: det(Kt) = {detKt}\n"
        else:
            evals = sc.linalg.eigvalsh(self.Kt)
            abs_evals = np.abs(evals)
            idx = np.argmin(abs_evals)
            detKt = evals[idx]  # we need to get the sign back
            if num_eigen > 0:
                detKt = []
                for cnt in range(num_eigen):
                    detKt.append(evals[idx])      # we need to get the sign back
                    abs_evals[idx] = 1.0e20
                    idx = np.argmin(abs_evals)
                detKt.sort()
                msg = f"\n ** Stability check: (smallest {num_eigen} eigenvalues of Kt)\n"
                for k, val in enumerate(detKt):
                    msg += f"\t\t\tmode {k}:{val:12.2f}\n"
            else:
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

            Use **System.plotBucklingMode()** for plotting of the eigenmode (see :doc:`../Domain/System_class`).
            This plotting method will be replaced by :code:`System.plot(...)` in a future release.

        :return: the eigenvalue, :math:`\lambda_{\mathtt{mode}}`
        """
        if not isinstance(mode,int) or mode < 0 or mode >= self.Kt.shape[0]:
            raise TypeError(f"mode out of range: must be an int between 0 and the number of d.o.f.s")

        w, v = sc.linalg.eigh(self.Kt, subset_by_index=[mode, mode])
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

    def showKt(self, filename="", **kwargs):

        plt.figure()
        plt.spy(self.Kt, marker='.',mec='b',mfc='b', **kwargs)
        if filename:
            plt.savefig(filename)
        else:
            plt.show()

    def getNodalReactions(self, dofs=None, cut_off=1.0e-6):

        self.assemble(force_only=True)

        R = []

        # assemble loads
        for node in self.nodes:

            reaction = np.zeros(3)

            if 'ux' in node.dofs:
                reaction[0] = self.R[node.lead.start + node.getIdx4DOFs(dofs=['ux'])]
            if 'uy' in node.dofs:
                reaction[1] = self.R[node.lead.start + node.getIdx4DOFs(dofs=['uy'])]
            if 'rz' in node.dofs:
                reaction[2] = self.R[node.lead.start + node.getIdx4DOFs(dofs=['rz'])]

            if np.linalg.norm(reaction) <= cut_off:
                reaction = np.zeros_like(reaction)

            R.append(reaction)

        return R


    def getNodalLoads(self, dofs=None, cut_off=1.0e-6):
        R = {}

        # collect nodal loads
        for node in self.nodes:
            R[node] = node.getLoad(dof_list=dofs, apply_load_factor=True)

        # Element Loop: assemble element forces and stiffness
        for element in self.elements:
            Pe = element.getLoad()      # Element State Update occurs here
            for (i,ndI) in enumerate(element.nodes):
                if isinstance(Pe[i], np.ndarray):
                    R[ndI] += self.loadfactor * Pe[i]

        Lds = []
        for node in self.nodes:
            force = R[node]
            if np.linalg.norm(force) <= cut_off:
                force = np.zeros_like(force)
            Lds.append(force)

        return Lds

