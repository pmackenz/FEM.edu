from .Solver import *

class LinearSolver(Solver):
    r"""
    A linear system solver

    This solver implies :math:`\{{\bf P}\} = [{\bf K}] \{{\bf u}\}`
    """

    def __init__(self):
        super().__init__()

    def solve(self, **kwargs):
        r"""
        Solves the system assuming the given load is the total load
        and the obtained displacement is the **total** displacement.

        .. note::

           This method will not verify wether or not the linear assumption is correct.
           The resulting forces may be out of equilibrium if the system experiences
           nonlinear behavior under the given load.

        """
        self.resetDisplacements()
        self.assemble()
        errorNorm = self.solveSingleStep()
        # recover residual force vector
        self.assemble(force_only=True)

        # time to add the information to the recorded data
        if self.record:
            self.recordThisStep()

        return errorNorm

    def solveSingleStep(self):

        # are we doing displacement control?
        if self.hasConstraint:

            # solve for displacement update: a single Newton step
            dQ = np.linalg.solve(self.Kt, np.stack([self.R, self.P]).T)

            if self.useArcLength:
                # arc-length control
                delU  = self.sysU - self.lastConverged['U']
                dload = self.loadfactor - self.lastConverged['lambda']
                # g = self.arclength2 - delU@delU - self.alpha * dload*dload * self.P@self.P

                denum = 2.*np.dot(dQ[:,1], delU) + 2.*self.alpha * dload * self.P@self.P

                dlam = ( self.g - 2.*np.dot(dQ[:,0], delU) ) / denum

            else:
                # displacement control
                # g = (self.targetU - np.dot(self.sysU, self.en))

                dlam = self.g - np.dot(dQ[:,0], self.en)
                dlam /= np.dot(dQ[:,1], self.en)

            self.loadfactor += dlam
            dU = dQ[:,0] + dlam * dQ[:,1]

        else:
            # solve for displacement update: a single Newton step
            dU = np.linalg.solve(self.Kt, self.R)

        self.U = dU

        # update nodal displacements
        for node in self.nodes:
            idxK = node.getIdx4DOFs()
            node._updateDisp(dU[idxK])


    def assemble(self, force_only=False):
        r"""
        inherited from :code:`Solver` class.
        """
        super(LinearSolver, self).assemble(force_only=force_only)

    def pushState(self, state):
        r"""
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
        super(LinearSolver, self).pushState(state)

        if 'lam1' in state:
            self.loadfactor = state['lam1']
