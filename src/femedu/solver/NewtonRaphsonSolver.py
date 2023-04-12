import numpy as np

from ..solver.Solver import Solver

class NewtonRaphsonSolver(Solver):
    """
    An iterative solver for load-controlled nonlinear analysis.
    """

    def __init__(self):
        super(NewtonRaphsonSolver, self).__init__()

    def solve(self, max_steps=10, verbose=False, **kwargs):
        """
        :param max_step: maximum number of iterations (int)
        :param verbose: set to :code:`True` for additional information
        """
        TOL = self.TOL

        if 'tol' in kwargs:
            TOL = kwargs['tol']
        if 'tolerance' in kwargs:
            TOL = kwargs['tolerance']

        for k in range(max_steps):

        # compute force vector and tangent stiffness
            self.assemble()

            # we are using , force_only=False and simply reuse self.Kt from that run (!)
            normR = self.checkResiduum(verbose, force_only=False)

            if normR < TOL:
                # we achieved convergence
                #
                # now broadcast that we can switch to this converged state
                self.on_converged()

                # time to add the information to the recorded data
                if self.record:
                    self.recordThisStep()

                break

            # Solve for equilibrium
            self.solveSingleStep()

        print('+')

        if normR > TOL:
            # we failed to converge
            #
            # back to safety: revert to the last converged step
            self.revert()

        return normR

    def solveSingleStep(self):
        """
        Helper function performing a single solution of the linearized system

        Called by **solve()**. (internal use only)
        """

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

        # update nodal displacements
        for node in self.nodes:
            idxK = node.lead.start + np.arange(node.ndofs)
            node._updateDisp(dU[idxK])

    def assemble(self, force_only=False):
        """
        inherited from :code:`Solver` class.
        """
        super(NewtonRaphsonSolver, self).assemble(force_only=force_only)

    def getDisplacements(self):
        U = self.sysU.copy()
        U.shape = (self.nNodes, self.ndof)
        return U

    def getForces(self):
        P = self.P.copy() * self.loadfactor
        P.shape = (self.nNodes, self.ndof)
        return P

    def getResiduum(self):
        """
        **NEEDS REDESIGN TO WORK WITH SMART NODES**
        """

        # R = self.R.copy()
        # R.shape = (self.nNodes, self.ndof)
        # return R

        return self.R.copy()

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
        super(NewtonRaphsonSolver, self).pushState( state)

        if 'lam1' in state:
            self.loadfactor = state['lam1']
