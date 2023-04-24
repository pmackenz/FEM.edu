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
                dload = self.loadfactor - self.lastConverged['lambda']
                # g = self.arclength2 - delU@delU - self.alpha * dload*dload * self.P@self.P

                denum = 2.*self.alpha * dload * self.P @ self.P
                numerator = 0.0

                for node in self.nodes:
                    idx = node.getIdx4DOFs()
                    delU = node.getDeltaU()
                    denum     += 2.*np.dot(dQ[idx,1], delU)
                    numerator += 2.*np.dot(dQ[idx,0], delU)

                dlam = ( self.g - numerator ) / denum

            else:
                # displacement control
                # g = (self.targetU - np.dot(self.sysU, self.en))
                idx = self.control_node.getIdx4DOFs(dofs=[self.control_dof])[0]

                dlam = self.g - dQ[idx,0]
                dlam /= dQ[idx,1]

            self.loadfactor += dlam
            dU = dQ[:,0] + dlam * dQ[:,1]

        else:
            # solve for displacement update: a single Newton step
            dU = np.linalg.solve(self.Kt, self.R)

        # update nodal displacements
        for node in self.nodes:
            idxK = node.getIdx4DOFs()
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

    # arc-length control

    def initArcLength(self, load_increment=1., alpha=0.0, tolerance=1.0e-12):
        """
        Initializes parameters for the arc-length constraint.

        .. math::

           g({\\bf u}, \lambda)) := \alpha ||\\bar {\bf P}|| (\lambda-\lambda_n)^2 + ({\\bf u} - {\\bf u}_n)({\\bf u} - {\\bf u}_n) - \Delta s^2 = 0


        :param load_increment:   load increment used to calibrate the constraint
        :param alpha:            load contribution factor
        :param tolerance:        convergence tolerance
        """

        # store analysis parameter(s)
        self.alpha = alpha
        self.TOL   = tolerance

        # make sure we start at an equilibrium point
        self.solve()

        # store current configuration as last converged
        for node in self.nodes:
            node.pushU()

        self.lastConverged = {'lambda':0.0, 'ds':0}

        # use load control to solve for the new equilibrium state
        self.hasConstraint = False          # this forces load control
        self.loadfactor += load_increment   # add reference load level
        self.solve()                        # find equilibrium configuration for given load level

        # compute the arc-length for that step and store as target arc length
        g = self.alpha * load_increment**2 * self.P@self.P
        for node in self.nodes:
            g += node.getNormDeltaU2()
        self.arclength2 = g                 # target arc-length

        # # store the current point as last converged point
        # self.previousConverged = self.lastConverged
        # self.lastConverged = {'lambda':0.0, 'ds':np.sqrt(self.arclength2)}

        # set solver parameters
        self.hasConstraint = True
        self.useArcLength  = True

    def stepArcLength(self, verbose=False):

        if not self.hasConstraint:
            # this method makes no sense for load control
            return

        # set solver parameters
        self.useArcLength = True

        # set suitable trial state
        for node in self.nodes:
            node.setTrialState()

        # self.sysU *= 2.0
        # self.sysU -= self.previousConverged['U']
        # self.loadfactor *= 2.0
        # self.loadfactor -= self.previousConverged['lambda']
        #
        # # store current state as "last converged solution"
        # self.previousConverged = self.lastConverged
        # self.lastConverged = Un

        if verbose:
            print(self.lastConverged)

        # solve for next point on the equilibrium path
        self.solve(verbose=True)
