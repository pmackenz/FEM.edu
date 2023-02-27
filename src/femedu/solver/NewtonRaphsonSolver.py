import numpy as np
import scipy as sp

from ..solver.Solver import Solver

class NewtonRaphsonSolver(Solver):
    """
    An iterative solver for load-controlled nonlinear analysis.
    """

    def __init__(self):
        super(NewtonRaphsonSolver, self).__init__()

    def solve(self, max_steps=10, verbose=False):
        """
        :param max_step: maximum number of iterations (int)
        :param verbose: set to :code:`True` for additional information
        """

        for k in range(max_steps):

        # compute force vector and tangent stiffness
            self.assemble()

            normR = self.checkResiduum(verbose)

            if normR<self.TOL:
                break

            # Solve for equilibrium
            self.solveSingleStep()

        if self.record:
            self.recordThisStep()

        print('+')

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

        # update the system displacements
        self.sysU += dU

    def assemble(self):
        """
        inherited from :code:`Solver` class.
        """
        super(NewtonRaphsonSolver, self).assemble()

    def setLoadFactor(self, lam):
        """
        Set the target load factor to **lam**
        """
        self.loadfactor = lam

    def getDisplacements(self):
        U = self.sysU.copy()
        U.shape = (self.nNodes, self.ndof)
        return U

    def getForces(self):
        P = self.P.copy() * self.loadfactor
        P.shape = (self.nNodes, self.ndof)
        return P

    def getResiduum(self):
        R = self.R.copy()
        R.shape = (self.nNodes, self.ndof)
        return R
