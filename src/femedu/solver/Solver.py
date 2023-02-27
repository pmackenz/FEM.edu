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

        # numeric iteration tolerance
        self.TOL = 1.0e-6

        self.hasConstraint = False
        self.targetU       = 0.0

        self.useArcLength = False
        self.lastConverged = {'U':self.sysU.copy(), 'lambda':0.0, 'ds':1.}

    def connect(self, nodes, elems):
        self.nodes = nodes
        self.elements = elems

    def fetch(self):
        """
        Fetch the current :code:`state` of the solver.

        .. list-table:: **state** is defined as a dictionary with the following contents:

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

        return state

    def push(self, state):
        """
        Pushes :code:`state` to the solver.
        The solver will use that data to update it's internal state.

        .. list-table:: **state** is defined as a dictionary with the following contents:

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
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def assemble(self):
        # initialize new residuum and tangent stiffness matrix
        Kt = np.zeros((self.sdof,self.sdof))
        R  = self.loadfactor * self.P.copy()  # initialize R to applied load vector

        # copy and reshape system displacement vector to a list of nodal displacement vectors
        U = self.sysU.copy()
        U.shape = (self.nNodes,self.ndof)

        # loop through elements
        for thisElem in self.elements:

            elem = thisElem['element']   # this is a pointer to the current member
            idxI = thisElem['i']   # this is the node index, not the system dof index
            idxJ = thisElem['j']   # this is the node index, not the system dof index

            #sidxI = np.arange(idxI * self.ndof, (idxI + 1) * self.ndof)  # system dofs for node I
            #sidxJ = np.arange(idxJ * self.ndof, (idxJ + 1) * self.ndof)  # system dofs for node J

            sidxI = elem.dofMap + idxI * self.ndof  # system dofs for node I
            sidxJ = elem.dofMap + idxJ * self.ndof  # system dofs for node J

            # update element displacements
            elem.setDisp(U[idxI],U[idxJ])

            # add element force to system forces
            (fi, fj) = elem.getForce()
            R[sidxI] -= fi         # subtract the resisting (internal) force added by member elem
            R[sidxJ] -= fj         # subtract the resisting (internal) force added by member elem

            # add element stiffness to system stiffness
            KTe = elem.getKt()  # this is the nodal stiffness, not the entire element stiffness matrix
            Kt[sidxI[:,np.newaxis],sidxI] += KTe[0][0]
            Kt[sidxI[:,np.newaxis],sidxJ] += KTe[0][1]
            Kt[sidxJ[:,np.newaxis],sidxI] += KTe[1][0]
            Kt[sidxJ[:,np.newaxis],sidxJ] += KTe[1][1]

        # apply boundary conditions
        for idx in self.fixities:
            Kt[:,idx]   = 0.0
            Kt[idx,:]   = 0.0
            Kt[idx,idx] = 1.0e+1
            R[idx]      = 0.0

        self.Kt = Kt
        self.R  = R


    def solve(self):
        """

        :return:
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def initialize(self):
        """

        :return:
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def reset(self):
        """

        :return:
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def checkStability(self):
        if self.sdof < 10:
            detKt = np.linalg.det(self.Kt)
        else:
            detKt = np.min(np.abs(sp.linalg.eigvals(self.Kt)))
        return detKt

    def getBucklingMode(self):
        w, v = sp.linalg.eigh(self.Kt)
        lam = np.abs(w).min()
        idx = np.argwhere(np.abs(w) == lam)
        lam = w[idx[0]]
        vec = v[:,idx[0]]
        return (lam, vec)

    def checkResiduum(self, report=False):

        # compute residual force and tangent stiffness
        self.assemble()

        normR = np.dot(self.R, self.R)

        # Add constriant violation in case we are using displacement control
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
            print('-')
            if self.hasConstraint:
                print(f"Constraint violation: {self.g:12.4e}")

            R = self.getResiduum()
            for i in range(len(R)):
                print("R({}): {:12.4e} {:12.4e}".format(i,*R[i]))

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

