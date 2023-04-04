import numpy as np
import scipy.sparse as scs
import scipy.sparse.linalg as spla

from ..solver.NewtonRaphsonSolver import NewtonRaphsonSolver

class SparseSolver(NewtonRaphsonSolver):
    """
    A sparse solver version for load-controlled nonlinear analysis.
    """

    def __init__(self):
        super(SparseSolver, self).__init__()

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
            ##dQ = np.linalg.solve(self.Kt, np.stack([self.R, self.P]).T)# LU factorization of Kt
            solve = spla.factorized(self.Kt)

            # solve for displacement update: a single Newton step
            dq0 = solve(self.R)
            dq1 = solve(self.P)

            if self.useArcLength:
                # arc-length control
                delU  = self.sysU - self.lastConverged['U']
                dload = self.loadfactor - self.lastConverged['lambda']
                # g = self.arclength2 - delU@delU - self.alpha * dload*dload * self.P@self.P

                denum = 2.*np.dot(dq1, delU) + 2.*self.alpha * dload * self.P@self.P

                dlam = ( self.g - 2.*np.dot(dq0, delU) ) / denum

            else:
                # displacement control
                # g = (self.targetU - np.dot(self.sysU, self.en))

                dlam = self.g - np.dot(dq0, self.en)
                dlam /= np.dot(dq1, self.en)

            self.loadfactor += dlam
            dU = dq0 + dlam * dq1

        else:
            dU = spla.spsolve(self.Kt, self.R)

        # update nodal displacements
        for node in self.nodes:
            idxK = node.start + np.arange(node.ndofs)
            node._updateDisp(dU[idxK])

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

        for constraint in self.constraints:
            constraint.setStart(ndof)
            ndof += constraint.countConditions()

        self.sdof = ndof  # number of system d.o.f.s

        Rsys = np.zeros(ndof)

        rows = []
        cols = []
        data = []

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
                        KIJ = element.Kt[i][j]
                        II,JJ = KIJ.shape
                        for p in range(II):
                            for q in range(JJ):
                                rows.append(idxK[p])
                                cols.append(idxM[q])
                                data.append(KIJ[p][q])

        # apply boundary conditions
        if not force_only:
            for node in self.nodes:
                for dof in node.dofs:
                    if node.isFixed(dof):
                        idx = node.start + node.dofs[dof]
                        Rsys[idx]      = 0.0
                        rows.append(idx)
                        cols.append(idx)
                        data.append(1.0e30)

        self.R  = Rsys

        KtS = scs.coo_array((data, (rows, cols)), shape=(self.sdof,self.sdof))
        self.Kt = scs.csc_array(KtS)
