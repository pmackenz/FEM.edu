from .Solver import *

class LinearSolver(Solver):
    """
    A linear system solver

    This solver implies :math:`\{{\\bf P}\} = [{\\bf K}] \{{\\bf u}\}`
    """

    def __init__(self):
        super().__init__()

    def solve(self):
        self.resetDisplacements()
        self.assemble()
        errorNorm = self.solveSingleStep()
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

        # update the system displacements
        self.sysU += dU


    def originalSolve(self):
        """
        Solve system of equations and find state of deformation for the given load.
        """

        # compute size parameters
        ndof = 0
        for node in self.nodes:
            node.setStart(ndof)
            ndof += node.ndofs
        Rsys = np.zeros(ndof)
        Ksys = np.zeros((ndof, ndof))

        # assemble loads
        for node in self.nodes:
            if node.hasLoad():
                K = node.start + np.arange(node.ndofs)
                Rsys[K] += node.getLoad()

        # Element Loop: assemble element forces and stiffness
        for element in self.elements:
            Fe = element.getForce()     # Element State Update occurs here
            for (i,ndI) in enumerate(element.nodes):
                K = ndI.start + np.arange(ndI.ndofs)
                Rsys[K] -= Fe[i]
                for (j,ndJ) in enumerate(element.nodes):
                    M = ndJ.start + np.arange(ndJ.ndofs)
                    Ksys[K[:, np.newaxis], M] += element.Kt[i][j]

        # apply boundary conditions
        for node in self.nodes:
            for dof in node.dofs:
                if node.isFixed(dof):
                    idx = node.start + node.dofs[dof]
                    Rsys[idx]      = 0.0
                    Ksys[:, idx]   = np.zeros(ndof)
                    Ksys[idx, :]   = np.zeros(ndof)
                    Ksys[idx, idx] = 1.0

        # stability check for system matrix
        (vals, vecs) = np.linalg.eig(Ksys)
        for (lam, v) in zip(vals, vecs.T):
            if np.abs(lam) < 1.0e-2:
                print(f"lambda = {lam:16.12e}")
                print(v)

        # solve for displacements
        U = np.linalg.solve(Ksys, Rsys)

        # update nodal displacements
        for node in self.nodes:
            K = node.start + np.arange(node.ndofs)
            node._updateDisp(U[K])

        # recompute residual force
        Rsys = np.zeros(ndof)
        for element in self.elements:
            Fe = element.getForce()     # Element State Update occurs here
            for (i,ndI) in enumerate(element.nodes):
                K = ndI.start + np.arange(ndI.ndofs)
                Rsys[K] -= Fe[i]

        #
        self.Rsys = Rsys
        self.disp = U

