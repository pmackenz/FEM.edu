from scipy.sparse import coo_array
from scipy.sparse.linalg import spsolve
from scipy.sparse.csgraph import reverse_cuthill_mckee

from .NewtonRaphsonSolver import *

class NewtonRaphsonSolverSparse(NewtonRaphsonSolver):
    r"""
    An iterative solver for general nonlinear analysis.

    This implementation only uses sparse data models and methods provided by :code:`scipy`.
    The implementation is somewhat harder to read than its close relative :code:`NewtonRaphsonSolver`.
    The benefit is a significantly lower demand for memory and a faster equation solver.
     """

    def __init__(self):
        super(NewtonRaphsonSolverSparse, self).__init__()

    def solveSingleStep(self):
        r"""
        Helper function performing a single solution of the linearized system

        Called by **solve()**. (internal use only)
        """

        # are we doing displacement control?
        if self.hasConstraint:

            # solve for displacement update: a single Newton step
            # dQ = np.linalg.solve(self.Kt, np.stack([self.R, self.P]).T)
            dQ = spsolve(self.Kt, np.stack([self.R, self.P]).T)   # from scipy

            if self.useArcLength:
                # arc-length control
                dload = self.loadfactor - self.loadfactor_n
                denum = 2.*self.alpha * dload * self.P @ self.P
                numerator = self.g

                for node in self.nodes:
                    idx = node.getIdx4DOFs()
                    delU = node.getDeltaU()
                    denum     += 2.*np.dot(dQ[idx,1], delU)
                    numerator -= 2.*np.dot(dQ[idx,0], delU)

                dlam = numerator / denum

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

            # print(self.R)
            dU = spsolve(self.Kt, self.R)   # from scipy
            # print(dU)

        # update nodal displacements
        for node in self.nodes:
            idxK = node.getIdx4DOFs()
            node._updateDisp(dU[idxK])


    def assemble(self, force_only=False):
        r"""
        A general assembler for mixed element types.

        This method will build the out-of-balance force vector (residuum :math:`{\bf R}`)
        and the tangent stiffness matrix (:math:`{\bf K}_t`) used by most solvers.

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

        # initialize arrays
        Psys = np.zeros(ndof)  # reference load vector (without load factor)
        Fsys = np.zeros(ndof)  # system internal force vector

        rows = []        # row indices for coo_array
        cols = []        # column indices for coo_array
        Ksys_data = []   # system tangent stiffness matrix

        # displacement control or arc-length control
        #
        # ### needs work ### #
        #
        if self.hasConstraint:
            if self.useArcLength:
                # arc-length control
                dload = self.loadfactor - self.loadfactor_n
                self.g = self.arclength2 - self.alpha * dload*dload * self.P@self.P
                for node in self.nodes:
                    self.g -= node.getNormDeltaU2()
            else:
                self.g = self.control_node.getDisp(self.control_dof) - self.targetU

        # assemble loads
        for node in self.nodes:
            if node.isLead() and node.hasLoad():
                idx = node.getIdx4DOFs()
                Psys[idx] += node.getLoad()

        # Element Loop: assemble element forces and stiffness
        for element in self.elements:

            Fe = element.getForce()     # Element State Update occurs here
            Pe = element.getLoad()      # Element State Update occurs here

            if not force_only:
                Ke = element.getStiffness() # fetch element stiffness matrix as array of nodal matrices

            for (i,ndI) in enumerate(element.nodes):
                # fetch dof mapping for node I
                idxK = ndI.getIdx4Element(element)

                # system reference load vector
                if isinstance(Pe[i], np.ndarray):
                    # assemble the load vector
                    Psys[idxK] += Pe[i]

                # system residual force vector
                Fsys[idxK] += Fe[i]

                # system tangent stiffness matrix
                if not force_only:
                    for (j,ndJ) in enumerate(element.nodes):
                        # fetch dof mapping for node J
                        idxM = ndJ.getIdx4Element(element)
                        # add to system matrix
                        rows += idxK.repeat(len(idxM)).flatten().tolist()
                        cols += idxM.tolist()*len(idxK)
                        Ksys_data += Ke[i][j].flatten().tolist()

        # system residual force vector
        self.P = Psys
        self.R = self.loadfactor * Psys - Fsys

        # apply boundary conditions
        if not force_only:

            maxKij = np.array(Ksys_data).max()

            for node in self.nodes:
                for dof in node.dofs:
                    if node.isFixed(dof):
                        #idx = node.lead.start + node.dofs[dof]
                        idx = node.getIdx4DOFs(dofs=[dof])[0]

                        ubar  = node.getFixedDisp(dof, local=True)
                        ubar -= node.getDisp(dof, local=True)

                        self.R[idx]   += ubar * 1.0e20 * maxKij
                        rows.append(idx)
                        cols.append(idx)
                        Ksys_data.append(1.0e20 * maxKij)

            rows = np.array(rows, dtype=int)
            cols = np.array(cols, dtype=int)
            Ksys_data = np.array(Ksys_data)
            self.Kt = coo_array((Ksys_data, (rows,cols))).tocsr()
