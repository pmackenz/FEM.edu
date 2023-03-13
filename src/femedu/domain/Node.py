import numpy as np
from .Transformation import *
from ..recorder.Recorder import Recorder

class Node():
    """
    class: representing a single Node
    """

    def __init__(self, x0, y0, z0=None):
        """

        :param x0: Initial position (List)

        """
        if isinstance(z0, (int, float)):
            self.pos = np.array([x0, y0, z0])
        else:
            self.pos = np.array([x0, y0])

        self.index       = -1
        self.disp        = None   # active current displacement vector
        self.disp_n      = None   # previously converged displacement vector
        self.disp_pushed = None   # stored displacement vector (see pushU() and popU())
        self.disp_mode   = None   # stored displacement representing a mode shape
        self.dofs        = {}
        self.ndofs       = 0
        self.start       = None
        self.elements    = []
        self.fixity      = []
        self.loads       = {}
        self._hasLoad    = False
        self.transform   = None    # nodal transformation object

        self.setRecorder(None)

        self.setLoadFactor(1.0)

    def __str__(self):
        s  = "Node_{}:\n    x:    {}".format(self.index, self.pos)
        if self.fixity:
            s += f"\n    fix:  {self.fixity}"
        load = self.getLoad()
        if isinstance(load, np.ndarray) and not np.isclose(np.linalg.norm(load), 0.0):
            s += f"\n    P:    {load}"
        s += f"\n    u:    {self.disp}"
        return s

    def __repr__(self):
        return "Node_{}(x={}, u={})".format(self.index, self.pos, self.disp)

    def getID(self):
        """
        :returns: the node ID (``str``)
        """
        return "Node_{}".format(self.index)

    def request(self, dof_list, caller):
        """
        send list or individual dof code. Common codes:

        .. list-table::
            :header-rows: 1

            * - code
              - description
            * - **ux**
              - displacement in x-direction
            * - **uy**
              - displacement in y-direction
            * - **uz**
              - displacement in z-direction
            * - **rx**
              - rotation about x-axis
            * - **ry**
              - rotation about y-axis
            * - **rz**
              - rotation about z-axis


        :param dof_list:  list of dof-codes required by calling element
        :param caller:  pointer to calling element (usually sent as self)
        """
        dof_idx = []
        for dof in dof_list:
            if dof not in self.dofs:
                self.dofs[dof] = self.ndofs
                self.ndofs += 1
            dof_idx.append(self.dofs[dof])

        if caller not in self.elements:
            self.elements.append(caller)

        return tuple(dof_idx)

    def fixDOF(self, *dofs):
        """
        provide a list of dof codes that shall be restrained

        see also: :code:`request()`

        :param dofs:
        """
        for dof in dofs:
            if isinstance(dof, str):
                if dof not in self.fixity:
                    self.fixity.append(dof)
            elif isinstance(dof,list) or isinstance(dof,tuple):
                for item in dof:
                    self.fixDOF(item)
            else:
                raise TypeError


    def __floordiv__(self, other):
        """

        :param other:
        :return: self
        """
        self.fixDOF(other)
        return self

    def isFixed(self, dof):
        """

        :param dof: dof code as defined in :code:`request()`
        """
        return (dof in self.fixity)

    def areFixed(self):
        """
        To be used by the assembly routines.

        return a list of indices pointing to fixed dofs in this node.
        Indices are local to this node: :code:`0..num_dofs`
        """
        idx = []
        for i,dof in enumerate(self.fixity):
            idx.append(i)
        return np.array(idx, dtype=int)

    def setStart(self, startInt):
        self.start = startInt

    def setDisp(self, U, dof_list=None, modeshape=False):
        """
        use for prescribed displacements

        **NEEDS TO BE IMPLEMENTED**

        :param U:
        :param modeshape: set to True if U represents a mode shape (BOOL)
        :param dof_list:
        """
        if isinstance(U,list) or isinstance(U,tuple):
            U = np.array(U)

        if modeshape:
            self.disp_mode = U
        else:
            self.disp = U

    def _updateDisp(self, dU):
        """
        For internal use by solvers only

        :param dU: displacement correction from last iteration step.
        """
        self.disp += dU

    def pushU(self):
        """
        Store the current displacement vector for later restore using :code:`popU()`.
        """
        self.disp_pushed = self.disp

    def popU(self):
        """
        Restore a previously pushed displacement vector (using :code:`pushU()`).
        """
        if isinstance(self.disp_pushed, np.ndarray):
            self.disp = self.disp_pushed
        else:
            raise TypeError("no pushed displacement data available")

    def getDisp(self, caller=None, dofs=None, **kwargs):
        """
        return a vector (nd.array) of (generalized) displacements.

        If a :code:`caller` is given, the element-specific d.o.f.s will be returned as a sequence
        in the same order as given by the element's :code:`request()` call.

        If a :code:`dof_list` is given, d.o.f.s will be returned as the sequence given by that list.
        A zero value will be returned for all d.o.f.s that do not exist at this node.

        .. note::

            A single d.o.f., e.g., "ux", can be requested using :code:`getDisp(dofs=('ux',))`.
            Do not forget the :code:`,` to indicate the definition of a single-element tuple.



        :param caller: pointer to element
        :param dofs: tuple or list of d.o.f. keys
        :return: nodal displacement vector
        """
        if 'modeshape' in kwargs and kwargs['modeshape']:
            if not isinstance(self.disp_mode, np.ndarray):
                self.disp_mode = np.zeros(self.ndofs)
            U = self.disp_mode
        else:
            if not isinstance(self.disp, np.ndarray):
                self.disp = np.zeros(self.ndofs)
            U = self.disp

        if dofs:
            ans = []
            for dof in dofs:
                if dof in self.dofs:
                    ans.append(U[self.dofs[dof]])
                else:
                    ans.append(0.0)
            return np.array(ans)
        else:
            return U

    def getPos(self):
        """
        :return: initial position vector
        """
        return self.pos

    def getDeformedPos(self, caller=None, factor=1.0, **kwargs):
        """
        Return deformed position :math:`{\\bf x} = {\\bf X} + f \\: {\\bf u}`

        If a caller is specified, the node will adjust the output to the d.o.f.s specified by the
        element's request() call. (called during initialization.)

        :param caller: pointer to the calling element.
        :param factor: deformation magnification factor, :math:`f`.
        :return: deformed position vector, :math:`{\\bf x}`.
        """
        if 'modeshape' in kwargs and kwargs['modeshape']:
            if not isinstance(self.disp_mode, np.ndarray):
                self.disp_mode = np.zeros(self.ndofs)

            return self.pos + factor * self.disp_mode
        else:
            if not isinstance(self.disp, np.ndarray):
                self.disp = np.zeros(self.ndofs)

            return self.pos + factor * self.disp

    def addTransformation(self, T):
        """
        Attach a transformation object to this node.
        The transformation defines a local coordinate system.
        If a transformation is given, all loads and prescribed displacements are assumed in that local coordinate system.
        Furthermore, all nodal displacements, velocity, or acceleration will be reported in that local coordinate system.
        """
        if T and isinstance(T, Transformation):
            self.transform = T

    def addLoad(self, loads, dofs):
        """

        :param loads:
        :param dofs:
        """
        # Check tuple type and if the dof exists (warn and continue)
        for (load, dof) in zip(loads, dofs):
            if dof in self.loads:
                self.loads[dof] += load
            else:
                self.loads[dof] = load
        self._hasLoad = True

    def setLoad(self, loads, dofs):
        """

        :param loads:
        :param dofs:
        """
        # Check tuple type and if the dof exists (warn and continue)
        for (load, dof) in zip(loads, dofs):
            self.loads[dof] = load
        self._hasLoad = True

    def resetLoad(self):
        """

        """
        self.loads = {}
        self._hasLoad = False

    def getLoad(self, dof_list=None, apply_load_factor=False):
        """
        :returns: nodal load vector (ndarray)
        """
        force = np.zeros(self.ndofs)
        for dof in self.loads:
            if dof in self.dofs:
                if apply_load_factor:
                    force[self.dofs[dof]] = self.loads[dof] * self.loadfactor
                else:
                    force[self.dofs[dof]] = self.loads[dof]
        return force

    def hasLoad(self):
        """
        :returns: **True** if this node has **any** loads (bool)
        """
        return self._hasLoad

    def resetDisp(self):
        """

        """
        self.disp = np.zeros(len(self.dofs))

    def resetAll(self):
        """

        """
        self.resetDisp()
        self.resetLoad()

    def setLoadFactor(self, lam):
        """
        Set the target load factor to **lam**

        .. warning::

            This method should not be called by the user.
            **USE** :code:`System.setLoadFactor(lam)` instead!

        The entered load pattern is considered a reference load,
        to be multiplied by a scalar load factor, :math:`\lambda`.

        If no load factor is set explicitly, a factor of 1.0 is assumed, i.e., the
        entire entered load is applied in full.
        """
        self.loadfactor = lam

    def setRecorder(self, recorder):
        if isinstance(recorder, Recorder):
            self.recorder = recorder
        else:
            self.recorder = None

    def recordThisStep(self, load_level):
        """
        record current state of the system
        """
        if self.recorder and self.recorder.isActive():
            data = {'lam':self.load_level}
            self.recorder.addData(data)


if __name__ == "__main__":
    # testing the Node class
    node = Node(2.0, 3.5)
    node.index = 42
    node.setLoad(1.2, 3.4)
    node.addLoad(5.6, 7.8)
    node.setDisp(0.1234, -4.321)
    node.fixDOF('uy')   # fixes y-direction

    print(repr(node))
    print(node)
