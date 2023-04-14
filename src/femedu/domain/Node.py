import numpy as np
from .Transformation import *
from ..recorder.Recorder import Recorder

class Node():
    """
    class: representing a single Node
    """
    COUNT = 0

    def __init__(self, x0, y0, z0=None):
        """

        :param x0: Initial position (List)

        """
        self.ID = Node.COUNT
        Node.COUNT += 1

        if isinstance(z0, (int, float)):
            self.pos = np.array([x0, y0, z0])
        else:
            self.pos = np.array([x0, y0])

        self.is_lead     = True   # is this a lead node?  Will be set to follower (is_lead = False) if tied
        self.lead        = self   # following yourself
        self.followers   = []     # list of following nodes

        self.disp        = None   # active current displacement vector
        self.disp_n      = None   # previously converged displacement vector
        self.disp_pushed = None   # stored displacement vector (see pushU() and popU())
        self.disp_mode   = None   # stored displacement representing a mode shape
        self.dofs        = {}
        self.ndofs       = 0
        self.start       = None
        self.elements    = []
        self._fixity     = []
        self.loads       = {}
        self._hasLoad    = False
        self.transform   = None    # nodal transformation object
        self.dof_maps    = {}      # dof_idx maps for attached elements

        self.setRecorder(None)

        self.setLoadFactor(1.0)

    def __str__(self):
        s  = "Node_{}:\n    x:    {}".format(self.ID, self.pos)
        if not self.is_lead:
            s +=  "\n    following {}".format(self.lead.getID())
        if self._fixity:
            s += f"\n    fix:  {self._fixity}"
        load = self.getLoad()
        if isinstance(load, np.ndarray) and not np.isclose(np.linalg.norm(load), 0.0):
            s += f"\n    P:    {load}"
        s += f"\n    u:    {self.disp}"
        return s

    def __repr__(self):
        return "Node_{}(x={}, u={})".format(self.ID, self.pos, self.disp)

    def getID(self):
        """
        :returns: the node ID (``str``)
        """
        return "Node_{}".format(self.ID)

    def getLeadID(self):
        """
        :returns: the node ID of the **lead** node (``str``)
        """
        return self.lead.getID()

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
        if self.is_lead:
            dof_idx = []
            for dof in dof_list:
                if dof not in self.dofs:
                    self.dofs[dof] = self.ndofs
                    self.ndofs += 1
                dof_idx.append(self.dofs[dof])

            if caller not in self.elements:
                self.elements.append(caller)

            # remember the dof_idx map for this element for future interaction
            self.dof_maps[caller] = dof_idx

            return tuple(dof_idx)

        else:
            return self.lead.request(dof_list=dof_list, caller=caller)

    def fixDOF(self, *dofs):
        """
        provide a list of dof codes that shall be restrained

        see also: :code:`request()`

        :param dofs: fix the dof defined by key(s)
        """
        if self.is_lead:
            for dof in dofs:
                if isinstance(dof, str):
                    if dof not in self._fixity:
                        self._fixity.append(dof)
                elif isinstance(dof,list) or isinstance(dof,tuple):
                    for item in dof:
                        self.fixDOF(item)
                else:
                    raise TypeError
        else:
            self.lead.fixDOF(*dofs)

    def __floordiv__(self, other):
        """

        :param other:
        :return: self
        """
        if self.is_lead:
            self.fixDOF(other)
            return self
        else:
            return self.lead.__floordiv__(other)

    def isFixed(self, dof):
        """

        :param dof: dof code as defined in :code:`request()`
        """
        if self.is_lead:
            return (dof in self._fixity)
        else:
            return self.lead.isFixed(dof)

    def areFixed(self):
        """
        To be used by the assembly routines.

        return a list of indices pointing to fixed dofs in this node.
        Indices are local to this node: :code:`0..num_dofs`
        """
        if self.is_lead:
            idx = []
            for i,dof in enumerate(self._fixity):
                idx.append(i)
            return np.array(idx, dtype=int)
        else:
            return self.lead.areFixed()


    def getFixedDofs(self):
        """
        :returns: a list of fixed dofs by dof-code strings.
        """
        if self.is_lead:
            return self._fixity[:]
        else:
            return self.lead.getFixedDofs()

    def setStart(self, startInt):
        if self.is_lead:
            self.start = startInt
        else:
            msg = "Illegal attempt to add {this_node} to the global system: {this_node} is a follower.".format(this_node=self.getID())
            raise TypeError(msg)

    def setDisp(self, U, dof_list=None, modeshape=False):
        """
        use for prescribed displacements

        **NEEDS TO BE IMPLEMENTED**

        :param U:
        :param modeshape: set to True if U represents a mode shape (BOOL)
        :param dof_list:
        """
        if self.is_lead:
            if isinstance(U,list) or isinstance(U,tuple):
                U = np.array(U)

            if dof_list:
                for ui, dof in zip(U, dof_list):
                    if dof in self.dofs:
                        if modeshape:
                            self.disp_mode[self.dofs[dof]] = ui
                        else:
                            self.disp[self.dofs[dof]] = ui
                    else:
                        msg = f"requested dof:{dof} not present at current node.  Available dofs are {self.dofs.keys()}"
                        raise TypeError(msg)
            else:
                if modeshape:
                    self.disp_mode = U
                else:
                    self.disp = U

        else:
            self.lead.setDisp(U, dof_list=dof_list, modeshape=modeshape)

    def _updateDisp(self, dU):
        """
        For internal use by solvers only

        :param dU: displacement correction from last iteration step.
        """
        if self.is_lead:
            self.disp += dU

        """
        Do not forward that call to the lead node or that increment will be duplicated.
        """

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

    def getDisp(self, dofs=None, caller=None, **kwargs):
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
        if self.is_lead:
            if 'modeshape' in kwargs and kwargs['modeshape']:
                if not isinstance(self.disp_mode, np.ndarray):
                    self.disp_mode = np.zeros(self.ndofs)
                U = self.disp_mode
            else:
                if not isinstance(self.disp, np.ndarray):
                    self.disp = np.zeros(self.ndofs)
                U = self.disp

            if caller:
                # we know the calling element.
                # ... ignoring dofs and using dof list from element map

                if caller not in self.dof_maps:
                    msg = "caller not registered with this node"
                    raise TypeError(msg)

                idx = self.dof_maps[caller]
                ans = U[idx]
                return np.array(ans)

            else:
                # we do not know who is requesting displacements, so provide all requested or ALL if no dofs were specified.
                if dofs:
                    if isinstance(dofs,str):
                        dofs = [dofs]
                    ans = []
                    for dof in dofs:
                        if dof in self.dofs:
                            ans.append(U[self.dofs[dof]])
                        else:
                            ans.append(0.0)
                else:
                    ans = U

                return np.array(ans)

        else:
            self.lead.getDisp(dofs=dofs, caller=caller, **kwargs)

    def getPos(self, caller=None, **kwargs):
        """
        :param caller: pointer to calling element
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
        if self.is_lead:
            if self.pos.shape[0] == 1:
                my_dofs = ('ux',)
            elif self.pos.shape[0] == 2:
                my_dofs = ('ux','uy')
            elif self.pos.shape[0] == 3:
                my_dofs = ('ux','uy','uz')
            else:  # don't know what that would be
                msg = f"Don't know how to interprete position: {self.pos}"
                raise TypeError(msg)

            if 'modeshape' in kwargs and kwargs['modeshape']:
                if not isinstance(self.disp_mode, np.ndarray):
                    self.disp_mode = np.zeros(self.ndofs)

                return self.pos + factor * self.getDisp(caller=None, dofs=my_dofs, modeshape=1)
            else:
                if not isinstance(self.disp, np.ndarray):
                    self.disp = np.zeros(self.ndofs)

                return self.pos + factor * self.getDisp(caller=None, dofs=my_dofs)

        else:
            return self.lead.getDeformedPos(caller=caller, factor=factor, **kwargs)

    def getIdx4Element(self, elem):
        if self.is_lead:
            if elem in self.dof_maps:
                return np.array(self.dof_maps[elem], dtype=int)
            else:
                msg = f"Element {elem} not in dof_map for node {self.ID}"
                raise TypeError(msg)
        else:
            return self.lead.getIdx4Element(elem)

    def getIdx4DOFs(self, dofs=[]):
        if self.is_lead:
            idx = []
            for dof in dofs:
                if dof in self.dofs:
                    idx.append(self.dofs[dof])
                else:
                    msg = f"dof {dof} not present at node {self.ID}"
                    raise TypeError(msg)

            return np.array(idx, dtype=int)

        else:
            return self.lead.getIdx4DOFs(dofs=dofs)

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
        Nodes are applied as components energetically linked to any degree of freedom provided at this node.

        For example, in a 2D-problem with the x-axis to the right and the y-axis pointing up,
        applying a vertical downward force, P, and a counter-clockwise moment, M, at a node, xn:

        .. code::

            xn.addLoad([-P, M],['uy','rz'])

        :param loads:
        :type loads: list of floats
        :param dofs:
        :type dofs: list of dof-codes
        """
        if self.is_lead:
            # Check tuple type and if the dof exists (warn and continue)
            for (load, dof) in zip(loads, dofs):
                if dof in self.loads:
                    self.loads[dof] += load
                else:
                    self.loads[dof] = load
            self._hasLoad = True
        else:
            self.lead.addLoad(loads, dofs)

    def setLoad(self, loads, dofs):
        """

        :param loads:
        :param dofs:
        """
        if self.is_lead:
            # Check tuple type and if the dof exists (warn and continue)
            for (load, dof) in zip(loads, dofs):
                self.loads[dof] = load
            self._hasLoad = True
        else:
            self.lead.setLoad(loads, dofs)

    def resetLoad(self):
        """

        """
        if self.is_lead:
            self.loads = {}
            self._hasLoad = False
        else:
            self.lead.resetLoad()

    def getLoad(self, dof_list=None, apply_load_factor=False):
        """
        :returns: nodal load vector (ndarray)
        """
        if self.is_lead:
            force = np.zeros(self.ndofs)
            for dof in self.loads:
                if dof in self.dofs:
                    if apply_load_factor:
                        force[self.dofs[dof]] = self.loads[dof] * self.loadfactor
                    else:
                        force[self.dofs[dof]] = self.loads[dof]
        else:
            force = self.lead.getLoad(dof_list=dof_list, apply_load_factor=apply_load_factor)

        return force

    def hasLoad(self):
        """
        :returns: **True** if this node has **any** loads (bool)
        """
        if self.is_lead:
            return self._hasLoad
        else:
            return self.lead.hasLoad()

    def resetDisp(self):
        """

        """
        if self.is_lead:
            self.disp = np.zeros(len(self.dofs))
        else:
            self.lead.resetDisp()

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

        if not self.is_lead:
            self.lead.setLoadFactor(lam)

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

    def on_converged(self):
        pass

    def revert(self):
        pass

    def isLead(self):
        return self.is_lead

    def make_follower(self, lead):
        """
        this command tells the current node (**self**) to "follow" whatever the **lead** node is doing.

        :param lead: pointer to the lead node
        :type lead: Node
        """

        self.lead = lead

        if self == lead:
            self.is_lead = True
        else:
            self.is_lead = False
            lead.addFollower(self)

        # transfer element maps
        for elem in self.dof_maps:
            elem_dof_map = []
            for dof in self.dofs:
                idx = self.dofs[dof]
                if idx in self.dof_maps[elem]:
                    elem_dof_map.append(dof)
            lead.request(elem_dof_map, elem)

        # transfer nodal loads
        if self._hasLoad:
            self.lead.addLoad(self.loads.values(), self.loads.keys())
            self.loads = {}
            self._hasLoad = False

        # transfer fixities
        self.lead.fixDOF(self._fixity)

    def addFollower(self, follower):
        if follower in self.followers:
            msg = "{} is already following {}. Circular tie?".format(follower.getID(), self.getID())
            raise TypeError(msg)
        else:
            self.followers.append(follower)



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
