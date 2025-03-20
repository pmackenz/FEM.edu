from copy import deepcopy

import numpy as np
from collections import deque

from ..elements import Element
from ..recorder.Recorder import Recorder

class Node():
    r"""
    class: representing a single Node

    :param x0: Initial position (List)
    """
    COUNT = 0

    def __init__(self, x0, y0=None, z0=None):
        self.ID = Node.COUNT
        Node.COUNT += 1

        if isinstance(z0, (int, float)):
            self.pos = np.array([x0, y0, z0])
        elif isinstance(y0, (int, float)):
            self.pos = np.array([x0, y0])
        else:
            self.pos = np.array([x0])

        self.is_lead     = True   # is this a lead node?  Will be set to follower (is_lead = False) if tied
        self.lead        = self   # following yourself
        self.followers   = []     # list of following nodes

        self.disp          = None   # active current displacement vector
        self.disp_n        = None   # previously converged displacement vector
        self.loadfactor_n  = 0.0    # load factor for previously converged state
        self.disp_nn       = None   # two steps back converged displacement vector
        self.loadfactor_nn = 0.0    # load factor for two steps back converged state
        self.disp_mode     = None   # stored displacement representing a mode shape
        self.disp_pushed   = deque()   # stored displacement vector (see pushU() and popU())

        self.dofs        = {}
        self.ndofs       = 0
        self.start       = None
        self.elements    = []
        self._fixity     = []  # list for dof keys for fixed dofs
        self._setU       = {}  # prescribed displacement parameters u0 and u1: u[dof] = u0 + loadfactor*u1
        self.loads       = {}
        self._hasLoad    = False
        self._transform  = None    # nodal transformation object
        self.dof_maps    = {}      # dof_idx maps for attached elements

        self._resetGaussPointMap()
        self.setRecorder(None)

        self.setLoadFactor(1.0)

    def __str__(self):
        float_formatter = "{:.3f}".format
        np.set_printoptions(formatter={'float_kind': float_formatter})

        s  = "Node_{}:\n    x:    {}".format(self.ID, self.pos)
        if not self.is_lead:
            s +=  "\n    following {}".format(self.lead.getID())
        if self._transform:
            T = self._transform.getT()
            s += f"\n    local: x={T[:,0]}, y={T[:,1]}"
            if T.shape[1]>2:
                s += f", z={T[:,2]}"
        if self._fixity:
            s += f"\n    fix:  {self._fixity}"
        load = self.getLoad()
        if isinstance(load, np.ndarray) and not np.isclose(np.linalg.norm(load), 0.0):
            s += f"\n    P:    {load}"
        s += f"\n    u:    {self.disp}"

        float_formatter = "{:g}".format
        np.set_printoptions(formatter={'float_kind': float_formatter})

        return s

    def __repr__(self):
        return "Node_{}(x={}, u={})".format(self.ID, self.pos, self.disp)

    def getID(self):
        r"""
        :returns: the node ID (``str``)
        """
        return "Node_{}".format(self.ID)

    def getLead(self):
        r"""
        :returns: pointer to the **lead** node (``Node``)
        """
        return self.lead

    def getLeadID(self):
        r"""
        :returns: the node ID of the **lead** node (``str``)
        """
        return self.lead.getID()

    def request(self, dof_list, caller):
        r"""
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

            # let any attached transformation know that the nodal dof-list has changed.
            if self._transform:
                self._transform.refreshMaps()

            return tuple(dof_idx)

        else:
            return self.lead.request(dof_list=dof_list, caller=caller)

    def fixDOF(self, *dofs):
        r"""
        provide a list of dof codes that shall be restrained

        :param dofs: fix the dof defined by key(s)

        See also:
           * :py:meth:`femedu.domain.Node.Node.setDOF`
           * :py:meth:`femedu.domain.Node.Node.request`
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

    def hasDOF(self, dof):
        """
        :param dof: a dof string ID
        :return: True(False) if dof in(not in) the node's dof-list
        """
        return (dof in self.dofs)

    def setDOF(self, dofs=[], values=[]):
        r"""
        alternative to the `fixDOF()` method that
        allows to set a prescribed value other than `0.0`
        for each degree of freedom.

        :param dofs: list of dofs, defined by key(s), for which displacement values shall be prescribed.
        :param values: list of prescribed values for the given list of dofs.
                       Values may be given as a float or as a list of two values.

                       * If only a single value is given, the respective dof will be set to that value.
                       * If a list is given, the dof will be set to :code:`value[0] + loadfactor*value[1]`

        See also:
           * :py:meth:`femedu.domain.Node.Node.fixDOF`
           * :py:meth:`femedu.domain.System.System.setLoadFactor`
           * :py:meth:`femedu.domain.Node.Node.request`
        """
        if (not dofs) or (len(dofs) != len(values)):
            msg = "setDOF requires matching non-empty lists of dofs and associated values"
            raise TypeError(msg)

        if self.is_lead:
            for dof, val in zip(dofs,values):
                self.fixDOF(dof)
                if isinstance(val,float) or isinstance(val,int):
                    self._setU[dof] = (val, 0.0)
                elif isinstance(val,list) or isinstance(val,tuple):
                    if len(val) == 1:
                        self._setU[dof] = (val[0], 0.0)
                    elif len(val) == 2:
                        self._setU[dof] = val
                    else:
                        raise TypeError(f"cannot interpret input values for {dof}")
                else:
                    raise TypeError
        else:
            self.lead.setDOF(dofs, values)

    def __floordiv__(self, other):
        r"""
        This is a short form for :py:meth:`Node.fixDOF(other)`

        :param other:
        :return: self
        """
        if self.is_lead:
            self.fixDOF(other)
            return self
        else:
            return self.lead.__floordiv__(other)

    def isFixed(self, dof):
        r"""

        :param dof: dof code as defined in :code:`request()`
        """
        if self.is_lead:
            return (dof in self._fixity)
        else:
            return self.lead.isFixed(dof)

    def areFixed(self):
        r"""
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
        r"""
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
        r"""
        use for prescribed displacements

        :param U: list or tuple of prescribed values
        :param dof_list: list or tuple of dof-codes associated with values in U
        :param modeshape: set to True if U represents a mode shape (BOOL)
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
        r"""
        For internal use by solvers only

        :param dU: displacement correction from last iteration step.
        """
        if self.is_lead:

            if self._transform:
                self.disp += self.v2g(dU, self)
            else:
                self.disp += dU

        """
        Do not forward that call to the lead node or that increment will be duplicated.
        """

    def pushU(self):
        r"""
        Store the current displacement vector for later restore using :code:`popU()`.
        """
        self.disp_pushed.append({'U':self.disp.copy(), 'lam':self.loadfactor})
        if len(self.disp_pushed)>2:
            self.disp_pushed.popleft()

    def popU(self):
        r"""
        Restore a previously pushed displacement vector (using :code:`pushU()`).
        """
        if len(self.disp_pushed):
            state = self.disp_pushed.pop()
            self.disp       = state['U']
            self.loadfactor = state['lam']
        else:
            raise TypeError("no pushed displacement data available")

    def getDisp(self, dofs=None, caller=None, **kwargs):
        r"""
        return a vector (`nd.array`) of (generalized) displacements.

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
                    self.disp_n = np.zeros(self.ndofs)
                    self.disp_nn = np.zeros(self.ndofs)
                U = self.disp

            # so far, U is in global coordinates.
            # see if local coordinates were requested
            if 'local' in kwargs and kwargs['local'] == True:
                # check if we have a nodal transformation attached
                # if isinstance(self._transform, Transformation):
                if self._transform:
                    #
                    # see Element.getLoad() and Element.getForce() methods
                    #
                    U = self.v2l(U, self)

            # *** the following code segment was moved to Node.getFixedDisp(...)
            # *** this will be used inside Solver and classes derived from it.

            # # apply prescribed displacements
            # for dof in self._setU:
            #     if dof in self.dofs:
            #         # the index of dof in self.U
            #         idx = self.dofs[dof]
            #         # the prescribed displacement value
            #         ubar = self._setU[dof][0] + self.loadfactor * self._setU[dof][1]
            #         # set prescribed value in nodal U-vector
            #         U[idx] = ubar

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
                    if isinstance(dofs, str):
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
            return self.lead.getDisp(dofs=dofs, caller=caller, **kwargs)

    def getFixedDisp(self, dofs=None, caller=None, **kwargs):
        r"""
        return a vector (`nd.array`) of **user-prescribed** (generalized) displacements.

        If a :code:`caller` is given, the element-specific d.o.f.s will be returned as a sequence
        in the same order as given by the element's :code:`request()` call.

        If a :code:`dof_list` is given, d.o.f.s will be returned as the sequence given by that list.
        A zero value will be returned for all d.o.f.s that do not exist at this node.

        .. note::

            A single d.o.f., e.g., "ux", can be requested using :code:`getFixedDisp(dofs=('ux',))`.
            Do not forget the :code:`,` to indicate the definition of a single-element tuple.

        :param caller: pointer to element
        :param dofs: tuple or list of d.o.f. keys
        :return: nodal displacement vector
        """
        if self.is_lead:
            U = np.zeros(self.ndofs)

            # apply prescribed displacements
            for dof in self._setU:
                if dof in self.dofs:
                    # the index of dof in self.U
                    idx = self.dofs[dof]
                    # the prescribed displacement value
                    ubar = self._setU[dof][0] + self.loadfactor * self._setU[dof][1]
                    # set prescribed value in nodal U-vector
                    U[idx] = ubar

            # so far, U is in global coordinates.
            # see if local coordinates were requested
            if 'local' in kwargs and kwargs['local'] == True:
                # check if we have a nodal transformation attached
                # if isinstance(self._transform, Transformation):
                if self._transform:
                    #
                    # see Element.getLoad() and Element.getForce() methods
                    #
                    U = self.v2l(U, self)

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
                    if isinstance(dofs, str):
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
            return self.lead.getFixedDisp(dofs=dofs, caller=caller, **kwargs)


    def getDeltaU(self, previous_step=False):
        r"""
        :return: delta u = (current u) - (last converged u)
        """
        if self.is_lead:
            if previous_step:
                dU = self.disp_n - self.disp_nn
            else:
                dU = self.disp - self.disp_n
            return dU
        else:
            return self.lead.getDeltaU(previous_step=previous_step)

    def getNormDeltaU2(self, previous_step=False):
        r"""
        :return: norm of delta u = (current u) - (last converged u)
        :return: 0.0 if node is a "follower"
        """
        if self.is_lead:
            dU = self.getDeltaU(previous_step)
            return dU @ dU
        else:
            return 0.0

    def getPos(self, caller=None, **kwargs):
        r"""
        :param caller: pointer to calling element
        :return: initial position vector
        """
        return self.pos

    def getDeformedPos(self, caller=None, factor=1.0, **kwargs):
        r"""
        Return deformed position :math:`{\bf x} = {\bf X} + f \: {\bf u}`

        If a caller is specified, the node will adjust the output to the d.o.f.s specified by the
        element's request() call. (called during initialization.)

        :param caller: pointer to the calling element.
        :param factor: deformation magnification factor, :math:`f`.
        :return: deformed position vector, :math:`{\bf x}`.
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
                    self.disp    = np.zeros(self.ndofs)
                    self.disp_n  = np.zeros(self.ndofs)
                    self.disp_nn = np.zeros(self.ndofs)

                return self.pos + factor * self.getDisp(caller=None, dofs=my_dofs)

        else:
            if self.pos.shape[0] == 1:
                my_dofs = ('ux',)
            elif self.pos.shape[0] == 2:
                my_dofs = ('ux','uy')
            elif self.pos.shape[0] == 3:
                my_dofs = ('ux','uy','uz')
            else:  # don't know what that would be
                msg = f"Don't know how to interprete position: {self.pos}"
                raise TypeError(msg)

            #return self.lead.getDeformedPos(caller=caller, factor=factor, **kwargs)
            disp = self.lead.getDisp(dofs=my_dofs, caller=None, factor=factor, **kwargs)

            return self.pos + factor * disp

    def getIdx4Element(self, elem):
        r"""
        :return: an index list to nodal dofs associated with the attached **elem** within the global dof list
        """
        if self.is_lead:
            if elem in self.dof_maps:
                if self._transform:
                    # Transformation may expand the reach of an element to other components, hence, use all available dofs
                    return self.getIdx4DOFs()
                else:
                    # use the subset of dofs used by this element
                    return self.start + np.array(self.dof_maps[elem], dtype=int)
            else:
                msg = f"Element {elem} not in dof_map for node {self.ID}"
                raise TypeError(msg)
        else:
            return self.lead.getIdx4Element(elem)

    def getIdx4DOFs(self, dofs=[], local=False):
        r"""
        :param dofs:  a list of named dofs for which an index is requested.
                      The returned list will be in the same order as listed in **dofs**.
        :param local: by default return index to the global list of dofs.
                      Set :code:`local=True` to return index to this node's list of dofs.

        :return: an index list to this node's dofs
        """
        if self.is_lead:

            if not dofs:
                dofs = self.dofs.keys()

            idx = []
            for dof in dofs:
                if dof in self.dofs:
                    idx.append(self.dofs[dof])
                else:
                    msg = f"dof {dof} not present at node {self.ID}"
                    raise TypeError(msg)

            ans = np.array(idx, dtype=int)
            if not local:
                ans += self.start

            return ans

        else:
            return self.lead.getIdx4DOFs(dofs=dofs, local=local)

    def getLocalTransformationMap(self, dofs=[]):
        r"""
        :return: Transformation matrix
        """
        #if isinstance(self._transform, Transformation):
        if self._transform:
            T = self._transform.getT()
        else:
            T = None

        return T

    def addTransformation(self, T):
        r"""
        Attach a :py:meth:`femedu.domain.Transformation` object to this node.
        The transformation defines a local coordinate system.

        * If a transformation is given, all loads and prescribed displacements are assumed in that local coordinate system.
        * Furthermore, all nodal displacements, velocity, or acceleration will be reported in that local coordinate system.
        """
        from .Transformation import Transformation
        if isinstance(T, Transformation):
            T.registerClient(self) # register this Node with the transformation
            self._transform = T

    def addLoad(self, loads, dofs):
        r"""
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
        r"""
        :param loads: list of force components
        :param dofs:  associated list of DOFs to which respective loads are to be applied
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
        Resets the load vectors
        """
        if self.is_lead:
            self.loads = {}
            self._hasLoad = False
        else:
            self.lead.resetLoad()

    def getLoad(self, dof_list=None, apply_load_factor=False):
        r"""
        :param dof_list: list (or tuple) if dof-keys for which nodal loads are requested. Fill missing dofs by 0.0.
        :param apply_load_factor: defaults to False -> no factors applied.
        :returns: nodal load vector (ndarray)
        """
        if self.is_lead:
            if dof_list:
                force = np.zeros(len(dof_list))
                for dof in dof_list:
                    if (dof in self.loads) and (dof in self.dofs):
                            idx = dof_list.index(dof)
                            if apply_load_factor:
                                force[idx] = self.loads[dof] * self.loadfactor
                            else:
                                force[idx] = self.loads[dof]
            else:
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
        r"""
        :returns: **True** if this node has **any** loads (bool)
        """
        if self.is_lead:
            return self._hasLoad
        else:
            return self.lead.hasLoad()

    def resetDisp(self):
        """
        Resets the displacement vector.
        """
        if self.is_lead:
            self.disp_nn = np.zeros(len(self.dofs))
            self.disp_n  = np.zeros(len(self.dofs))
            self.disp    = np.zeros(len(self.dofs))
        else:
            self.lead.resetDisp()

    def resetAll(self):
        """
        Resets load and displacement vectors.
        """
        self.resetDisp()
        self.resetLoad()

    def setLoadFactor(self, lam):
        r"""
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

    def isClose(self, x, TOL=1.0e-5):
        r"""
        :param x: :py:obj:`nd.array`
        :return: **True** if x is within **TOL** from this node.
        """
        pos = self.pos
        X = np.array(x)

        this_dim = pos.shape[0]
        x_dim = X.shape[0]

        if x_dim != this_dim:
            msg = f"target dimension ({x_dim}) incompatible with Node dimension ({this_dim})"
            raise TypeError(msg)

        is_close = (np.abs(pos - X) <= TOL).all()

        return is_close

    def distanceTo(self, x):
        r"""
        :param x: :py:obj:`nd.array`
        :return: scalar distance from x
        """
        pos = self.pos
        X = np.array(x)

        this_dim = pos.shape[0]
        x_dim = X.shape[0]

        if x_dim != this_dim:
            msg = f"target dimension ({x_dim}) incompatible with Node dimension ({this_dim})"
            raise TypeError(msg)

        dist = np.linalg.norm(pos - X)

        return dist

    def setRecorder(self, recorder):
        if isinstance(recorder, Recorder):
            self.recorder = recorder
        else:
            self.recorder = None

    def recordThisStep(self, load_level):
        r"""
        record current state of the system
        """
        if self.recorder and self.recorder.isActive():
            # data = {'lam':load_level}
            data = {}
            for var in self.recorder.getVariables():
                value = self.getDisp(var)[0]
                data[var] = value
            self.recorder.addData(data)

    def startRecorder(self):
        if self.recorder:
            self.recorder.enable()

    def pauseRecorder(self):
        if self.recorder:
            self.recorder.disable()

    def stopRecorder(self):
        if self.recorder:
            self.recorder.disable()

    def on_converged(self):
        r"""
        This method is called every time a solver signals a converged solution.
        """
        if self.is_lead:
            # rotate states (n)->(n-1) and current->(n)
            if isinstance(self.disp_n, np.ndarray):
                self.disp_nn = self.disp_n.copy()
            if isinstance(self.disp, np.ndarray):
                self.disp_n  = self.disp.copy()
            self.loadfactor_nn = self.loadfactor_n
            self.loadfactor_n  = self.loadfactor

        # do not duplicate this step for the load node !!

    def revert(self):
        pass

    def isLead(self):
        return self.is_lead

    def make_follower(self, lead):
        r"""
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

    def setTrialState(self):
        r"""
        This will set a suitable trial position for the
        arc-length method to

        .. math::

           {\bf u}^{(0)}_{n+1} = 2{\bf u}^{(\infty)}_{n} - {\bf u}^{(\infty)}_{n-1}

        This should be used by a solution algorithm but not by regular
        user input.
        """
        self.disp       = 2.0 * self.disp_n - self.disp_nn
        self.loadfactor = 2.0 * self.loadfactor_n - self.loadfactor_nn

    #
    # Transformation methods
    #

    def hasTransform(self):
        return (self._transform != None)

    def v2l(self, U, caller=None):
        """
        transform a nodal vector from global to local coordinates

        This is commonly used by Elements to transform nodal forces to the Node's local
        system such that the :code:`Solver.assemble` method receives forces in local coordinates.

        :param U: vector to be transformed
        :param caller: pointer to the calling element
        :return: the transformed vector
        """
        if self.is_lead:

            if self._transform:
                # build a full length nodal vector
                if isinstance(caller, Element):
                    Unode = np.zeros(self.ndofs)                           # initialize to zeros
                    map = self.getIdx4DOFs(caller.getDofs(), local=True)   # get a map from element dogs to node dofs
                    Unode[map] = U                                         # fill in element vector to nodal vector
                else:
                    Unode = deepcopy(U)

                Ulocal = self._transform.v2l(Unode, self)
            else:
                Ulocal = deepcopy(U)

            return Ulocal

        else:
            return self.lead.v2l(U, caller=caller)

    def v2g(self, U, caller=None):
        """
        transform a nodal vector from local to global coordinates

        This is commonly used by the solver to transform nodal displacements from the Node's local to the global
        system such that the :code:`Node` class can store displacements internally in global coordinates.

        :param U: vector to be transformed
        :param caller: pointer to the calling element
        :return: the transformed vector
        """
        if self.is_lead:


            if self._transform:
                # build a full length nodal vector
                if isinstance(caller, Element):
                    Unode = np.zeros(self.ndofs)  # initialize to zeros
                    map = self.getIdx4DOFs(caller.getDofs(), local=True)  # get a map from element dogs to node dofs
                    Unode[map] = U  # fill in element vector to nodal vector
                else:
                    Unode = deepcopy(U)

                Uglobal = self._transform.v2g(Unode, self)
            else:
                Uglobal = deepcopy(U)

            return Uglobal

        else:
            return self.lead.v2g(U, caller=caller)

    def m2l(self, M, caller=None):
        """
        transform a nodal matrix from global to local coordinates

        This is commonly used by Elements to transform nodal forces to the Node's local
        system such that the :code:`Solver.assemble` method receives stiffness in local coordinates.

        :param M: vector to be transformed
        :param caller: pointer to the calling element
        :return: the transformed matrix
        """
        if self.is_lead:

            if self._transform:
                # build a full length nodal vector
                if isinstance(caller, Element):
                    Mnode = np.zeros((self.ndofs,self.ndofs))              # initialize to zeros
                    map = self.getIdx4DOFs(caller.getDofs(), local=True)   # get a map from element dogs to node dofs
                    Mnode[map[:, np.newaxis],map] = M                                     # fill in element matrix to nodal vector
                else:
                    Mnode = deepcopy(M)

                Mlocal = self._transform.m2l(Mnode, self)
            else:
                Mlocal = deepcopy(M)

            return Mlocal

        else:
            return self.lead.m2l(M, caller=caller)

    #
    # Mapping functions for Gauss-point values

    def _resetGaussPointMap(self, var=None):
        r"""
        Initialize or reset variable used for
        mapping gauss-point values to nodes
        """
        self._mapped_variable = var
        self._weighted_value  = 0.0
        self._weight          = 0.0

        if not self.is_lead:
            self.lead._resetGaussPointMap(var=var)

    def _addToMap(self, weight, weighted_value):
        """
        create nodal weighted sum of values and sum of weights
        for mapping of gauss-point values to nodes

        :param weight:
        :param weighted_value:
        """
        self._weight         += weight
        self._weighted_value += weighted_value

        if not self.is_lead:
            self.lead._addToMap(weight, weighted_value)

    def _getMappedValues(self, var=None, ignore_lead=False):
        r"""
        Returns the mapped and weighted nodal value.

        :param var:  optional. If a string identifying the variable is given, this function will return zero i fthe current variable does not match
                     the once provided by the caller.
        :return:     the weighted nodal average for the requested variable.
        """
        if self.is_lead or ignore_lead:
            if var and var != self._mapped_variable:
                ans = 0.0
            else:
                if self._weight > 0.0:
                    ans = self._weighted_value / self._weight
                else:
                    ans = 0.0
            return ans
        else:
            return self.lead._getMappedValues(var)

    def getMappedValue(self, var=None, ignore_lead=False):
        r"""
        Returns the mapped and weighted nodal value.

        :param var:  optional. If a string identifying the variable is given, this function will return zero i fthe current variable does not match
                     the once provided by the caller.
        :return:     the weighted nodal average for the requested variable.
        """
        if self.is_lead or ignore_lead:
            if var and var != self._mapped_variable:
                # reset mapped values
                self._resetGaussPointMap(var)
                # compute map
                for elem in self.elements:
                    # request map to the current node only
                    elem.mapGaussPoints(var, target_node=self)

            if self._weight > 0.0:
                ans = self._weighted_value / self._weight
            else:
                ans = 0.0
            return ans
        else:
            return self.lead.getMappedValue(var)
