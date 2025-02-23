import numpy as np

from ..Element import *
from ...domain.Node import *

class BeamSolidLink(Element):
    r"""
    class: linking a node from a 2d model (follower) to a frame node (lead)

    **Assumptions**

    * The frame and the 2d-plate are in the X-Y-plane.
    * Small displacements and moderate rotations
    * Navier's assumptions (plane sections)

    The element is using dofs
    * :math:`u` (:code:`ux`), :math:`v` (:code:`uy`) and :math:`\theta` (:code:`rz`) on the frame node
    * :math:`u` (:code:`ux`) and :math:`v` on the 2d-plate node
    * one unique lagrange multiplier per link (no generic variable name)

    .. list-table:: Internal variables

        * - self.nodes
          - nodes i and j (tuple)
        * - self.material
          - pointer to Material object
        * - self.force
          - internal force (list of arrays for axial force, _f_, shear, _V_, and moment, _M_)
        * - self.Forces
          - nodal force vectors (list of np.arrays)
        * - self.Kt
          - tangent stiffness (list of np.arrays)

    """

    COUNT = 0

    def __init__(self, frame_node, plate_node, label=""):
        """
        :param frame_node: pointer to the lead node on the frame model
        :param plate_node: pointer to the follower node on the 2d-plate
        """
        super().__init__((frame_node, plate_node), None, label=label)
        self.element_type = DrawElement.LINE

        dim = plate_node.getPos().size
        if dim != 2:
            raise TypeError("spatial dimension of nodes must be 2")

        # genrate a unique Lagrange Multiplier
        lm1_name = self.getUniqueLMName()
        lm2_name = self.getUniqueLMName()

        dof_list_lead = ('ux', 'uy', 'rz')
        dof_idx = frame_node.request(dof_list_lead, self)
        self.dof_idx[frame_node] = dof_idx  # marked for removal

        dof_list_follow = ('ux', 'uy', lm1_name, lm2_name)
        dof_idx = plate_node.request(dof_list_follow, self)
        self.dof_idx[plate_node] = dof_idx  # marked for removal

        self._dof_list = (dof_list_lead, dof_list_follow)

        # initialize element load to zero
        self.distributed_load = 0.0

        rf            = self.nodes[1].getPos() - self.nodes[0].getPos()  # xf - xl
        self.L0       = np.linalg.norm(rf)
        self.force    = 0.0
        self.Forces   = [np.zeros(3), np.zeros(4)]

        Kll = np.zeros((3, 3))

        Klf = np.zeros((3, 4))
        Klf[0,2] = -1.
        Klf[1,3] = -1.
        # (rx)   (hx)   ( ry hz         )   [  0   0  ry ](hx)
        # (ry) x (hy) = (-rx hz         ) = [  0   0 -rx ](hy)
        # ( 0)   (hz)   ( rx hy - ry hx )   [-ry  rx   0 ](hz)
        Klf[2,2] =  rf[1]
        Klf[2,3] = -rf[0]

        Kff = np.zeros((4, 4))
        Kff[0,2] = 1.
        Kff[1,3] = 1.
        Kff[2,0] = 1.
        Kff[3,1] = 1.

        self.Kt = [[ Kll  , Klf ],
                   [ Klf.T, Kff ]]

        Glead = [ [-1., 0.], [ 0.,-1.], [rf[1],-rf[0]] ]
        Gfollower = [ [1.,0.], [0.,1.], [0.,0.], [0.,0.] ]
        self.G  = [np.array(Glead), np.array(Gfollower) ]

        self.internal_forces = {'fx plate':0.0, 'fy plate':0.0, 'fx frame':0.0, 'fy frame':0.0, 'Mz frame':0.0}

    def __str__(self):
        s = super(BeamSolidLink, self).__str__()
        s += "\n    internal forces: fx plate={fx plate:.2f} fy plate={fy plate:.2f} fx frame={fx frame:.2f} fy frame={fy frame:.2f} Mz frame={Mz frame:.2f}".format(**self.internal_forces)
        return s

    def getDofs(self, node=None):
        r"""
        returns the dof-codes for this element in a list
        """
        if not node:
            msg = "BeamSoliLink.getDofs requires information about the calling Node. None given"
            raise TypeError(msg)

        if node == self.nodes[0]:
            return self._dof_list[0]
        elif node == self.nodes[1]:
            return self._dof_list[1]
        else:
            msg = "Calling Node object ({}) not part of this BeamSoliLink (linking {} and {})"
            raise TypeError(msg.format(node.getID(), self.nodes[0].getID(), self.nodes[1].getID()))

    def getUniqueLMName(self):
        name = f"LM_BS_{BeamSolidLink.COUNT:d}"
        BeamSolidLink.COUNT += 1
        return name

    def getInternalForce(self, variable=''):
        """
        Dummy implementation since this element is not printable.

        :returns: tuple (None, None)
        """
        # self.updateState()

        return (None, None)

    def updateState(self):

        # 0 ... ux
        # 1 ... uy
        # 2 ... theta
        dispi = self.getDisp(0)

        # 0 ... ux_follower
        # 1 ... uy_follower
        # 2 ... LM_x (Lagrange multiplier)
        # 3 ... LM_y (Lagrange multiplier)
        dispj = self.getDisp(1)
        LM = dispj[2:4]

        # ** build element load vector and element tangent stiffness

        # .. internal force
        self.Fl = self.G[0] @ LM
        self.Ff = self.G[1] @ LM
        self.Forces = [self.Fl, self.Ff]

        # internal forces at nodes
        # .. plate is the follower
        # .. frame is the lead
        self.internal_forces = {
            'fx plate': self.Ff[0],
            'fy plate': self.Ff[1],
            'fx frame': self.Fl[0],
            'fy frame': self.Fl[1],
            'Mz frame': self.Fl[2]
        }