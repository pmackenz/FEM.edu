"""
==============================================
A 1d linear spring element
==============================================

The spring is defined as the horizontal connection between
two nodes :py:obj:`ndi` and :py:obj:`ndj`.

The mechanics of the element is described by length change

.. math::

   \delta = U_{j} - U_{i}

and internal force

.. math::

   f = c\:\delta

where :math:`c` is the spring stiffness in force per length change.


"""
import numpy as np
from ..Element import Element
from ...materials import Material

class Spring(Element):
    """
    A 1d linear elastic spring element.

    This element only uses the x-coordinate and the displacement in x-direction ('ux').
    """

    def __init__(self, ndi, ndj, c=1, label=None):
        """
        :param ndi:  node object
        :param ndj:  node object
        :patam c:    spring stiffness (force/length change)
        """
        super(Spring, self).__init__((ndi, ndj), Material(params={'E':1.0, 'A':c}), label=label)
        self.element_type = Element.LINE
        self.c = c

        self._requestDofs(['ux'])

        # initialization
        self.force  = 0.0           # internal force
        self.delta  = 0.0           # length change
        self.Forces = [0.0,0.0]     # nodal forces
        self.Kt = [[c, -c],[-c, c]] # spring stiffness matrix

    def __str__(self):
        """
        nice format for printing this spring
        """
        s = "Spring {}: {} to {} with c={}".format(self.getID(), self.nodes[0].getID(), self.nodes[1].getID(), self.c)
        s += f"\n    length change:  delta = {self.delta} "
        s += f"\n    internal force: force = {self.force}"
        return s

    def __repr__(self):
        """
        nice format for debugging
        """
        s = "Spring({}, {}, c={})".format(self.getID(), self.nodes[0].getID(), self.nodes[1].getID(), self.c)
        return s

    def updateState(self):
        """
        called to compute internal force and current nodal reactions for equilibrium test

        :return: no return values needed (inherited from Element class)
        """
        Ui = self.getDisp(0)[0]  # this functions ALWAYS returns a list. We need only the first element
        Uj = self.getDisp(1)[0]  # this functions ALWAYS returns a list. We need only the first element

        # length change (elongation positive)
        self.delta = Uj - Ui

        # internal force (tension positive)
        self.force = self.c * self.delta

        # nodal forces (global coordinates)
        self.Forces = [-self.force, self.force]

    def getInternalForce(self, variable=''):
        """
        Helper function used by the plotting functions.

        :return: :py:obj:`(s,f)` , two lists containing start and end values of normalized position and internal force,
                respectively, at both the start and the end point of the line of action.
        """
        s = np.array( [0.0, 1.0] )
        f = np.array( [-self.Forces[0], self.Forces[1]] )
        return (s, f)

