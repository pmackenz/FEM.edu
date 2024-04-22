import numpy as np
import os
import sys

from .Element import *
from ..recorder.Recorder import Recorder
from .Face2D import *
from .Face3D import *


class LinearElement(Element):
    """
    abstract class: representing a single generic linear element

    The LinearElement adds functions to emulate a nonlinear element
    within the given framework.  Functions include
    * computing the internal force vector for a linear formulation
    *

    """

    def __init__(self, nodes, material, label=None):
        """

        :param nodes:
        :param material:
        :param label:
        """
        super(LinearElement, self).__init__(nodes, material, label=label)

    def initialize(self, type=DrawElement.UNKNOWN, dofs=[] ):
        """

        :param type: define `DrawElement.TYPE` for plotting. See :ref:`DrawElement_class` for available options.
        :param dofs: list of dofs used by this element. Example: `["ux","uy"]` for 2d-displacements.
        """
        self.element_type = type

        if dofs:
            self._requestDofs(dofs)
        else:
            raise TypeError("mandatory list of dofs is missing for element {}".format(self))

        self.createFaces()

        nnodes = len(self.nodes)
        ndofs  = len(dofs)

        self.distributed_load = [ 0.0 for i in range(nnodes) ]
        self.force    = 0.0
        self.Forces   = [ np.zeros(ndofs) for k in range(len(self.nodes)) ]
        self.Kt       = [ [ np.zeros((ndofs,ndofs)) for k in range(nnodes) ] for m in range(nnodes) ]

    def getForce(self):
        """
        Request the internal force vector (stress driven force only; **no applied element loads**)

        :return:
        """
        self.updateState()
        return self.Forces


