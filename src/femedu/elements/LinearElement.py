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

    def getForce(self):
        """
        Request the internal force vector (stress driven force only; **no applied element loads**)

        :return:
        """
        self.updateState()
        return self.Forces


