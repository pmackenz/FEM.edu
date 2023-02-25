import sys
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm


class AbstractPlotter():

    def __init__(self):
        self.plot3D = False

        self.elements  = []  # element pointer list
        self.nodes     = []  # node pointer list
        self.factor    = 0.0 # displacement factor
        self.variable  = ''  # code for state variable to plot

    def __str__(self):
        return "{}() object".format(self.__class__.__name__)

    def __repr__(self):
        return str(self)

    def setMesh(self, nodes, elems):
        """
        Link the nodes and elements so the ElementPlotter can get information from them.

        :param nodes:    list of node pointers
        :param elements: list of element pointer
        """
        self.nodes    = nodes
        self.elements = elems

    def setDisplacements(self, disp):
        """

        :param disp:
        """
        print("** WARNING ** {}.{} deprecated".format(self.__class__.__name__, sys._getframe().f_code.co_name))
        raise DeprecationWarning

    def setValues(self, vals):
        """

        :param vals:
        """
        print("** WARNING ** {}.{} deprecated".format(self.__class__.__name__, sys._getframe().f_code.co_name))
        raise DeprecationWarning

    def setReactions(self, R):
        """

        :param R: array of nodal force vectors
        """
        print("** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name))
        raise NotImplementedError

    def displacementPlot(self, file=None):
        """
        Create a deformed system plot

        If **file** is given, store the plot to that file.
        Use proper file extensions to indicate the desired format (.png, .pdf)

        :param file: filename (str)
        """
        print("** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name))
        raise NotImplementedError

    def valuePlot(self, variable_name='', deformed=False, file=None):
        """
        Create a plot using colors to identify magnitude of internal force.

        If **file** is given, store the plot to that file.
        Use proper file extensions to indicate the desired format (.png, .pdf)

        :param deformed: True | **False**
        :param file: filename (str)
        """
        print("** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name))
        raise NotImplementedError

    def addForces(self, axs):
        """
        add nodal forces to the plot shown in **axs**

        :param axs: axis on which to plot
        """
        print("** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name))
        raise NotImplementedError

