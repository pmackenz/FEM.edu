
import sys
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm


class AbstractPlotter():
    """
    Defines all supported plot features.
    Implementation is required through derived classes.
    """

    def __init__(self):
        self.plot3D = False

        self.elements  = []  # element pointer list
        self.nodes     = []  # node pointer list
        self.factor    = 0.0 # displacement factor
        self.variable  = ''  # code for state variable to plot
        self.loads     = []  # list of nodal applied loads (from nodes and elements)
        self.reactions = []  # list of nodal reactions (from residual)

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
        msg = "** WARNING ** {}.{} marked deprecated".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise DeprecationWarning(msg)

    def setValues(self, vals):
        """

        :param vals:
        """
        msg = "** WARNING ** {}.{} marked deprecated".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise DeprecationWarning(msg)

    def setReactions(self, R):
        """

        :param R: list or tuple of nodal force vectors
        """
        self.reactions = R

        # msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        # raise NotImplementedError(msg)

    def setNodalLoads(self, P):
        """

        :param P: list or tuple of nodal force vectors
        """
        self.loads = P

        # msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        # raise NotImplementedError(msg)

    def displacementPlot(self, factor=1.0, file=None):
        """
        Create a deformed system plot

        If **file** is given, store the plot to that file.
        Use proper file extensions to indicate the desired format (.png, .pdf)

        :param factor: displacement magnification factor
        :param file: filename (str)
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def valuePlot(self, variable_name='', deformed=False, file=None):
        """
        Create a plot using colors to identify magnitude of internal force.

        If **file** is given, store the plot to that file.
        Use proper file extensions to indicate the desired format (.png, .pdf)

        :param variable_name: string code for variable
        :param deformed: True | **False**
        :param file: filename (str)
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def beamValuePlot(self, variable_name='', factor=0.0, file=None):
        """
        Create a traditional beam value plot, i.e., moment and shear diagrams.

        If **file** is given, store the plot to that file.
        Use proper file extensions to indicate the desired format (.png, .pdf)

        :param variable_name: string code for variable
        :param factor: displacement scaling factor
        :param file: filename (str)
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def addForces(self, axs):
        """
        add nodal forces to the plot shown in **axs**

        :param axs: axis on which to plot
        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)

    def set_axes_equal(self, ax):
        '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
        cubes as cubes, etc..  This is one possible solution to Matplotlib's
        ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

        :param ax: a matplotlib axis, e.g., as output from plt.gca().

        Source of this method from: https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
        '''

        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()

        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)

        # The plot bounding box is a sphere in the sense of the infinity
        # norm, hence I call half the max range the plot radius.
        plot_radius = 0.5 * max([x_range, y_range, z_range])

        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
