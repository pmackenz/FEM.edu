import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm


class Plotter():
    """
    class: representing a Plotter object
    """

    def __init__(self):
        self.vertices  = []
        self.lines     = []
        self.disp      = []
        self.values    = []
        self.reactions = []

    def __str__(self):
        return "Plotter() object"

    def __repr__(self):
        return str(self)

    def setMesh(self, vert, lines):
        """

        :param vert:
        :param lines:
        """
        self.vertices = np.array(vert)
        self.lines    = np.array(lines)

    def setDisplacements(self, disp):
        """

        :param disp:
        """
        self.disp = np.array(disp)

    def setValues(self, vals):
        """

        :param vals:
        """
        self.values = np.array(vals)

    def setReactions(self, R):
        """

        :param R: array of nodal force vectors
        """
        self.reactions = np.array(R)

    def displacementPlot(self, file=None):
        """
        Create a deformed system plot

        If **file** is given, store the plot to that file.
        Use proper file extensions to indicate the desired format (.png, .pdf)

        :param file: filename (str)
        """
        if len(self.vertices[0]) == 3:
            fig = plt.figure(figsize=(10, 10))
            axs = fig.gca(projection='3d')

            # plot the undeformed lines
            for line in self.lines:
                vert0 = self.vertices[line[0]]
                vert1 = self.vertices[line[1]]
                x = [vert0[0], vert1[0]]
                y = [vert0[1], vert1[1]]
                z = [vert0[2], vert1[2]]
                axs.plot(x, y, z, '-k', lw=2)

            # plot the deformed lines
            if len(self.disp) == len(self.vertices):
                for line in self.lines:
                    vert0 = self.vertices[line[0]].copy()
                    vert1 = self.vertices[line[1]].copy()
                    vert0 += self.disp[line[0]]
                    vert1 += self.disp[line[1]]
                    x = [vert0[0], vert1[0]]
                    y = [vert0[1], vert1[1]]
                    z = [vert0[2], vert1[2]]
                    axs.plot(x, y, z, '-r', lw=3)

            if self.reactions != []:
                self.addForces(axs)

            self.set_axes_equal(axs)
            #axs.set_aspect('equal')
            #axs.set_axis_off()

        else:
            fig, axs = plt.subplots()

            # plot the undeformed lines
            for line in self.lines:
                vert0 = self.vertices[line[0]]
                vert1 = self.vertices[line[1]]
                x = [vert0[0], vert1[0]]
                y = [vert0[1], vert1[1]]
                axs.plot(x, y, '-k', lw=2)

            # plot the deformed lines
            if len(self.disp) == len(self.vertices):
                for line in self.lines:
                    vert0 = self.vertices[line[0]].copy()
                    vert1 = self.vertices[line[1]].copy()
                    vert0 += self.disp[line[0]]
                    vert1 += self.disp[line[1]]
                    x = [vert0[0], vert1[0]]
                    y = [vert0[1], vert1[1]]
                    axs.plot(x, y, '-r', lw=3)

            if self.reactions != []:
                self.addForces(axs)

            axs.set_aspect('equal')
            axs.set_axis_off()

        plt.show()

    def valuePlot(self, deformed=False, file=None):
        """
        Create a plot using colors to identify magnitude of internal force.

        If **file** is given, store the plot to that file.
        Use proper file extensions to indicate the desired format (.png, .pdf)

        :param deformed: True | **False**
        :param file: filename (str)
        """
        fig, axs = plt.subplots()

        # plot the lines
        segments = []

        if len(self.disp) == len(self.vertices):
            for line in self.lines:
                vert0 = self.vertices[line[0]].copy()  # we need a copy since we will be modifying this in the lines below
                vert1 = self.vertices[line[1]].copy()  # we need a copy since we will be modifying this in the lines below
                if deformed:
                    vert0 += self.disp[line[0]]   # it's this += that modifies the vertices if we don't use a copy
                    vert1 += self.disp[line[1]]   # it's this += that modifies the vertices if we don't use a copy
                #x = [vert0[0], vert1[0]]
                #y = [vert0[1], vert1[1]]
                #axis.plot(x,y,'-r',lw=3)
                segments.append(np.array([vert0, vert1]))

        # Create a continuous norm to map from data points to colors
        lc = LineCollection(np.array(segments), cmap='rainbow')
        # Set the values used for colormapping
        lc.set_array(self.values)
        lc.set_linewidth(3)
        line = axs.add_collection(lc)
        fig.colorbar(line, ax=axs)

        if self.reactions != []:
            self.addForces(axs)

        axs.set_aspect('equal')
        axs.set_axis_off()

        plt.autoscale(enable=True, axis='x', tight=False)
        plt.autoscale(enable=True, axis='y', tight=False)
        plt.show()

    def addForces(self, axs):
        """
        add nodal forces to the plot shown in **axs**

        :param axs: axis on which to plot
        """
        if len(self.reactions) == len(self.vertices):
            Fx = []
            Fy = []
            X = []
            Y = []
            for (point, force) in zip(self.vertices, self.reactions):
                if np.linalg.norm(force) > 1.0e-3:
                    X.append(point[0])
                    Y.append(point[1])
                    Fx.append(-force[0])
                    Fy.append(-force[1])

            axs.quiver(X,Y, Fx, Fy, color='green')

    def set_axes_equal(self, ax):
        '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
        cubes as cubes, etc..  This is one possible solution to Matplotlib's
        ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

        Input
          ax: a matplotlib axis, e.g., as output from plt.gca().

        Taken from: https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
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



if __name__ == "__main__":
    # testing the plotter
    plotter = Plotter()
