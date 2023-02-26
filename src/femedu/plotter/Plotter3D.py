from .Plotter import *

class Plotter3d(Plotter):
    """
    class: representing a Plotter object
    """

    def __init__(self):
        super(Plotter3d, self).__init__()

    def displacementPlot(self, file=None):
        """
        Create a deformed 3d system plot

        If **file** is given, store the plot to that file.
        Use proper file extensions to indicate the desired format (.png, .pdf)

        :param file: filename (str)
        """
        fig, axs = plt.subplots()

        # plot the undeformed lines
        for line in self.lines:
            vert0 = self.vertices[line[0]]
            vert1 = self.vertices[line[1]]
            x = [vert0[0], vert1[0]]
            y = [vert0[1], vert1[1]]
            z = [vert0[2], vert1[2]]
            axs.plot(x,y,z,'-k',lw=2)

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
                axs.plot(x,y,z,'-r',lw=3)

        if self.reactions != []:
            self.addForces(axs)

        #axs.set_aspect('equal')
        #axs.set_axis_off()

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




if __name__ == "__main__":
    # testing the plotter
    plotter = Plotter3d()
