from .AbstractPlotter import *


class Plotter(AbstractPlotter):
    """
    class: representing a Plotter object
    """

    def __init__(self):
        super(Plotter, self).__init__()
        self.vertices  = []
        self.lines     = []

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

        :param R: list or tuple of nodal force vectors
        """
        self.reactions = R

    def displacementPlot(self, file=None, **kwargs):
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

    def valuePlot(self, deformed=False, file=None, **kwargs):
        r"""
        Create a plot using colors to identify magnitude of internal force.

        .. list-table:: known variable_name
            :header-rows: 1

            * - keyword
              - description
            * - 'ux', 'uy', 'uz'
              - displacement component in the `x`, `y` or `z` direction
            * - 'rx', 'ry', 'rz'
              - rotation around the `x`, `y` or `z` axis
            * - 'sxx', 'syy', 'szz'
              - normal stress :math:`\sigma_{xx}`, :math:`\sigma_{yy}`, :math:`\sigma_{zz}`
            * - 'sxy', 'syz', 'szx'
              - shear stress :math:`\sigma_{xy}`, :math:`\sigma_{yz}`, :math:`\sigma_{zx}`
            * - 'epsxx', 'epsyy', 'epszz'
              - normal strain :math:`\varepsilon_{xx}`, :math:`\varepsilon_{yy}`, :math:`\varepsilon_{zz}`
            * - 'epsxy', 'epsyz', 'epszx'
              - engineering shear strain :math:`\gamma_{xy}=2\varepsilon_{xy}`,
                :math:`\gamma_{yz}=2\varepsilon_{yz}`,
                :math:`\gamma_{zx}=2\varepsilon_{zx}`
            * - 'T'
              - temperature
            * - 'qx', 'qy', 'qz'
              - `x`, `y` or `z` component of the temperature gradient, :math:`q_i = \partial T/\partial x_i`


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
    plotter = Plotter()
