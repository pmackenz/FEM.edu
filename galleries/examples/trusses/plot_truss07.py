"""
================================================================
Post-buckling of a truss column using load Control.
================================================================

We shall be using load stepping to illustrate the limitation of this control technique.

Author: Peter Mackenzie-Helnwein
"""

# %%
# Setup

from femedu.examples.Example import *

from femedu.domain.System import *
from femedu.domain.Node import *
from femedu.elements.finite.Truss import *
from femedu.materials.FiberMaterial import *
from femedu.solver.NewtonRaphsonSolver import *

# %%
# Create the example by subclassing the :py:class:`Example`

class ExampleTruss07(Example):

    # sphinx_gallery_start_ignore
    # sphinx_gallery_thumbnail_number = -2

    def docString(self):
        s = """
        Buckling of a truss column.

        We shall be using load stepping to illustrate the limitation of this control technique.
        
        Author: Peter Mackenzie-Helnwein 
        """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):
        # units
        mm = 1.
        m = 1000. * mm
        N  = 1.0
        kN = 1000. * N
        MN = 1000. * kN
        Pa = N / m**2
        MPa = 1.e6 * Pa
        GPa = 1000. * MPa

        # mesh parameters
        H = 5.00 * m
        W = H / 20

        Px = 0.05 * kN

        #Ny = 11
        Ny = 20

        params_vert = dict(
            E = 200 * GPa,
            A = 10 * mm**2
        )

        params_diag = dict(
            E = 200 * GPa,
            A = 50 * mm**2
        )

        params_brace = dict(
            E = 200 * GPa,
            A = 50 * mm**2
        )

        EI = params_vert['E'] * params_vert['A'] * W**2 / 2
        Py = EI * np.pi**2 / (2*H)**2

        # initialize a system model
        model = System()
        model.setSolver(NewtonRaphsonSolver())

        # create floor nodes
        h = 0.0
        nd0 = Node(0.0, h)
        nd1 = Node(  W, h)

        nd_00 = nd0
        nd_01 = nd1

        model.addNode(nd0, nd1)

        for layer in range(Ny):
            h += H / Ny
            nd_10 = Node(0.0, h)
            nd_11 = Node(  W, h)

            model.addNode(nd_10, nd_11)

            # create elements
            model += Truss(nd_00, nd_10, FiberMaterial(params_vert))
            model += Truss(nd_01, nd_11, FiberMaterial(params_vert))
            model += Truss(nd_00, nd_11, FiberMaterial(params_diag))
            model += Truss(nd_10, nd_11, FiberMaterial(params_brace))

            # prep for the next level
            nd_00 = nd_10
            nd_01 = nd_11

        # apply boundary conditions
        nd0.fixDOF(['ux','uy'])
        nd1.fixDOF(['ux','uy'])

        # build reference load
        nd_10.addLoad([-Py/2],['uy'])
        nd_11.addLoad([-Py/2],['uy'])
        nd_10.addLoad([ Px/2],['ux'])
        nd_11.addLoad([ Px/2],['ux'])

        # write out report
        model.report()

        # create plots
        model.setLoadFactor(1.0)    # needed to show the reference load
        model.plot(factor=0., filename="truss07_undeformed.png", show_loads=1, show_bc=1, title="Undeformed System")

        #
        # performing the analysis
        #
        model.resetDisp()

        # setting target load levels
        levels = np.linspace(0.0, 1.50, 15)
        levels = np.linspace(0.0, 2.00, 30)

        # set up data collection
        data_list = []

        # reset the analysis
        model.resetDisp()

        # apply all load steps
        for lam in levels:

            model.setLoadFactor(lam)
            model.solve()

            # collect data
            data_list.append(nd_11.getDisp())

        # plot the deformed shape
        model.plot(factor=1.0, filename="truss07_deformed.png", show_loads=False, show_reactions=False)

        data = np.array(data_list)

        plt.figure()
        plt.plot(data, levels, '--o')
        plt.grid(True)
        plt.xlabel('displacements $ u_i $')
        plt.ylabel('load factor $ \lambda $')
        plt.legend(['$ u_x $','$ u_y $'])
        plt.savefig("truss07_deformation_history.png")
        plt.show()

# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleTruss07()
    ex.run()

