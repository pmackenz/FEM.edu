"""
================================================================
Post-buckling of a truss column using displacement Control.
================================================================

This is the same structure as in problem :ref:`sphx_glr_auto_examples_trusses_plot_truss07.py` but
using a mix of load control and displacement control.

Author: Peter Mackenzie-Helnwein
"""

# %%
# Setup
import numpy as np
import matplotlib.pyplot as plt

from femedu.examples import Example

from femedu.domain import System, Node
from femedu.elements.finite import Truss
from femedu.materials import FiberMaterial
from femedu.solver import NewtonRaphsonSolver

# %%
# Create the example by subclassing the :py:class:`Example`

class ExampleTruss08(Example):

    # sphinx_gallery_start_ignore
    # sphinx_gallery_thumbnail_number = -1
    def docString(self):
        s = """
        Buckling of a truss column.

        This is the same structure as in problem :ref:`plot_truss07.py` but 
        using a mix of load control and displacement control.
        
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
        model.plot(factor=0., filename="truss08_undeformed.png", show_loads=1, show_bc=1, title="Undeformed System")

        #
        # performing the analysis
        #
        model.resetDisp()
        model.setLoadFactor(0.0)

        # setting target load levels
        levels = np.linspace(0.0, 2.00, 30)

        # set up data collection
        load_list = []
        data_list = []

        # reset the analysis
        model.resetDisp()

        # apply all load steps
        for lam in levels:

            model.setLoadFactor(lam)
            model.solve()

            # collect data
            load_list.append(lam)
            data_list.append(nd_11.getDisp())

            # stop load control once lateral displacement of the top node exceeds 75 mm
            if nd_11.getDisp('ux')[0] > 75 * mm:
                break

        # plot the deformed shape
        model.plot(factor=1.0,
                   title="Deformed Sytem at Transition to Displacement Control",
                   filename="truss08_deformed_end_load_control.png",
                   show_loads=False, show_reactions=False)

        #
        # switching to displacement control
        #

        # remember displacement at which we switched to displacement control
        disp_switch = nd_11.getDisp('ux')[0]

        # let's start the displacement control at the current level to verify
        # functionality.
        target = disp_switch

        while True:

            model.setDisplacementControl(nd_11, 'ux', target)
            model.solve()

            # collect data
            lam = model.loadfactor
            load_list.append(lam)
            data_list.append(nd_11.getDisp())

            # increase the target displacement by 200 mm
            target += 200 * mm

            # stop displacement control once lateral displacement of the top node exceeds 3500 mm
            if nd_11.getDisp('ux')[0] > 3500 * mm:
                break

        model.plot(factor=1.0,
                   title="Deformed Sytem at Transition back to Load Control",
                   filename="truss08_deformed_end_disp_control.png",
                   show_loads=False, show_reactions=False)

        #
        # returning to load control
        #

        # remember displacement at which we switched to displacement control
        lam_switch = model.loadfactor

        # let's start the load control at the current level to verify
        # functionality.
        lam = lam_switch

        while True:

            model.setLoadFactor(lam)
            model.solve()

            # collect data
            load_list.append(lam)
            data_list.append(nd_11.getDisp())

            # increase the target load level by $ \Delta \lambda = 0.10 $
            lam += 0.10

            # stop load control once the load level exceeds 2.0
            if lam > 2.0:
                break

        # plot the deformed shape
        model.plot(factor=1.0,
                   title="Deformed Sytem at $ \\lambda={:.2f} $".format(lam),
                   filename="truss08_deformed_end_load_2_control.png",
                   show_loads=False, show_reactions=False)

        levels = np.array(load_list)
        data   = np.array(data_list)

        plt.figure()
        plt.plot(data, levels, '--o')

        plt.plot([disp_switch, disp_switch], [0.0,1.1],'-r')
        plt.text(2.0*disp_switch, 1.05, "transition to\ndisplacement control", rotation=90.)

        plt.plot([0.5*target,1.1*target], [lam_switch, lam_switch], '-g')
        plt.text(0.5*target, 1.05 * lam_switch, "transition to\nload control" , ha='left')

        plt.grid(True)
        plt.xlabel('displacements $ u_i $')
        plt.ylabel('load factor $ \lambda $')
        plt.legend(['$ u_x $','$ u_y $'])
        plt.savefig("truss08_deformation_history.png")
        plt.show()

# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleTruss08()
    ex.run()

