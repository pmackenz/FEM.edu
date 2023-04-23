"""
==========================
Simple triangular truss.
==========================

Study of snap-through behavior using finite deformation truss elements.

We shall be using displacement control to trace the unstable portion of the static equilibrium path.

Author: Peter Mackenzie-Helnwein
"""

# sphinx_gallery_start_ignore
# sphinx_gallery_thumbnail_number = -2
# sphinx_gallery_end_ignore

import numpy as np

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

class ExampleTruss06(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = """
        Simple triangular truss. 
        
        Study of snap-through behavior using finite deformation truss elements.

        We shall be using displacement control to trace the unstable portion of the static equilibrium path.
        
        Author: Peter Mackenzie-Helnwein 
        """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):
        # initialize a system model
        model = System()
        model.setSolver(NewtonRaphsonSolver())

        # create notes
        x1=Node(0.0,0.0)
        x2=Node(5.5,0.5)
        x3=Node(9.5,0.0)

        model.addNode(x1,x2,x3)

        params = dict(
            E = 2100.,   # MOE
            A = 1.       # cross section area
        )

        # create elements
        elemA = Truss(x1,x2, FiberMaterial(params))
        elemB = Truss(x3,x2, FiberMaterial(params))

        model += elemA
        model += elemB

        # apply boundary conditions
        x1.fixDOF(['ux','uy'])
        x3.fixDOF(['ux','uy'])

        # build reference load
        x2.addLoad([-1.],['uy'])

        # write out report
        model.report()

        # create plots
        model.plot(factor=1., filename="truss05_deformed.png")

        #
        # performing the analysis
        #
        model.resetDisp()

        # setting target displaement levels
        disps = np.linspace(0.0, 1.1, 24)

        # set up data collection
        load_list = []   # will hold load factors
        data_list = []   # will hold displacements

        # reset the analysis
        model.resetDisp()
        model.setLoadFactor(0.0)

        # apply all load steps
        for u_bar in disps:

            #model.setLoadFactor(lam)
            model.setDisplacementControl(x2, 'uy', -u_bar)
            model.solve()

            # collect data
            load_list.append(model.loadfactor)
            data_list.append(x2.getDisp())

            # plot the deformed shape
            model.plot(factor=1.0, show_loads=False, show_reactions=False)

        load = np.array(load_list)
        data = np.array(data_list)

        plt.figure()
        plt.plot(data, load)
        plt.grid(True)
        plt.xlabel('displacements $ u_i $')
        plt.ylabel('load factor $ \lambda $')
        plt.legend(['$ u_x $','$ u_x $'])
        plt.show()

        #plt.figure()
        fig, (ax0,ax1) = plt.subplots(1,2)

        ax0.plot(load, data[:,0])
        ax0.grid(True)
        ax0.set_xlabel('load factor $ \lambda $')
        ax0.set_ylabel('displacements $ u_i $')
        ax0.legend(['$ u_x $'])

        ax1.plot(load, data[:,1])
        ax1.grid(True)
        ax1.set_xlabel('load factor $ \lambda $')
        ax1.set_ylabel('displacements $ u_i $')
        ax1.legend(['$ u_y $'])

        plt.show()

# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleTruss06()
    ex.run()

