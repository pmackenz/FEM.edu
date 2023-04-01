"""
======================================================
Buckling of a braced building frame
======================================================

modeled using a 2D frame element

.. list-table:: setting given parameters

    * - N  = 2
      - number of elements
    * - L  = 100.0
      - column length
    * - EA = 2000000.0
      - axial stiffness
    * - EI = 21000.0
      - flexural stiffness
    * - w  = 0.1
      - applied lateral load

Author: Peter Mackenzie-Helnwein
"""
import matplotlib.pyplot as plt
from femedu.examples.Example import *

from femedu.domain import *
from femedu.solver.NewtonRaphsonSolver import *
from femedu.elements.Truss import *
from femedu.elements.Frame2D import *
from femedu.materials.FiberMaterial import *
from femedu.materials.ElasticSection import *


class ExampleMixed01(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = """
    Buckling of a building frame

    degrees of freedom:
        ux ... horizontal displacement
        uy ... vertical displacement
        rz ... rotation, theta
        
    Author: Peter Mackenzie-Helnwein 
    """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):
        # initialize a system model

        N  = 8     # number of elements

        B = 240.
        H = 200.

        E  = 29000.0     # steel MOE

        A = 20.0         # frame area
        I = 10.0         # frame moment of inertia
        Ab = 1.0         # brace area

        w = 1.0

        Ph = 0.01      # additional horizontal load per floor
        Ph = 0.10      # additional horizontal load per floor
        Ph = 1.00      # additional horizontal load per floor
        Ph = 0.00      # additional horizontal load per floor

        # ========== setting global parameters ==============

        target_load_level = 27
        max_steps = 10
        load_levels = np.linspace(0, target_load_level, max_steps)

        # ========= build your structural model =============

        model = System()
        model.setSolver(NewtonRaphsonSolver())

        x0 = 0.0
        x1 = B / 3
        x2 = 2 * B / 3
        x3 = B

        y0 = 0.0
        y1 = H / 4
        y2 = 2 * H / 4
        y3 = 3 * H / 4
        y4 = H

        X10 = Node(x0, y0)
        X11 = Node(x0, y1)

        X20 = Node(x1, y0)
        X21 = Node(x1, y1)

        model.addNode(X10,X11)
        model.addNode(X20,X21)

        # columns

        frameParams = {'E': E, 'A': A, 'I': I}
        C11 = Frame2D(X10, X11, ElasticSection(frameParams))
        C21 = Frame2D(X20, X21, ElasticSection(frameParams))
        model.addElement(C11,C21)

        # floors

        params = {'E': E, 'A': A, 'I': 8*I}
        F11 = Frame2D(X11, X21, ElasticSection(params))
        model.addElement(F11)

        # braces

        braceParams = {'E': E, 'A': Ab}
        model += Truss(X10, X21, FiberMaterial(braceParams))

        # fixities
        X10.fixDOF('ux','uy','rz')   # fixed
        X20.fixDOF('ux','uy','rz')   # fixed

        # reference load
        #Pcr = np.pi**2 * EI / L**2
        model.resetLoad()            # size load vector and initialize
        #model.addLoad(Xn, -Pcr, dof=0) # add a horizontal force (first dof only) ; remember C-style indexing: 0,1,...,(n-1)

        # floor loading as distributed loads ...
        F11.setDistLoad(-w)

        # wind load ...
        X11.addLoad([Ph],['ux'])   # horizontal load

        # show model information
        model.report()

        print("\n==== perform the analysis ===\n")

        # * apply the load in multiple smaller load steps

        # set up data recorder
        model.initRecorder()
        model.trackStability(True)

        # initialize the analysis:
        model.resetDisplacements()   # set U to all zeros
        model.setLoadFactor(0.0)     # define a known equilibrium solution

        model.startRecorder()

        detKt   = []
        lambdas = []

        # solve for all load_levels
        for loadfactor in load_levels:

            # define node X2 as the controled node; downward direction is prescribed:
            model.setLoadFactor(loadfactor)
            model.solve(verbose=True)

            # stability check
            lambdas.append(model.loadfactor)
            detKt.append(model.solver.checkStability())

            # report results
            print('+')
            #model.report()

            print("\n=== next load level ===\n")

        #
        # ==== create some nice plots ===
        #

        model.report()

        model.plot(factor=1.0, filename="frame3_deformed.png")

        fig, ax = plt.subplots()

        ax.plot(lambdas,detKt,'--*r')
        ax.grid(True)
        ax.set_xlabel('Load factor, $ \lambda $')
        ax.set_ylabel("Stability index, $ {det}\: {\\bf K}_t $")

        fig.savefig("frame3_stability.png")
        fig.show()

        model.beamValuePlot("F", filename="frame3_force.png")
        model.beamValuePlot("V", filename="frame3_shear.png")
        model.beamValuePlot("M", filename="frame3_moment.png")

        model.plotBucklingMode(factor=25., filename="frame3_buckling_mode0.png")


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleMixed01()
    ex.run()



