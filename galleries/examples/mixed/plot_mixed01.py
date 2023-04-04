"""
======================================================
Braced building frame
======================================================

modeled using a 2D frame element for the main structure and a truss element for the brace


Author: Peter Mackenzie-Helnwein
"""
import matplotlib.pyplot as plt
from femedu.examples.Example import *

from femedu.domain import *
from femedu.solver.NewtonRaphsonSolver import *
from femedu.elements.finite.Truss import *
from femedu.elements.finite.Frame2D import *
from femedu.materials.FiberMaterial import *
from femedu.materials.ElasticSection import *


class ExampleMixed01(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = """
    Braced building frame

    degrees of freedom:
        ux ... horizontal displacement
        uy ... vertical displacement
        rz ... rotation, theta
        
    The brace has no flexural stiffness and only uses 'ux' and 'uy'
        
    Author: Peter Mackenzie-Helnwein 
    """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):
        # initialize a system model

        B = 80.
        H = 50.

        E  = 29000.0   # steel MOE

        A = 20.0       # frame area
        I = 10.0       # frame moment of inertia
        Ab = 1.0       # brace area

        w = 0.50       # uniform load on floor beam

        Ph = 20.00      # additional horizontal load per floor

        # ========== setting global parameters ==============

        target_load_level = 10
        max_steps = 10
        load_levels = np.linspace(0, target_load_level, max_steps)

        # ========= build your structural model =============

        model = System()
        model.setSolver(NewtonRaphsonSolver())

        X10 = Node(0.0, 0.0)
        X11 = Node(0.0, H)

        X20 = Node(B, 0.0)
        X21 = Node(B, H)

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

        # initialize the analysis:
        model.resetDisplacements()   # set U to all zeros
        model.setLoadFactor(0.0)     # define a known equilibrium solution

        model.plot(factor=0.0, title="undeformed system", filename="mixed01_undeformed.png")

        model.startRecorder()

        # solve for all load_levels
        for loadfactor in load_levels:

            # define node X2 as the controled node; downward direction is prescribed:
            model.setLoadFactor(loadfactor)
            model.solve(verbose=True)

            model.recordThisStep()

            print("\n=== next load level ===\n")

        #
        # ==== create some nice plots ===
        #

        model.report()

        model.plot(factor=10.0, filename="mixed01_deformed.png")
        model.beamValuePlot("F", filename="mixed01_force.png")
        model.beamValuePlot("V", filename="mixed01_shear.png")
        model.beamValuePlot("M", filename="mixed01_moment.png")


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleMixed01()
    ex.run()



