"""
======================================================
Buckling of a building frame
======================================================

modeled using
* a 2D frame element
* inclined boundary condition (using Transformation)
* element loading

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
import numpy as np

from femedu.examples.Example import *

from femedu.domain import *
from femedu.solver.NewtonRaphsonSolver import *
from femedu.solver.LinearSolver import *
from femedu.elements.linear.Frame2D import *
from femedu.domain.Frame2dTransformation import Transformation
from femedu.materials.ElasticSection import *


class ExampleFrame04(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = """
    Modeled a plane frame system using 
    * a 2D frame element
    * inclined boundary condition (using Transformation)
    * element loading

    degrees of freedom:
        ux ... horizontal displacement
        uy ... vertical displacement
        rz ... rotation, theta
        
    Author: Peter Mackenzie-Helnwein 
    """
        return s

    # sphinx_gallery_end_ignore

    def createMesh(self, nelem=4):
        """
        Creates a mesh using `nelem` elements for each of the straight segments.
        """
        if (not isinstance(nelem, int)) or nelem % 4:
            msg = "nelem must be an integer multiple of 4 (4,8,12,...)"
            raise TypeError(msg)

        # units
        inch = 1
        ft = 12 * inch
        kip = 1
        kips = kip
        lb = kip / 1000
        psi = lb / inch ** 2
        ksi = kip / inch ** 2
        degrees = np.pi / 180.

        # problem parameters
        MOE =  29000 * ksi
        EI  = 162133 * kips * inch ** 2
        EA  =  30400 * kips
        L   =     10 * ft
        w0  =  0.333 * kip / ft

        s = np.linspace(0.0, 1.0, nelem // 2 + 1)

        params = dict(
            E = MOE,
            A = EA/MOE,
            I = EI/MOE,
        )

        material = ElasticSection(params)

        # nodes
        nodes  = [ Node(0.0, si * L) for si in s]
        nodes += [ Node(si * L * np.cos(np.radians(30.)), L + si * L * np.sin(np.radians(30.))) for si in s[1:]]

        # elements
        elements = [Frame2D(nodes[i], nodes[i+1], material) for i in range(nelem)]

        # fixities
        # ... the first node
        nodes[0].fixDOF(['ux','uy','rz'])
        # ... the last node
        nvec = nodes[-1].getPos() - nodes[-2].getPos()   # vector parallel to the member axis
        svec = np.array([[0, -1],[1, 0]]) @ nvec       # vector perpendicular to the member axis

        try:
            transform = Frame2dTransformation(nvec, svec)  # an in-plane rotation
            nodes[-1].addTransformation(transform)         # defining a local frame for the last node
        except:
            print("no transformation class found")
        nodes[-1].fixDOF(['uy', ])                     # fixing the LOCAL y-direction

        # load the top half of the vertical member
        # and the first half of the inclined member
        for i in range(nelem//4, 3*nelem//4):
            elements[i].setDistLoad(-w0)

        return (nodes, elements)

    def problem(self):
        # initialize a system model

        N  = 16     # number of elements

        # ========== setting global parameters ==============

        target_load_level = 1.

        # ========= build your structural model =============

        model = System()
        model.setSolver(NewtonRaphsonSolver())
        # model.setSolver(LinearSolver())

        nodes, elements = self.createMesh(N)

        model.addNode(*nodes)
        model.addElement(*elements)

        # show model information
        print(model)

        print("\n==== perform the analysis ===\n")

        # solve
        model.setLoadFactor(target_load_level)
        model.solve(verbose=True, max_steps=10)

        #
        # ==== create some nice plots ===
        #

        model.report()

        model.plot(factor=10.0, filename="frame4_deformed.png", show_bc=1)

        model.beamValuePlot("F", filename="frame4_force.png")
        model.beamValuePlot("V", filename="frame4_shear.png")
        model.beamValuePlot("M", filename="frame4_moment.png")


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleFrame04()
    ex.run()


