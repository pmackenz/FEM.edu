"""
=========================================================================================================
Patch test for mixed mesh of triangular and quadrilateral plate elements under in-plane loading
=========================================================================================================

PatchMesher test for the mixed mesh of triangular and quadrilateral plate elements.

"""
import math
import numpy as np

from femedu.examples import Example

from femedu.domain import System, Node
from femedu.solver import NewtonRaphsonSolver
#from femedu.elements.linear import Quad, Quad9
from femedu.elements.finite import Quad, Quad9
from femedu.materials import PlaneStress
from femedu.mesher import *


class ExamplePlate09b(Example):

    # sphinx_gallery_start_ignore
    # sphinx_gallery_thumbnail_number = 2
    def docString(self):
        s = """
    ## Patch test for mixed mesh triangular and quadrilateral plate under in-plane loading

    The patch test is an empirical minimum test which every finite element has to pass to ensure convergence with mesh refinement.

    It consists of a problem for which a known homogeneous solution exists.  
    For plates, we commonly use a rectangular plate subject to homogeneous edge loading, 
    e.g., constant tension in the x-direction, or constant shear, etc.

    The mesh must contain distorted elements and at least one element not attached to any node on the boundary.

    Author: Peter Mackenzie-Helnwein 
    """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):
        # ========== setting mesh parameters ==============

        Nx = 3  # number of elements in the mesh
        Ny = 3  # number of elements in the mesh
        Lx = 100.0  # length of plate in the x-direction
        Ly = 80.0  # length of plate in the y-direction

        # ========== setting material parameters ==============

        params = dict(
            E=20000.,  # Young's modulus
            nu=0.250,  # Poisson's ratio
            t=1.00  # thickness of the plate
        )

        # ========== setting load parameters ==============

        px = 10.0  # uniform load normal to x=const
        py = 0.0  # uniform load normal to y=const
        pxy = 0.0  # uniform shear load on x=const and y=const

        # ========== setting analysis parameters ==============

        target_load_level = 1.00  # reference load
        max_steps = 1  # number of load steps: 1 -> [0.0, 1.0]

        # define a list of target load levels
        load_levels = np.linspace(0, target_load_level, max_steps+1)

        #
        # ==== Build the system model ====
        #

        model = System()
        model.setSolver(NewtonRaphsonSolver())

        # create reference points
        pt0 = (0, 0)
        pt1 = (Lx, 0.0)
        pt2 = (Lx, Ly)
        pt3 = (0, Ly)
        pt4 = (0.4 * Lx, 0)
        pt5 = (Lx, 0.5 * Ly)
        pt6 = (0.5 * Lx, Ly)
        pt7 = (0.0, 0.4 * Ly)
        pt8 = (0.6 * Lx, 0.55 * Ly)

        mesher = PatchMesher(model,
                              pt0, pt1, pt2, pt3,  # corner nodes
                              pt4, pt5, pt6, pt7,  # mid-side nodes
                              pt8)  # center node
        #nodes, elements = mesher.quadMesh(Nx, Ny, Quad, PlaneStress(params))
        nodes, elements = mesher.quadMesh(Nx, Ny, Quad9, PlaneStress(params))

        # define support(s)

        ## find nodes at y==0 and x==0

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[0], 0.0):
                node.fixDOF('ux')  # horizontal support left side
            if math.isclose(X[1], 0.0):
                node.fixDOF('uy')  # vertical support at y==0

        # ==== complete the reference load ====

        for elem, face in model.findFacesAlongLine((Lx, 0.0), (0.0, 1.0)):
            face.setLoad(px, 0.0)

        # model.report()

        model.plot(factor=0, title="undeformed system", show_bc=1, show_loads=1)

        model.setLoadFactor(1.0)
        model.solve(verbose=1)

        model.report()

        model.plot(factor=25.)
        model.valuePlot('sxx', show_mesh=1)
        model.valuePlot('syy', show_mesh=1)
        model.valuePlot('sxy', show_mesh=1)


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExamplePlate09b()
    ex.run()


