"""
==========================================================
Profiling the code
==========================================================

Using PatchMesher to model a quarter of the plate

"""

# sphinx_gallery_thumbnail_path = '_static/profiler.png'

SPARSE = False

import math
import sys
import numpy as np

from femedu.examples import Example

from femedu.domain import System, Node
from femedu.solver import NewtonRaphsonSolver, SparseSolver
from femedu.elements.linear import Triangle
from femedu.materials import PlaneStress
from femedu.mesher import *


class ExamplePlate14(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = """
    ## Pulling a plate with a circular hole
    
    The full plate has a dimension of (2*Lx)*(2*Ly) , Ly < Lx, 
    and a circular hole with radius R < Ly.

    Using PatchMesher to model a quarter of the plate,
    applying symmetry conditions, and face loading on the 
    free surface x=a.
        
    Author: Peter Mackenzie-Helnwein 
    """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):
        # ========== setting mesh parameters ==============

        Nx = 60        # number of elements in the mesh
        Ny = 40        # number of elements in the mesh

        Lx = 120.0    # length of plate in the x-direction
        Ly =  80.0    # length of plate in the y-direction
        R  = Ly / 2.

        # ========== setting material parameters ==============

        params = dict(
            E  = 200.,      # Young's modulus
            nu = 0.450,     # Poisson's ratio
            t  = 1.00       # thickness of the plate
        )

        # ========== setting load parameters ==============

        px  = 20.0         # uniform load normal to x=const
        py  =  0.0         # uniform load normal to y=const
        pxy =  0.0         # uniform shear load on x=const and y=const

        # ========== setting analysis parameters ==============

        target_load_level = 1.00     # reference load
        max_steps = 2                # number of load steps: 2 -> [0.0, 1.0]

        # define a list of target load levels
        load_levels = np.linspace(0, target_load_level, max_steps)

        #
        # ==== Build the system model ====
        #

        model = System()
        if SPARSE:
            model.setSolver(SparseSolver())
        else:
            model.setSolver(NewtonRaphsonSolver())

        # create nodes

        #  3---------2
        #  |         |
        #  |         |
        #  |         |
        #  1---------1

        pts = (
            ( 0, 0),
            (Lx, 0),
            (Lx, Ly),
            (0, Ly),
        )

        mesher = PatchMesher(model, pts[0], pts[1], pts[2], pts[3])
        nodes, elements = mesher.triangleMesh(Nx, Ny, Triangle, PlaneStress(params))

        # define support(s)

        ## find nodes at y==0 and x==0

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[0], 0.0):
                node.fixDOF('ux')    # horizontal support left side
            if math.isclose(X[1], 0.0):
                node.fixDOF('uy')    # vertical support at y==0

        # ==== build the reference load ====

        # the section at the right end
        for elem, face in model.findFacesAlongLine((Lx, 0.0), (0.0, 1.0), orientation=+1):
            print('+', elem.getID(), face.id, ":", face.area)
            face.setLoad(px, pxy)

        model.plot(factor=0, title="undeformed system", filename="plate14_undeformed.png", show_bc=1, show_loads=1)

        model.setLoadFactor(1.0)
        model.solve()

        #model.plot(factor=1., filename="plate14_deformed.png")

        #model.solver.showKt(filename="plate14_spy_Kt.png")
        #np.save("plate14_Kt.npy",model.solver.Kt)


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#
# This time, we are running the example in the profiler, writing profiling data to file.
#

# sphinx_gallery_start_ignore
if __name__ == "__main__":
# sphinx_gallery_end_ignore

    import cProfile

    ex = ExamplePlate14()

    if SPARSE:
        cProfile.run('ex.run()','profile_data_sparse.txt')
    else:
        cProfile.run('ex.run()','profile_data_full.txt')

# %%
# Now it's time to process the profiling data
#

# sphinx_gallery_start_ignore
if __name__ == "__main__":
# sphinx_gallery_end_ignore

    import pstats
    from pstats import SortKey

    if SPARSE:
        p = pstats.Stats('profile_data_sparse.txt')
        p.strip_dirs() #.sort_stats(-1).print_stats()
        p.sort_stats(SortKey.NAME)
        #p.print_stats()

        p.sort_stats(SortKey.CUMULATIVE).print_stats(20)
        p.sort_stats(SortKey.TIME).print_stats(20)

    else:
        p = pstats.Stats('profile_data_full.txt')
        p.strip_dirs() #.sort_stats(-1).print_stats()
        p.sort_stats(SortKey.NAME)
        #p.print_stats()

        p.sort_stats(SortKey.CUMULATIVE).print_stats(20)
        p.sort_stats(SortKey.TIME).print_stats(20)

