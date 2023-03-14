"""
Example: plate04

PatchMesher test for the LinearTriangle (Constant Strain Triangle)

"""
from ...examples.Example import *

from ...domain import *
from ...solver.NewtonRaphsonSolver import *
from ...elements.LinearTriangle import *
from ...materials.PlaneStress import *
from ...mesher import *


class ExamplePlate04(Example):

    def docString(self):
        s = """
    ## Patch test for triangular plate under in-plane loading

    The patch test is an empirical minimum test which every finite element has to pass to ensure convergence with mesh refinement.

    It consists of a problem for which a known homogeneous solution exists.  For plates, we commonly use a rectangular plate subject to homogenous edge loading, e.g., constant tension in the x-direction, or constant shear, etc.

    The mesh must contain distorted elements and at least one element not attached to any node on the boundary.
        
    Author: Peter Mackenzie-Helnwein 
    """
        return s

    def problem(self):
        # ========== setting mesh parameters ==============

        Nx = 6        # number of elements in the mesh
        Ny = 5        # number of elements in the mesh
        Lx = 100.0    # length of plate in the x-direction
        Ly =  80.0    # length of plate in the y-direction

        # ========== setting material parameters ==============

        params = dict(
            E  = 20000.,    # Young's modulus
            nu = 0.250,     # Poisson's ratio
            t  = 1.00       # thickness of the plate
        )

        # ========== setting load parameters ==============

        px  = 10.0         # uniform load normal to x=const
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
        model.setSolver(NewtonRaphsonSolver())

        # create nodes

        mesher1 = PatchMesher(model,
                             (0,0), (Lx,0), (Lx,Ly), (0,Ly),                                     # corner nodes
                             (Lx/2, -0.05*Ly), (1.2*Lx, Ly/2), (Lx/2, 0.90*Ly), (0.05*Lx, Ly/2), # mid-side nodes
                             (0.55*Lx, 0.45*Ly))                                                 # center node
        nodes1, elements1 = mesher1.triangleMesh(Nx, Ny, LinearTriangle, PlaneStress(params))

        mesher2 = TriPatchMesher(model,
                                 (0,0), (Lx,0), (Lx/2,Ly),                                       # corner nodes
                                 (Lx/2,0.15*Ly), (0.85*Lx,0.8*Ly), (0.20*Lx,0.6*Ly) # mid-side nodes
                                 )
        mesher2.shift(1.25*Lx, Ly/2)
        #nodes2, elements2 = mesher2.triangleMesh(Nx, LinearTriangle, PlaneStress(params))
        nodes2, elements2 = mesher2.quadMesh(Nx, LinearTriangle, PlaneStress(params))

        nodes    = nodes1    + nodes2
        elements = elements1 + elements2

        # define support(s)

        fix_x = (0,)
        fix_y = (0,4)

        for idx in fix_x:
            nodes[idx].fixDOF('ux')    # horizontal support left end
        for idx in fix_y:
            nodes[idx].fixDOF('uy')          # vertical support right end

        # ==== complete the reference load ====

        # surface loads on the left side
        elements[ 0].setSurfaceLoad(2,px)
        elements[8].setSurfaceLoad(2,px)
        elements[16].setSurfaceLoad(2,px)

        # surface loads on the right side
        elements[ 7].setSurfaceLoad(2,px)
        elements[15].setSurfaceLoad(2,px)
        elements[23].setSurfaceLoad(2,px)

        # these are only nodal forces as part of the reference load
        # .. load only the upper node
        #print(model)

        model.report()

        model.plot(factor=10., title="undeformed system", filename="plate04_undeformed.png")

        model.setLoadFactor(10.0)
        model.solve()

        model.report()

        model.plot(factor=10., filename="plate04_deformed.png")

