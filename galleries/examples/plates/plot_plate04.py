"""
==========================================================
Patch test for triangular plate under in-plane loading
==========================================================

PatchMesher test for the LinearTriangle (Constant Strain Triangle)

"""
import math

from femedu.examples.Example import *

from femedu.domain import *
from femedu.solver.NewtonRaphsonSolver import *
from femedu.elements.linear.Triangle import *
from femedu.elements.linear.Quad import *
from femedu.materials.PlaneStress import *
from femedu.mesher import *


class ExamplePlate04(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = """
    ## Patch test for triangular plate under in-plane loading

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

        Nx = 7        # number of elements in the mesh
        Ny = 6        # number of elements in the mesh
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

        # create reference points
        pt0 = (0,0)
        pt1 = (0.5*Lx,0)
        pt2 = (Lx,Ly)
        pt3 = (0,Ly)
        pt4 = (0.25*Lx, 0.0)
        pt5 = (0.65*Lx, Ly/2)
        pt6 = (0.55*Lx, Ly)
        pt7 = (0.0, Ly/2)
        pt8 = (0.325*Lx, 0.55*Ly)
        pt9 = (Lx,0.0)

        mesher1 = PatchMesher(model,
                             pt0, pt1, pt2, pt3, # corner nodes
                             pt4, pt5, pt6, pt7, # mid-side nodes
                             pt8)                # center node
        #nodes1, elements1 = mesher1.triangleMesh(Nx, Ny, LinearTriangle, PlaneStress(params))
        nodes1, elements1 = mesher1.quadMesh(Nx, Ny, Quad, PlaneStress(params))

        mesher2 = TriPatchMesher(model,
                                 pt1, pt9, pt2, # corner nodes
                                 None, None, pt5, # mid-side nodes
                                 )
        #mesher2.shift(1.25*Lx, Ly/2)
        nodes2, elements2 = mesher2.triangleMesh(Ny, Triangle, PlaneStress(params))
        #nodes2, elements2 = mesher2.quadMesh(Ny, Triangle, PlaneStress(params))

        # tie the patches together
        mesher1.tie(mesher2)

        nodes    = nodes1    + nodes2
        elements = elements1 + elements2

        # define support(s)

        ## find nodes at y==0 and x==0

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[0], 0.0):
                node.fixDOF('ux')    # horizontal support left side
            if math.isclose(X[1], 0.0):
                node.fixDOF('uy')    # vertical support at y==0

        # ==== complete the reference load ====

        Xo = np.array([Lx, 0.0])
        No = np.array([1.0, 0.0])

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[0],Lx):
                print(node)
                for elem in node.elements:
                    print('+', elem)
                    for face in elem.faces:
                        for x, area in zip(face.pos, face.area):
                            if np.abs( (x - Xo) @ No ) < 1.0e-2 and  No @ area / np.linalg.norm(area):
                                face.setLoad(px, 0.0)


        #model.report()

        model.plot(factor=0, title="undeformed system", filename="plate04_undeformed.png", show_bc=1, show_loads=1)

        model.setLoadFactor(10.0)
        model.solve()

        #model.report()

        model.plot(factor=10., filename="plate04_deformed.png")


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExamplePlate04()
    ex.run()


