"""
==========================================================
Pulling a plate with a circular hole
==========================================================

Using PatchMesher to model a quarter of the plate

"""
import math
import sys
import numpy as np

from femedu.examples.Example import *

from femedu.domain import *
from femedu.solver.NewtonRaphsonSolver import *
from femedu.elements.linear.Quad import *
from femedu.materials.PlaneStress import *
from femedu.mesher import *


class ExamplePlate06(Example):

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

        Nx = 10        # number of elements in the mesh
        Ny = 8        # number of elements in the mesh
        Lx = 120.0    # length of plate in the x-direction
        Ly =  80.0    # length of plate in the y-direction
        R  = Ly / 2.

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

        #  5--------------6--------------7
        #  |             / \             |
        #  |            /   \            |
        #  |           /     \           |
        #  3----8----4        \          |
        #             \        \         |
        #              9        \        |
        #               \        \       |
        #                0--------1------2

        pts = (
            ( R, 0),
            (0.5*(R+Lx), 0),
            (Lx, 0),
            (0, R),
            (R*np.cos(np.radians(45.)), R*np.sin(np.radians(45.))),
            (0, Ly),
            (0.5*Lx, Ly),
            (Lx, Ly),
            (R*np.cos(np.radians(67.5)), R*np.sin(np.radians(67.5))),
            (R*np.cos(np.radians(22.5)), R*np.sin(np.radians(22.5))),
        )

        mesher1 = PatchMesher(model, pts[3], pts[4], pts[6], pts[5], pts[8])
        nodes1, elements1 = mesher1.quadMesh(Nx, Ny, Quad, PlaneStress(params))

        mesher2 = PatchMesher(model, pts[0], pts[1], pts[6], pts[4], None, None, None, pts[9])                                               # center node
        nodes2, elements2 = mesher2.quadMesh(Ny, Nx, Quad, PlaneStress(params))

        mesher3 = PatchMesher(model, pts[1], pts[2], pts[7], pts[6])
        nodes3, elements3 = mesher3.quadMesh(Nx, Nx, Quad, PlaneStress(params))

        mesher1.tie(mesher2)
        mesher2.tie(mesher3)

        nodes    = nodes1    + nodes2    + nodes3
        elements = elements1 + elements2 + elements3

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

        model.report()

        model.plot(factor=0, title="undeformed system", filename="plate05_undeformed.png", show_bc=1, show_loads=1)

        model.setLoadFactor(10.0)
        model.solve()

        model.report()

        model.plot(factor=10., filename="plate05_deformed.png", show_bc=1, show_loads=1, show_reactions=1)


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExamplePlate06()
    ex.run()


