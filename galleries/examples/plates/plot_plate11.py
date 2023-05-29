"""
==========================================================
Bending a cantilever beam
==========================================================

Using PatchMesher to model the beam

"""
import math
import sys
import numpy as np

from femedu.examples import Example

from femedu.domain import System, Node
from femedu.solver import NewtonRaphsonSolver
#from femedu.elements.linear import Quad
from femedu.elements.finite import Quad
from femedu.materials import PlaneStress
from femedu.mesher import *


class ExamplePlate11(Example):

    # sphinx_gallery_start_ignore
    # sphinx_gallery_thumbnail_number = 2
    def docString(self):
        s = """
    ## Bending a cantilever beam

    Using PatchMesher to model the beam
        
    Author: Peter Mackenzie-Helnwein 
    """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):
        # ========== setting mesh parameters ==============

        Nx = 24       # number of elements in the mesh
        Ny = 8        # number of elements in the mesh
        Lx = 120.0    # length of plate in the x-direction
        Ly =  20.0    # length of plate in the y-direction

        # ========== setting material parameters ==============

        params = dict(
            E  = 20000.,    # Young's modulus
            nu = 0.250,     # Poisson's ratio
            t  = 1.00       # thickness of the plate
        )

        # ========== setting load parameters ==============

        px  =  0.0         # uniform load normal to x=Lx
        py  =  0.0         # uniform load normal to y=Ly
        pxy =  1.5         # uniform shear load on x=L

        # ========== setting analysis parameters ==============

        target_load_level = 100.00    # reference load
        max_steps = 21                # number of load steps: 2 -> [0.0, 1.0]

        # target_load_level = 10.00    # reference load
        # max_steps = 3                # number of load steps: 2 -> [0.0, 1.0]

        # define a list of target load levels
        load_levels = np.linspace(0, target_load_level, max_steps)

        #
        # ==== Build the system model ====
        #

        model = System()
        model.setSolver(NewtonRaphsonSolver())

        # create nodes

        mesher = PatchMesher(model, (0.,0.),(Lx,0.),(Lx,Ly),(0.,Ly) )
        nodes, elements = mesher.quadMesh(Nx, Ny, Quad, PlaneStress(params))

        # define support(s)

        ## find nodes at y==0 and x==0

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[0], 0.0):
                node.fixDOF('ux','uy')    # fix left side

        # ==== complete the reference load ====

        Xo = np.array([Lx, 0.0])
        No = np.array([1.0, 0.0])

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[0],Lx):
                # locate the node at the centerline
                if math.isclose(X[1],Ly/2.):
                    end_node = node
                # load the end faces
                for elem in node.elements:
                    for face in elem.faces:
                        for x, area in zip(face.pos, face.area):
                            if np.abs( (x - Xo) @ No ) < 1.0e-2 and  No @ area / np.linalg.norm(area):
                                face.setLoad(px, -pxy)

        Xo = np.array([0.0, Ly])
        No = np.array([0.0, 1.0])

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[1],Ly):
                for elem in node.elements:
                    for face in elem.faces:
                        for x, area in zip(face.pos, face.area):
                            if np.abs( (x - Xo) @ No ) < 1.0e-2 and  No @ area / np.linalg.norm(area):
                                face.setLoad(-py, 0.0)

        #model.report()

        # set up a recorder
        model.initRecorder(variables=['ux','uy'], nodes=[end_node])
        model.startRecorder()

        model.plot(factor=0, title="undeformed system", filename="plate11_undeformed.png", show_bc=1, show_loads=1)

        for lf in np.linspace(0.0, target_load_level, max_steps):

            model.setLoadFactor(lf)
            model.solve(verbose=True)

            #model.report()


        model.plot(factor=1., filename=f"plate11_deformed_lf{lf:.2f}.png", show_bc=1, show_loads=1, show_reactions=1)

        model.valuePlot('ux', filename=f"plate11_ux_lf{lf:.2f}.png")
        model.valuePlot('uy', show_mesh=True, filename=f"plate11_uy_lf{lf:.2f}.png")

        # create a history plot for the end node

        model.historyPlot('lam', ['ux','uy'], nodes=[end_node,end_node])
        model.historyPlot(('ux',end_node), 'uy', node=end_node)


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExamplePlate11()
    ex.run()


