"""
==========================================================
Heat transfer through a thick cylinder
==========================================================

This problem demonstrates the use of prescribed temperature
on both sides of the wall.

Using

* :py:class:`mesher.PatchMesher` (see :ref:`patch_mesher_class`)
* :py:class:`diffusion.Triangle` (see :ref:`diffusion_triangle_class`)
* :py:class:`materials.Thermal`  (see :ref:`diffusion_material_classes`)

"""
import matplotlib.pyplot as plt

import math
import sys
import numpy as np

from femedu.examples.Example import *

from femedu.domain import *
from femedu.mesher import PatchMesher
from femedu.elements.diffusion import Triangle
from femedu.materials import Thermal


class ExampleThermal01(Example):

    # sphinx_gallery_start_ignore
    # sphinx_gallery_thumbnail_number = 2

    def docString(self):
        s = """ Heat transfer through a thick cylinder

    Using

    * mesher.PatchMesher
    * diffusion.Triangle
    * materials.Thermal

    Author: Peter Mackenzie-Helnwein 
    """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):
        # ========== setting mesh parameters ==============

        Nx = 16  # number of elements through the wall
        Ny = 8  # number of elements parallel to the wall
        Lx = 10.00  # wall thickness in m
        Ly =  5.00  # wall thickness in m
        Ri =  5.00
        Ro = Ri + Lx
        alpha = np.radians(45.0)

        # ========== setting material parameters ==============

        params = dict(
            E=20000.,  # Young's modulus
            nu=0.250,  # Poisson's ratio
            t=1.00     # thickness of the plate
        )

        # ========== setting load parameters ==============

        qn = 1.00  # uniform flux normal to x=const

        # ========== setting analysis parameters ==============

        target_load_level = 1.00  # reference load
        max_steps = 2  # number of load steps: 2 -> [0.0, 1.0]

        # define a list of target load levels
        load_levels = np.linspace(0, target_load_level, max_steps)

        #
        # ==== Build the system model ====
        #

        model = System()

        # create nodes

        #  2 -------- 3
        #  |          |
        #  |          |
        #  |          |          |
        #  0 -------- 1

        pts = (
            ( Ri,  0),  # 0
            ( Ro,  0),  # 1
            ( Ri*np.cos(alpha), Ri*np.sin(alpha)),  # 2
            ( Ro*np.cos(alpha), Ro*np.sin(alpha)),  # 3
            ( Ri*np.cos(alpha/2), Ri*np.sin(alpha/2)),  # 4
            ( Ro*np.cos(alpha/2), Ro*np.sin(alpha/2)),  # 5
        )

        mesher = PatchMesher(model, pts[0], pts[1], pts[3], pts[2], None, pts[5], None, pts[4])
        nodes, elements = mesher.triangleMesh(Nx, Ny, Triangle, Thermal(params))

        model.plot(factor=0.0,
                   title='Radial diffusion',
                   show_reactions=0, show_bc=0, show_loads=0)

        model.report()

        # boundary condition(s)

        ## find nodes at y==0 and x==0

        for node in nodes:
            X = node.getPos()
            R = np.linalg.norm(X)
            if math.isclose(R, Ri, rel_tol=0.02):
                node.setDOF(['T'],[200.])    # prescribed temperature at x=0.0
            if math.isclose(R, Ro, rel_tol=0.02):
                node.setDOF(['T'],[300.])    # prescribed temperature at x=0.0

        # perform the analysis
        model.setLoadFactor(1.0)
        model.solve()

        model.report()

        model.valuePlot('T', show_mesh=True)

        # creating a path plot

        R_list = []
        T_list = []

        for node in nodes:
            X = node.getPos()
            R = np.linalg.norm(X)
            T = node.getDisp('T')
            R_list.append(R)
            T_list.append(T)

        fig, axs = plt.subplots()
        axs.plot(R_list,T_list,'ro')
        axs.set_title('Nodal Temperature for ALL Nodes')
        axs.set_xlabel("Radial distance")
        axs.set_ylabel('T')
        axs.grid(True)
        plt.show()



# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleThermal01()
    ex.run()
