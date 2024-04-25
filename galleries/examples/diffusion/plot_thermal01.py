"""
==========================================================
Heat transfer through a wall
==========================================================

This problem demonstrates a combination of prescribed temperature
on the left and prescribed flux on the right side.

Using

* :py:class:`mesher.PatchMesher` (see :ref:`patch_mesher_class`)
* :py:class:`diffusion.Triangle` (see :ref:`diffusion_triangle_class`)
* :py:class:`materials.Thermal`  (see :ref:`diffusion_material_classes`)

Theory
---------
We shall consider a stationary heat transfer problem within a wall.
The inner surface of the wall, :math:`x=0~m`, is heated to :math:`200~K`,
the outer surface of the wall, :math:`x=10~m`, to :math:`300~K`.

The thermal equation for the uni-directional problem can be expressed as

.. math::
    \\Delta T = \\frac{\\partial^2 T}{\\partial r^2}  = 0

where :math:`\Delta` is the Laplace operator.

The analytic solution follows as

.. math::
    T(x) =  T_i  + q_n \: x

This solution will be compared against the finite element solution in the last figure.

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
        s = """ Heat transfer at the corner of a building

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

        Nx = 5  # number of elements through the wall
        Ny = 1  # number of elements parallel to the wall
        Lx = 10.00  # wall thickness in m
        Ly =  1.00  # wall thickness in m

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
            ( 0,  0),  # 0
            (Lx,  0),  # 1
            ( 0, Ly),  # 2
            (Lx, Ly),  # 3
        )

        mesher = PatchMesher(model, pts[0], pts[1], pts[3], pts[2])
        nodes, elements = mesher.triangleMesh(Nx, Ny, Triangle, Thermal(params))

        model.plot(factor=0.0,
                   title='Uni-directional diffusion',
                   show_reactions=0, show_bc=0, show_loads=0)

        model.report()

        # boundary condition(s)

        ## find nodes at y==0 and x==0

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[0], 0.0):
                node.fixDOF('T')    # prescribed temperature at x=0.0
            # if math.isclose(X[0], Lx):
            #     node.fixDOF('T')  # prescribed temperature at x=Lx

        # ==== complete the reference load ====

        Xo = np.array([Lx, 0.0])
        No = np.array([1.0, 0.0])

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[0], Lx):
                print(node)
                for elem in node.elements:
                    print('+', elem)
                    for face in elem.faces:
                        for x, area in zip(face.pos, face.area):
                            if np.abs((x - Xo) @ No) < 1.0e-2 and No @ area / np.linalg.norm(area) > 1.0e-2:
                                face.setFlux(qn)

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
            T = node.getDisp('T')
            R_list.append(X[0])
            T_list.append(T)

        # the analytic solution for comparison
        x = np.linspace(0, Lx, 21)
        T = 0.0 * (1 - x/Lx) + qn * x

        fig, axs = plt.subplots()
        axs.plot(x,T,'-b',label="analytic solution")
        axs.plot(R_list,T_list,'ro',label="FEM")
        axs.set_title('Nodal Temperature for ALL Nodes')
        axs.set_xlabel("X distance")
        axs.set_ylabel('T')
        axs.legend()
        axs.grid(True)
        plt.show()


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleThermal01()
    ex.run()
