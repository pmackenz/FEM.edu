r"""
==========================================================
Heat transfer through a thick cylinder
==========================================================

This problem demonstrates the use of prescribed temperature
on both sides of the wall.

Using

* :py:class:`mesher.PatchMesher` (see :ref:`patch_mesher_class`)
* :py:class:`diffusion.Triangle` (see :ref:`diffusion_triangle_class`)
* :py:class:`materials.Thermal`  (see :ref:`diffusion_material_classes`)

Theory
---------
We shall consider a stationary heat transfer problem on  a thick ring.
The inner surface of the ring, :math:`r_i`, is heated to :math:`200~K`,
the outer surface of the ring, :math:`r_o`, to :math:`300~K`.

The thermal equation for an axi-symmetric problem can be expressed as

.. math::
    \Delta T = \frac{1}{r} \: \frac{\partial }{\partial r}
    \left( r \frac{\partial T}{\partial r} \right) = 0

where :math:`\Delta` is the Laplace operator.

The analytic solution follows as

.. math::
    T(r) = {\displaystyle \frac{T_i \log(r_o/r) + T_o \log(r/r_i)}{\log(r_o/r_i)}}

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


class ExampleThermal04(Example):

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

        Nx = 8  # number of elements through the wall
        Ny = 4  # number of elements parallel to the wall
        Lx = 10.00  # wall thickness in m
        Ri =  5.00
        Ro = Ri + Lx
        alpha = np.radians(45.0)

        # ========== setting material parameters ==============

        params = dict(
            specific_heat = 1000,  # J/kg.K
            density       = 2500,  # kg/m3
            conductivity  =  100,  # W/m.K
            thickness     =  1.0   # m
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
        model.valuePlot('qx', show_mesh=True)
        model.valuePlot('qy', show_mesh=True)

        # creating a path plot

        R_list = []
        T_list = []

        for node in nodes:
            X = node.getPos()
            R = np.linalg.norm(X)
            T = node.getDisp('T')
            R_list.append(R)
            T_list.append(T)

        # the analytic solution for comparison
        r = np.linspace(Ri, Ro, 21)
        T = (200. * np.log(Ro/r) + 300. * np.log(r/Ri)) / np.log(Ro/Ri)

        fig, axs = plt.subplots()
        axs.plot(r,T,'-b',label="analytic solution")
        axs.plot(R_list,T_list,'ro',label="FEM")
        axs.set_title('Nodal Temperature for ALL Nodes')
        axs.set_xlabel("Radial distance")
        axs.set_ylabel('T')
        axs.legend()
        axs.grid(True)
        plt.show()



# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleThermal04()
    ex.run()
