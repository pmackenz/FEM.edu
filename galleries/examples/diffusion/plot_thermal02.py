"""
==========================================================
Heat transfer at the corner of a building
==========================================================

Using

* :py:class:`mesher.PatchMesher` (see :ref:`patch_mesher_class`)
* :py:class:`diffusion.Triangle` (see :ref:`diffusion_triangle_class`)
* :py:class:`materials.Thermal`  (see :ref:`diffusion_material_classes`)

"""
import math
import sys
import numpy as np

from femedu.examples.Example import *

from femedu.domain import *
from femedu.mesher import PatchMesher
from femedu.elements.diffusion import Triangle
from femedu.materials.Thermal import *


class ExampleThermal02(Example):

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

        Nx = 8  # number of elements through the wall
        Ny = 8  # number of elements parallel to the wall
        L = 0.30  # wall thickness in m
        alpha = 0.50  # refinement parameter for placing nodes 8-11

        # ========== setting material parameters ==============

        params = dict(
            E=20000.,  # Young's modulus
            nu=0.250,  # Poisson's ratio
            t=1.00  # thickness of the plate
        )

        # ========== setting load parameters ==============

        qn = 10.0  # uniform flux normal to x=const

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

        #  5 -------- 6 -------- 7
        #  |          |          |
        #  |         11          |
        #  |          |          |
        #  2----9---- 3 ---10--- 4
        #             |          |
        #             8          |
        #             |          |
        #             0 -------- 1

        pts = (
            (0, -L),  # 0
            (L, -L),  # 1
            (-L, 0),  # 2
            (0, 0),  # 3
            (L, 0),  # 4
            (-L, L),  # 5
            (0, L),  # 6
            (L, L),  # 7
            (0, -alpha * L),  # 8
            (-alpha * L, 0),  # 9
            (alpha * L, 0),  # 10
            (0, alpha * L),  # 11
        )

        mesher1 = PatchMesher(model, pts[2], pts[3], pts[6], pts[5], pts[9], pts[11], None, None)
        nodes1, elements1 = mesher1.triangleMesh(Ny, Nx, Triangle, Thermal(params))

        mesher2 = PatchMesher(model, pts[3], pts[4], pts[7], pts[6], pts[10], None, None, pts[11])  # center node
        nodes2, elements2 = mesher2.triangleMesh(Nx, Nx, Triangle, Thermal(params))

        mesher3 = PatchMesher(model, pts[0], pts[1], pts[4], pts[3], None, None, pts[10], pts[8])
        nodes3, elements3 = mesher3.triangleMesh(Nx, Ny, Triangle, Thermal(params))

        mesher1.tie(mesher2)
        mesher2.tie(mesher3)

        nodes = nodes1 + nodes2 + nodes3
        elements = elements1 + elements2 + elements3

        model.plot(factor=0.0,
                   title='Meshing the corner in a wall',
                   show_reactions=0, show_bc=0, show_loads=0)

        model.report()

        # define support(s)

        ## find nodes at y==0 and x==0

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[0], 0.0) and X[1] <= 0.0:
                node.fixDOF('T')  # horizontal support left side
            if math.isclose(X[1], 0.0) and X[0] <= 0.0:
                node.fixDOF('T')  # vertical support at y==0

        # ==== complete the reference load ====

        loaded_elements = []

        Xo = np.array([L, 0.0])
        No = np.array([1.0, 0.0])

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[0], L):
                print(node)
                for elem in node.elements:
                    print('+', elem)
                    for face in elem.faces:
                        for x, area in zip(face.pos, face.area):
                            if np.abs((x - Xo) @ No) < 1.0e-4 and No @ area / np.linalg.norm(area) > 1.0e-2:
                                face.setFlux(qn)
                                loaded_elements.append(elem)

        Xo = np.array([0.0, L])
        No = np.array([0.0, 1.0])

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[1], L):
                print(node)
                for elem in node.elements:
                    print('-', elem)
                    for face in elem.faces:
                        for x, area in zip(face.pos, face.area):
                            if np.abs((x - Xo) @ No) < 1.0e-4 and No @ area / np.linalg.norm(area) > 1.0e-2:
                                face.setFlux(qn)
                                loaded_elements.append(elem)

        # perform the analysis
        model.setLoadFactor(1.0)
        model.solve()

        # model.report()
        for elem in loaded_elements:
            print(elem)

        model.valuePlot('T', show_mesh=True)


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleThermal02()
    ex.run()
