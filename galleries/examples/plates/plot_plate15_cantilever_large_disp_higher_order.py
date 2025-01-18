r"""
==========================================================
Bending a cantilever beam using Quad9 elements
==========================================================

Using PatchMesher to model the beam

.. dropdown::  Background Theory

    This problem can be approximately validated using Bernoulli-Euler theory for
    small deformations. The given problem shall be modeled using

    .. list-table::
        :widths: 20 30 50
        :header-rows: 1

        * - parameter
          - value
          - description
        * - :math:`E`
          - 20000.
          - modulus of elasticity (in ksi)
        * - :math:`I`
          - 666.667
          - area moment of inertia (in :math:`inches^4`)
        * - :math:`L`
          - 120.
          - length of the cantilever (in inches)
        * - :math:`P`
          - 30.
          - force at :math:`x=L` (in kips)

    The general solution then yields

    .. math::
        v(x) = -\frac{P L^3}{6 EI}\left( \frac{x}{L} \right)^2\left( 3 - \frac{x}{L} \right)

    .. math::
        \theta(x) = \frac{d}{dx} v(x) = -\frac{P L^2}{2 EI}\left( \frac{x}{L} \right)\left( 2 - \frac{x}{L} \right)

    .. math::
        M(x) = EI \frac{d}{dx} \theta(x) = -\frac{P L}{6} \left( 1 - \frac{x}{L} \right)

    .. math::
        V(x) = \frac{d}{dx} M(x) = P

    The horizontal movement follows as (:math:`2^{nd}` order accurate)

    .. math::
        u(x) = \int\limits_{0}^{x} -\frac{1}{2} \theta^2(s) \: ds
             = -\frac{P^2 L^5}{120 (EI)^2}\left( \frac{x}{L} \right)^3\left( 20 - 15\:\frac{x}{L}+3 \left( \frac{x}{L} \right)^2 \right)



    .. list-table:: Reference values for a load factor of :math:`\lambda=1.0`
        :widths: 20 30 50
        :header-rows: 1

        * - variable
          - value
          - description
        * - :math:`u(L)`
          - -0.0083981
          - end displacement (in inches). :math:`u>0` means moving to the right.
        * - :math:`v(L)`
          - -1.296
          - end displacement (in inches). :math:`v>0` means moving up.



    .. list-table:: Reference values for a load factor of :math:`\lambda=10.0`
        :widths: 20 30 50
        :header-rows: 1

        * - variable
          - value
          - description
        * - :math:`u(L)`
          - -0.83981
          - end displacement (in inches). :math:`u>0` means moving to the right.
        * - :math:`v(L)`
          - -12.96
          - end displacement (in inches). :math:`v>0` means moving up.

"""

import numpy as np

from femedu.examples import Example

from femedu.domain import System
from femedu.solver import NewtonRaphsonSolver
from femedu.elements.linear import Quad9
#from femedu.elements.finite import Quad9
from femedu.materials import PlaneStress
from femedu.mesher import *


class ExamplePlate15(Example):

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

        Nx = 4  # number of elements in the mesh
        Ny = 4  # number of elements in the mesh
        Lx = 120.0  # length of plate in the x-direction
        Ly = 20.0  # length of plate in the y-direction

        # ========== setting material parameters ==============

        params = dict(
            E=20000.,  # Young's modulus
            nu=0.250,  # Poisson's ratio
            t=1.00  # thickness of the plate
        )

        # ========== setting load parameters ==============

        px = 0.0  # uniform load normal to x=Lx
        py = 0.0  # uniform load normal to y=Ly
        pxy = 1.5  # uniform shear load on x=L

        # ========== setting analysis parameters ==============

        target_load_level = 10.00  # reference load
        max_steps = 10  # number of load steps: 2 -> [0.0, 1.0]

        # define a list of target load levels
        load_levels = np.linspace(0, target_load_level, max_steps+1)

        #
        # ==== Build the system model ====
        #

        model = System()
        model.setSolver(NewtonRaphsonSolver())

        # create nodes

        mesher = PatchMesher(model, (0., 0.), (Lx, 0.), (Lx, Ly), (0., Ly))
        nodes, elements = mesher.quadMesh(Nx, Ny, Quad9, PlaneStress(params))

        # define support(s)

        ## find nodes at x==0
        for node, _ in model.findNodesAlongLine((0.0, 0.0), (0.0, 1.0)):
            node.fixDOF('ux', 'uy')

        # ==== complete the reference load ====

        # the section at the right end
        for _, face in model.findFacesAlongLine((Lx, 0.0), (0.0, 1.0), orientation=+1):
            face.setLoad(px, -pxy)

        # durface loading on the top face
        for _, face in model.findFacesAlongLine((0.0, Ly), (1.0, 0.0), orientation=-1):
            face.setLoad(-py, 0.0)

        # find the node on the beam axis (y==Ly/2) at the end of the beam (x==Lx)
        end_node, _ = model.findNodesAt((Lx, Ly / 2))[0]

        # set up a recorder
        model.initRecorder(variables=['ux', 'uy'], nodes=[end_node])
        model.startRecorder()

        model.plot(factor=0, title="undeformed system", filename="plate11_undeformed.png", show_bc=1, show_loads=1)

        for lf in load_levels:
            model.setLoadFactor(lf)

            for node, _ in model.findNodesAlongLine((Lx, 0.), (0., 1.)):
                print(node)

            for _, face in model.findFacesAlongLine((Lx, 0.), (0., 1.)):
                print(face)


            model.solve(verbose=True)

            # model.report()

        model.plot(factor=1., filename=f"plate11_deformed_lf{lf:.2f}.png", show_bc=1, show_loads=1, show_reactions=1)

        model.valuePlot('sxx', show_mesh=True)
        model.valuePlot('sxy', show_mesh=True)

        # create a history plot for the end node

        model.historyPlot('lam', ['ux', 'uy'], nodes=[end_node, end_node])
        model.historyPlot(('ux', end_node), 'uy', node=end_node)


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExamplePlate15()
    ex.run()


