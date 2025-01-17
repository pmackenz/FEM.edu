"""
==========================================================
Patch test for triangular plate under in-plane loading
==========================================================

The patch test is an empirical minimum test which every finite element has to pass to ensure convergence with mesh refinement.

It consists of a problem for which a known homogeneous solution exists.
For plates, we commonly use a rectangular plate subject to homogeneous edge loading,
e.g., constant tension in the x-direction, or constant shear, etc.

The mesh must contain distorted elements and at least one element not attached to any node on the boundary.

Author: Peter Mackenzie-Helnwein
"""
import numpy as np

from femedu.examples import Example

from femedu.domain import System, Node
from femedu.solver import NewtonRaphsonSolver
from femedu.elements.linear import Triangle
from femedu.materials import PlaneStress


class ExamplePlate03(Example):

    # sphinx_gallery_start_ignore
    # sphinx_gallery_thumbnail_number = 3
    def docString(self):
        s = """
    ## Patch test for triangular plate under in-plane loading

    The patch test is an empirical minimum test which every finite element has to pass to ensure convergence with mesh refinement.

    It consists of a problem for which a known homogeneous solution exists.  For plates, we commonly use a rectangular plate subject to homogenous edge loading, e.g., constant tension in the x-direction, or constant shear, etc.

    The mesh must contain distorted elements and at least one element not attached to any node on the boundary.
        
    Author: Peter Mackenzie-Helnwein 
    """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):
        # ========== setting mesh parameters ==============

        N = 8         # number of elements in the mesh
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

        nodes = (
            Node(0.0*Lx, 0.0*Ly),  # nd 0
            Node(0.2*Lx, 0.0*Ly),  # nd 1
            Node(0.5*Lx, 0.0*Ly),  # nd 2
            Node(0.7*Lx, 0.0*Ly),  # nd 3
            Node(1.0*Lx, 0.0*Ly),  # nd 4
            #
            Node(0.0*Lx, 0.2*Ly),  # nd 5
            Node(0.15*Lx,0.3*Ly),  # nd 6
            Node(0.5*Lx, 0.2*Ly),  # nd 7
            Node(0.8*Lx, 0.3*Ly),  # nd 8
            Node(1.0*Lx, 0.2*Ly),  # nd 9
            #
            Node(0.0*Lx, 0.6*Ly),  # nd 10
            Node(0.2*Lx, 0.5*Ly),  # nd 11
            Node(0.7*Lx, 0.7*Ly),  # nd 12
            Node(0.9*Lx, 0.6*Ly),  # nd 13
            Node(1.0*Lx, 0.8*Ly),  # nd 14
            #
            Node(0.0*Lx, 1.0*Ly),  # nd 15
            Node(0.3*Lx, 1.0*Ly),  # nd 16
            Node(0.55*Lx,1.0*Ly),  # nd 17
            Node(0.8*Lx, 1.0*Ly),  # nd 18
            Node(1.0*Lx, 1.0*Ly),  # nd 19
        )

        elements = (
            Triangle(nodes[0],nodes[1],nodes[5],PlaneStress(params)),  # elem 0
            Triangle(nodes[1],nodes[2],nodes[6],PlaneStress(params)),  # elem 1
            Triangle(nodes[2],nodes[3],nodes[7],PlaneStress(params)),  # elem 2
            Triangle(nodes[3],nodes[4],nodes[8],PlaneStress(params)),  # elem 3
            #
            Triangle(nodes[6],nodes[5],nodes[1],PlaneStress(params)),  # elem 4
            Triangle(nodes[7],nodes[6],nodes[2],PlaneStress(params)),  # elem 5
            Triangle(nodes[8],nodes[7],nodes[3],PlaneStress(params)),  # elem 6
            Triangle(nodes[9],nodes[8],nodes[4],PlaneStress(params)),  # elem 7
            #
            Triangle(nodes[5],nodes[6],nodes[10],PlaneStress(params)),  # elem 8
            Triangle(nodes[6],nodes[7],nodes[11],PlaneStress(params)),  # elem 9
            Triangle(nodes[7],nodes[8],nodes[12],PlaneStress(params)),  # elem 10
            Triangle(nodes[8],nodes[9],nodes[13],PlaneStress(params)),  # elem 11
            #
            Triangle(nodes[11],nodes[10],nodes[6],PlaneStress(params)),  # elem 12
            Triangle(nodes[12],nodes[11],nodes[7],PlaneStress(params)),  # elem 13
            Triangle(nodes[13],nodes[12],nodes[8],PlaneStress(params)),  # elem 14
            Triangle(nodes[14],nodes[13],nodes[9],PlaneStress(params)),  # elem 15
            #
            Triangle(nodes[10],nodes[11],nodes[15],PlaneStress(params)),  # elem 16
            Triangle(nodes[11],nodes[12],nodes[16],PlaneStress(params)),  # elem 17
            Triangle(nodes[12],nodes[13],nodes[17],PlaneStress(params)),  # elem 18
            Triangle(nodes[13],nodes[14],nodes[18],PlaneStress(params)),  # elem 19
            #
            Triangle(nodes[16],nodes[15],nodes[11],PlaneStress(params)),  # elem 20
            Triangle(nodes[17],nodes[16],nodes[12],PlaneStress(params)),  # elem 21
            Triangle(nodes[18],nodes[17],nodes[13],PlaneStress(params)),  # elem 22
            Triangle(nodes[19],nodes[18],nodes[14],PlaneStress(params)),  # elem 23
        )

        model.addNode(*nodes)
        model.addElement(*elements)

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
        elements[ 8].setSurfaceLoad(2,px)
        elements[16].setSurfaceLoad(2,px)

        # surface loads on the right side
        elements[ 7].setSurfaceLoad(2,px)
        elements[15].setSurfaceLoad(2,px)
        elements[23].setSurfaceLoad(2,px)

        # these are only nodal forces as part of the reference load
        # .. load only the upper node
        #print(model)

        model.report()

        model.plotSystem()

        model.setLoadFactor(10.0)
        model.solve()

        model.solver.showKt()

        model.report()

        model.plot(factor=25., filename="plate08_deformed.png")

        # tmpl = "{:12.4f} {:12.4f}"
        # for elem in model.elements:
        #     print(tmpl.format(elem.area, elem.material.getStress()['xx']))

        # requires femedu-1.0.25 or newer
        model.valuePlot('sxx', show_mesh=True)

        # tmpl = "{:>12} {:>5} {:12.4f} {:12.4f} {:12.4f}"
        # for node in model.nodes:
        #     print(tmpl.format(node.getID(), node._mapped_variable, node._weighted_value, node._weight, node.getMappedValue('sxx')))

        model.valuePlot('syy', show_mesh=True)
        model.valuePlot('sxy', show_mesh=True)


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExamplePlate03()
    ex.run()


