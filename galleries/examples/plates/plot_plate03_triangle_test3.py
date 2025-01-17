r"""
===============================================================
A square patch made of two triangular plate elements
===============================================================

Basic implementation test with applied loads.
    Testing the tangent stiffness computation
    for a :func:`Triangle6` (using quadratic shape functions).

Using

* :py:class:`elements.linear.Triangle6` (see :ref:`triangle6_class`)
* :py:class:`materials.PlaneStress`  (see :ref:`plane_stress_material_class`)

.. code::

    free   free
     ^     ^
     |     |
     3--6--2 -> free
     |\  b | >
     | \   | >
     7  8  5 > (w = 1.0)
     |   \ | >
     | a  \| >
     0--4--1 -> free

    width:  10.
    height: 10.

    Material parameters: St. Venant-Kirchhoff, plane stress
        E  = 10.0
        nu =  0.30
        t  =  1.0

    Element loads:
        node 0: [ 0.0, 0.0]
        node 1: [ 10./6, 0.0]
        node 2: [ 10./6, 0.0]
        node 3: [ 0.0, 0.0]
        node 4: [ 0.0, 0.0]
        node 5: [ 10.*2/3, 0.0]
        node 6: [ 0.0, 0.0]
        node 7: [ 0.0, 0.0]
        node 8: [ 0.0, 0.0]

    2nd Piola-Kirchhoff stress:
        S_XX =  w                  =  1.000
        S_YY = S_XY = S_YX = S_ZZ  =  0.000

    Green Lagrange strain:
        eps_XX = (1/E) ((1.000) - (0.30)(0.000)) =  0.100
        eps_YY = (1/E) ((0.000) - (0.30)(1.000)) = -0.030
        eps_XY = eps_YX                          =  0.000
        eps_ZZ = -nu * (eps_XX + eps_YY)         = -0.021

    Stretches:
        lam_X = sqrt(1 + 2 eps_XX) = 1.095
        lam_Y = sqrt(1 + 2 eps_YY) = 0.9695

    Displacements:
        ux = (lam_X - 1) * x, uy = (lam_Y - 1) * y
        node 0: [ 0.000,  0.000  ]
        node 1: [ 0.954,  0.000  ]
        node 2: [ 0.954, -0.305  ]
        node 3: [ 0.000, -0.305  ]
        node 4: [ 0.477,  0.000  ]
        node 5: [ 0.954, -0.1525 ]
        node 6: [ 0.477, -0.305  ]
        node 7: [ 0.954, -0.1525 ]
        node 8: [ 0.477, -0.1525 ]

Author: Peter Mackenzie-Helnwein
"""

from femedu.examples import Example

from femedu.domain import System, Node
from femedu.solver import NewtonRaphsonSolver
from femedu.elements.linear import Triangle6
from femedu.materials import PlaneStress


class ExamplePlate02b(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = r"""
    ## A square patch made of two triangular plate elements

    Basic implementation test with applied loads.
    Testing the tangent stiffness computation. 

    free   free
     ^     ^
     |     |
     3--6--2 -> free
     |\  b | >
     | \   | >
     7  8  5 > (w = 1.0)
     |   \ | >
     | a  \| >
     0--4--1 -> free

    width:  10.
    height: 10.

    Material parameters: St. Venant-Kirchhoff, plane stress
        E  = 10.0
        nu =  0.30
        t  =  1.0
        
    Author: Peter Mackenzie-Helnwein 
    """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):

        params = dict(
            E  = 10.,   # Young's modulus
            nu = 0.3,   # Poisson's ratio
            t  = 1.0,   # thickness of the plate
            fy = 1.e30  # yield stress
        )

        a = 10.     # length of the plate in the x-direction
        b = 10.     # length of the plate in the y-direction

        model = System()
        model.setSolver(NewtonRaphsonSolver())

        nd0 = Node( 0.0, 0.0)
        nd1 = Node(   a, 0.0)
        nd2 = Node(   a,   b)
        nd3 = Node( 0.0,   b)
        nd4 = Node( a/2, 0.0)
        nd5 = Node(   a, b/2)
        nd6 = Node( a/2,   b)
        nd7 = Node( 0.0, b/2)
        nd8 = Node( a/2, b/2)

        nd0.fixDOF('ux', 'uy')
        nd1.fixDOF('uy')
        nd3.fixDOF('ux')

        model.addNode(nd0, nd1, nd2, nd3, nd4, nd5, nd6, nd7, nd8)

        elemA = Triangle6(nd0, nd1, nd3, nd4, nd8, nd7, PlaneStress(params))
        elemB = Triangle6(nd2, nd3, nd1, nd6, nd8, nd5, PlaneStress(params))

        model.addElement(elemA, elemB)

        elemA.setSurfaceLoad(face=2, pn=1.0)
        elemB.setSurfaceLoad(face=2, pn=1.0)

        model.setLoadFactor(1.0)

        nd0.setDisp([0.0,  0.0])
        nd1.setDisp([5.0,  0.0])
        nd2.setDisp([5.0, -5.0])
        nd3.setDisp([0.0, -5.0])
        nd4.setDisp([2.5,  0.0])
        nd5.setDisp([5.0, -2.5])
        nd6.setDisp([2.5, -5.0])
        nd7.setDisp([0.0, -2.5])
        nd8.setDisp([2.5, -2.5])

        elemA.updateState()
        elemB.updateState()

        model.report()

        model.plot(factor=1.0, filename="plate02b_deformed.png")


        model.setLoadFactor(0.0)

        # model.solver.assemble()
        # model.solver.showKt()
        #
        model.solve()

        model.report()  # activate this line for lots of debug info
        model.plot(factor=0.0, title="Undeformed system", filename="plate02b_undeformed.png", show_bc=1)


        model.setLoadFactor(1.0)
        model.solve(verbose=1)
        model.plot(factor=1.0, filename="plate02b_deformed.png")

        model.report()

        # requires femedu-1.0.25 or newer
        model.valuePlot('sxx', show_mesh=True)
        model.valuePlot('syy', show_mesh=True)
        model.valuePlot('sxy', show_mesh=True)


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExamplePlate02b()
    ex.run()


