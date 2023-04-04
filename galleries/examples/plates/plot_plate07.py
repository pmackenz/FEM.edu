"""
===============================================================
A square patch made of one quadrilateral plate elements
===============================================================

Basic implementation test with applied loads.
    Testing the tangent stiffness computation.

.. code::

    free   free
     ^     ^
     |     |
     3-----2 -> free
     |     | >
     |  a  | > (w = 1.0)
     |     | >
     0-----1 -> free

    width:  10.
    height: 10.

    Material parameters: St. Venant-Kirchhoff, plane stress
        E  = 10.0
        nu =  0.30
        t  =  1.0

    Element loads:
        node 0: [ 0.0, 0.0]
        node 1: [ 5.0, 0.0]
        node 2: [ 5.0, 0.0]
        node 3: [ 0.0, 0.0]

Author: Peter Mackenzie-Helnwein
"""
import numpy as np

from femedu.examples.Example import *

from femedu.domain.System import *
from femedu.solver.NewtonRaphsonSolver import *
from femedu.elements.linear.Quad import *
from femedu.materials.PlaneStress import *


class ExamplePlate07(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = """
    ## A square patch made of one quadrilateral plate elements

    Basic implementation test with applied loads.
    Testing the tangent stiffness computation. 

    free   free
     ^     ^
     |     |
     3-----2 -> free
     |     | >
     |  a  | > (w = 1.0)
     |     | >
     0-----1 -> free
    
    width:  10.
    height: 10.
    
    Material parameters: St. Venant-Kirchhoff, plane stress
        E  = 10.0
        nu =  0.30
        t  =  1.0
        
    Element loads:
        node 0: [ 0.0, 0.0]
        node 1: [ 5.0, 0.0]
        node 2: [ 5.0, 0.0]
        node 3: [ 0.0, 0.0]
    
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
        node 0: [ 0.000,  0.000 ]
        node 1: [ 0.954,  0.000 ]
        node 2: [ 0.954, -0.305 ]
        node 3: [ 0.000, -0.305 ]
        
    Author: Peter Mackenzie-Helnwein 
    """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):

        params = dict(
            E  = 10., # Young's modulus
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

        # nd0.fixDOF('ux', 'uy')
        # nd1.fixDOF('uy')
        # nd3.fixDOF('ux')

        model.addNode(nd0, nd1, nd2, nd3)

        elemA = Quad(nd0, nd1, nd2, nd3, PlaneStress(params))

        model.addElement(elemA)

        elemA.setSurfaceLoad(face=1, pn=1.0)

        model.plot(factor=0.0, title="Undeformed system", filename="plate07_undeformed.png", show_bc=1)

        # %%
        # We can have a quick look at the stiffness mode shapes using the
        # buckling-mode plotter.  These are simply eigenvalues and eigenvectors of Kt
        # at the current load level (0.0)
        #
        model.setLoadFactor(0.0)
        model.solve()

        np.save('../../../Kplate.npy', model.solver.Kt)

        for k in range(8):
            name = f"plate07_mode{k:2d}.png"
            model.plotBucklingMode(mode=k,filename=name,factor=1.0)

        # %%
        # Note the three rigid body modes (lam=0.0). It can be shown that all three
        # are limear combinations of translations in x and y-directions and a
        # rigid body rotation.
        #

        # %%
        # Now it is time to add boundary conditions, apply loads
        # and check the convergence behavior.
        #

        nd0.fixDOF('ux', 'uy')
        nd1.fixDOF('uy')
        nd3.fixDOF('ux')

        model.setLoadFactor(1.0)
        model.solve()

        # %%
        # The output shows that we do have a quadratic rate of convergence.

        # %%
        # Let's finish off with a nice plot of the deformed system.

        model.plot(factor=1.0, filename="plate07_deformed.png")

        model.report()


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExamplePlate07()
    ex.run()


