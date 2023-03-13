"""
Example plate02
"""
from ...examples.Example import *

from ...domain.System import *
from ...solver.NewtonRaphsonSolver import *
from ...elements.LinearTriangle import *
from ...materials.PlaneStress import *


class ExamplePlate02(Example):

    def docString(self):
        s = """
    ## A square patch made of two triangular plate elements

    Basic implementation test with applied loads.
    Testing the tangent stiffness computation. 

    free   free
     ^     ^
     |     |
     3-----2 -> free
     |\  b | >
     | \   | >
     |  \  | > (w = 1.0)
     |   \ | >
     | a  \| >
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

        nd0.fixDOF('ux', 'uy')
        nd1.fixDOF('uy')
        nd3.fixDOF('ux')

        model.addNode(nd0, nd1, nd2, nd3)

        elemA = LinearTriangle(nd0, nd1, nd3, PlaneStress(params))
        elemB = LinearTriangle(nd2, nd3, nd1, PlaneStress(params))

        model.addElement(elemA, elemB)

        elemB.setSurfaceLoad(face=2, w=1.0)

        model.setLoadFactor(0.0)
        model.solve()
        #model.report()  # activate this line for lots of debug info
        model.plot(factor=1.0,
                   title="Undeformed system",
                   filename="plate02_undeformed.png")

        model.setLoadFactor(1.0)
        model.solve()
        model.plot(factor=1.0, filename="plate02_deformed.png")

        model.report()

