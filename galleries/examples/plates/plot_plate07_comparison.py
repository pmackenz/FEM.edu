r"""
======================================================
A square patch: element comparison
======================================================
Basic implementation test with all prescribed displacements.
    Testing the internal force computation.

.. code::

    v=-5   v=-5
     |     |
     v     v
     3-----2 -> u=5
     |\  b | >
     | \   | >
     |  \  | > (w = 1.0)
     |   \ | >
     | a  \| >
     0-----1 -> u=5

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

    Green Lagrange strain:
        eps_XX = 0.5 * ((1.5)^2 - 1)      =  0.625
        eps_YY = 0.5 * ((0.5)^2 - 1)      = -0.375
        eps_XY = eps_YX                   =  0.000
        eps_ZZ = - nu * (eps_XX + eps_YY) = -0.075

    2nd Piola-Kirchhoff stress:
        D = E t/(1 - nu^2) = 10.989
        S_XX = (10.989) * ((0.625) + (0.30)(-0.375)) =  5.632
        S_YY = (10.989) * ((-0.375) + (0.30)(0.625)) = -2.060
        S_XY = S_YX = S_ZZ                           =  0.000


Author: Peter Mackenzie-Helnwein

"""

from femedu.examples import Example

from femedu.domain import System, Node
from femedu.elements.linear import Triangle
from femedu.elements.linear import Triangle6
from femedu.elements.linear import Quad
from femedu.elements.linear import Quad9
from femedu.materials import PlaneStress


class ExamplePlate07(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = r"""
    ## A square patch made of two triangular plate elements

    Basic implementation test with all prescribed displacements.  
    Testing the internal force computation. 

    v=-5   v=-5
     |     |
     v     v
     3-----2 -> u=5
     |\  b | >
     | \   | >
     |  \  | > (w = 1.0)
     |   \ | >
     | a  \| >
     0-----1 -> u=5

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

    Green Lagrange strain:
        eps_XX = 0.5 * ((1.5)^2 - 1)      =  0.625
        eps_YY = 0.5 * ((0.5)^2 - 1)      = -0.375
        eps_XY = eps_YX                   =  0.000
        eps_ZZ = - nu * (eps_XX + eps_YY) = -0.075

    2nd Piola-Kirchhoff stress:
        D = E t/(1 - nu^2) = 10.989
        S_XX = (10.989) * ((0.625) + (0.30)(-0.375)) =  5.632
        S_YY = (10.989) * ((-0.375) + (0.30)(0.625)) = -2.060
        S_XY = S_YX = S_ZZ                           =  0.000

    Author: Peter Mackenzie-Helnwein 
    """
        return s

    # sphinx_gallery_end_ignore
    # sphinx_gallery_thumbnail_number = 2
    def problem(self):
        params = dict(
            E=10.,  # Young's modulus
            nu=0.3,  # Poisson's ratio
            t=1.0,  # thickness of the plate
            fy=1.e30  # yield stress
        )

        a = 10.  # length of the plate in the x-direction
        b = 10.  # length of the plate in the y-direction

        model = System()

        # made of Triangle elements
        """
        Nodes:
        03   02
        00   01
        """
        nd00 = Node(0.0, 0.0)
        nd01 = Node(a, 0.0)
        nd02 = Node(a, b)
        nd03 = Node(0.0, b)

        model.addNode(nd00, nd01, nd02, nd03)

        elemA0 = Triangle(nd00, nd01, nd03, PlaneStress(params))
        elemA1 = Triangle(nd02, nd03, nd01, PlaneStress(params))

        model.addElement(elemA0, elemA1)

        elemA1.setSurfaceLoad(face=2, pn=1.0)


        # made of Triangle6 elements
        xo=22.
        yo= 0.

        """
        Nodes:
        16   17   18
        13   14   15
        10   11   12
        """
        nd10 = Node(xo+0.0, yo+0.0)
        nd11 = Node(xo+a/2, yo+0.0)
        nd12 = Node(xo+a,   yo+0.0)
        nd13 = Node(xo+0.0, yo+b/2)
        nd14 = Node(xo+a/2, yo+b/2)
        nd15 = Node(xo+a,   yo+b/2)
        nd16 = Node(xo+0.0, yo+b  )
        nd17 = Node(xo+a/2, yo+b  )
        nd18 = Node(xo+a,   yo+b  )

        model.addNode(nd10,nd11,nd12,nd13,nd14,nd15,nd16,nd17,nd18)

        elemB0 = Triangle6(nd10, nd12, nd16, nd11, nd14, nd13, PlaneStress(params))
        elemB1 = Triangle6(nd16, nd12, nd18, nd14, nd15, nd17, PlaneStress(params))

        model.addElement(elemB0, elemB1)

        elemB1.setSurfaceLoad(face=1, pn=1.0)


        # made of Quad9 elements
        xo=22.
        yo=15.

        """
        Nodes:
        26   27   28
        23   24   25
        20   21   22
        """
        nd20 = Node(xo+0.0, yo+0.0)
        nd21 = Node(xo+a/2, yo+0.0)
        nd22 = Node(xo+a,   yo+0.0)
        nd23 = Node(xo+0.0, yo+b/2)
        nd24 = Node(xo+a/2, yo+b/2)
        nd25 = Node(xo+a,   yo+b/2)
        nd26 = Node(xo+0.0, yo+b  )
        nd27 = Node(xo+a/2, yo+b  )
        nd28 = Node(xo+a,   yo+b  )

        model.addNode(nd20,nd21,nd22,nd23,nd24,nd25,nd26,nd27,nd28)

        elemC0 = Quad9(nd20, nd22, nd28, nd26, nd21, nd25, nd27, nd23, nd24, PlaneStress(params))

        model.addElement(elemC0)

        elemC0.setSurfaceLoad(face=1, pn=1.0)


        # made of Quad elements
        xo=0.
        yo=15.

        """
        Nodes:
        32   33
        30   31
        """
        nd30 = Node(xo+0.0, yo+0.0)
        nd31 = Node(xo+a,   yo+0.0)
        nd32 = Node(xo+0.0, yo+b  )
        nd33 = Node(xo+a,   yo+b  )

        model.addNode(nd30,nd31,nd32,nd33)

        elemD0 = Quad(nd30, nd31, nd33, nd32, PlaneStress(params))

        model.addElement(elemD0)

        elemD0.setSurfaceLoad(face=1, pn=1.0)

        # show the model
        model.plot(factor=0, title="Undeformed system", filename="plate07_undeformed.png", show_bc=1)

        model.report()

        # analyze the model
        model.setLoadFactor(1.0)

        nd00.setDisp([0.0, 0.0])
        nd01.setDisp([5.0, 0.0])
        nd02.setDisp([5.0, -5.0])
        nd03.setDisp([0.0, -5.0])

        nd10.setDisp([0.0,  0.0])
        nd13.setDisp([0.0, -2.5])
        nd16.setDisp([0.0, -5.0])
        #
        nd11.setDisp([2.5,  0.0])
        nd14.setDisp([2.5, -2.5])
        nd17.setDisp([2.5, -5.0])
        #
        nd12.setDisp([5.0,  0.0])
        nd15.setDisp([5.0, -2.5])
        nd18.setDisp([5.0, -5.0])

        nd20.setDisp([0.0,  0.0])
        nd23.setDisp([0.0, -2.5])
        nd26.setDisp([0.0, -5.0])
        #
        nd21.setDisp([2.5,  0.0])
        nd24.setDisp([2.5, -2.5])
        nd27.setDisp([2.5, -5.0])
        #
        nd22.setDisp([5.0,  0.0])
        nd25.setDisp([5.0, -2.5])
        nd28.setDisp([5.0, -5.0])

        nd30.setDisp([0.0, 0.0])
        nd31.setDisp([5.0, 0.0])
        nd32.setDisp([0.0, -5.0])
        nd33.setDisp([5.0, -5.0])

        for elem in model.elements:
            elem.updateState()

        model.report()

        model.plot(factor=1.0)
        model.valuePlot('sxx', show_mesh=1)
        model.valuePlot('syy', show_mesh=1)
        model.valuePlot('sxy', show_mesh=1)


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExamplePlate07()
    ex.run()
