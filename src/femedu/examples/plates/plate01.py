"""
Example
"""
from ...examples.Example import *

from ...domain.System import *
from ...domain.Node import *
from ...elements.LinearTriangle import *
from ...materials.PlaneStress import *


class ExamplePlate01(Example):

    def problem(self):

        params = {'E': 10., 'A': 1., 'nu': 0.0, 'fy': 1.e30}

        nd0 = Node(  0.0,  0.0)
        nd1 = Node( 10.0, 10.0)
        nd2 = Node(  0.0, 20.0)

        mat = PlaneStress(params)

        elem1 = LinearTriangle(nd0, nd1, nd2, mat)

        nd0.setDisp( 0.0, 0.0)
        nd1.setDisp( 5.0, 5.0)
        nd2.setDisp( 0.0,-5.0)

        elem1.updateState()

        print(elem1)

