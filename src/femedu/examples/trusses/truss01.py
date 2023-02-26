"""
Example
"""
from ...examples.Example import *

from ...domain.System import *
from ...domain.Node import *
from ...elements.Truss import *
from ...materials.FiberMaterial import *


class ExampleTruss01(Example):

    def docString(self):
        s = """
        Simple triangular truss. 
        
        The system is statically determined and allows for easy validation of 
        calculated deformation, reactions and internal forces.
        
        Author: Peter Mackenzie-Helnwein 
        """
        return s

    def problem(self):
        # initialize a system model
        B = 6.0 * 12
        H = 3.0 * 12
        params = {'E': 10., 'A': 1., 'nu': 0.0, 'fy': 1.e30}

        model = System()

        # create nodes
        nd0 = Node(0.0, 0.0)
        nd1 = Node(  B, 0.0)
        nd2 = Node(0.5*B, H)

        model.addNode(nd0, nd1, nd2)

        # create elements
        model.addElement(Truss(nd0, nd1, FiberMaterial(params)))  # bottom 1
        model.addElement(Truss(nd0, nd2, FiberMaterial(params)))  # up right diag 1
        model.addElement(Truss(nd1, nd2, FiberMaterial(params)))  # up left diag 1

        # define support(s)
        nd0.fixDOF('ux', 'uy')    # pin support left end
        nd1.fixDOF('uy')            # roller support right end

        # add loads
        # .. load only the upper nodes
        nd2.setLoad([0.0, -1.0], ('ux', 'uy'))

        # analyze the model
        model.solve()

        # write out report
        model.report()

        # create plots
        model.plot(factor=1.)

        model.beamValuePlot('f')

