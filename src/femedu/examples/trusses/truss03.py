"""
Example
"""
from ...examples.Example import *

from ...domain.System import *
from ...domain.Node import *
from ...elements.Truss import *
from ...materials.FiberMaterial import *


class ExampleTruss03(Example):

    def docString(self):
        s = """
        Simple triangular truss. 
        
        This example is structurally identical to ExampleTruss01 but utilizes 
        the alternative input style for adding Nodes and Elements to the model.
        
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

        model += nd0
        model += nd1
        model += nd2

        # create elements
        model += Truss(nd0, nd1, FiberMaterial(params))  # bottom 1
        model += Truss(nd0, nd2, FiberMaterial(params))  # up right diag 1
        model += Truss(nd1, nd2, FiberMaterial(params))  # up left diag 1

        # define support(s)
        nd0.fixDOF('ux')    # horizontal support left end
        #nd0 //= 0
        nd0.fixDOF('uy')    # vertical support left end
        nd1.fixDOF('uy')    # vertical support right end

        # add loads
        # .. load only the upper nodes
        nd2.setLoad((0.0, -1.0), ('ux','uy'))

        # analyze the model
        model.solve()

        # write out report
        model.report()

        # create plots
        model.plot(factor=1., filename="truss03_deformed_a.png")

        # fix horizontal motion of node 1
        nd1.fixDOF('ux')

        # add loads: same load -- nothing to do

        # RE-analyze the model
        model.resetDisp()
        model.solve()

        # skip the report
        model.report()

        # create plots
        model.plot(factor=1., filename="truss03_deformed_b.png")

