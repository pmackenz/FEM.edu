"""
=================================================================
3d truss example demonstrating large deformation analysis.
=================================================================

The system is statically determined and allows for easy validation of
calculated deformation, reactions and internal forces.

Author: Tatsu Sweet
"""

# %%
# Setup

from femedu.examples.Example import *

from femedu.domain.System import *
from femedu.domain.Node import *
from femedu.elements.linear.Truss import *
from femedu.materials.FiberMaterial import *


class ExampleTruss04(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = """
        3d truss example demonstrating large deformation analysis.
        
        Author: Tatsu 
        """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):

        # initialize a system model
        params = {'E': 2100., 'A': 1., 'nu': 0.0, 'fy': 1.e30}

        model = System()

        # create nodes
        H = 5
        nd1 = Node(0.0, 5.0, 0.0)
        nd2 = Node(9.5, 5.0, 0.0)
        nd3 = Node(0.0, 0.0, 0.0)
        nd4 = Node(9.5, 0.0, 0.0)
        nd5 = Node(5.5, 3.75, H)
        nd6 = Node(5.5, 1.25, H)

        nodeList = [nd1, nd2, nd3, nd4, nd5, nd6]
        model.addNode(*nodeList)


        # create elements
        model.addElement(Truss(nd1, nd5, FiberMaterial(params)))  # bottom 1
        model.addElement(Truss(nd1, nd6, FiberMaterial(params)))  # up right diag 1
        model.addElement(Truss(nd2, nd5, FiberMaterial(params)))  # up left diag 1
        model.addElement(Truss(nd3, nd6, FiberMaterial(params)))  # bottom 1
        model.addElement(Truss(nd4, nd5, FiberMaterial(params)))  # up right diag 1
        model.addElement(Truss(nd4, nd6, FiberMaterial(params)))  # up left diag 1
        model.addElement(Truss(nd5, nd6, FiberMaterial(params)))  # bottom 1


        # define support(s)
        translation_dofs = ('ux', 'uy', 'uz')
        for node in [nd1, nd2, nd3, nd4]:
            node.fixDOF(*translation_dofs)

        # add loads
        nd5.setLoad((-100.0,), ('uz',))

        # analyze the model
        model.solve()

        # write out report
        model.report()

        # create plots
        model.plot(factor=1.)

# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleTruss04()
    ex.run()

