"""
=======================================================
Statically determinate truss bridge.
=======================================================

The system is statically determined and allows for relatively easy validation of
calculated deformation, reactions and internal forces.

Author: Peter Mackenzie-Helnwein
"""
from femedu.examples.Example import *

from femedu.domain.System import *
from femedu.solver.NewtonRaphsonSolver import NewtonRaphsonSolver
from femedu.domain.Node import *
from femedu.elements.linear.Truss import *
from femedu.materials.FiberMaterial import *


class ExampleTruss02(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = """
        Statically determinate truss bridge.
        
        The system is statically determined and allows for relatively easy validation of 
        calculated deformation, reactions and internal forces.
        
        Author: Peter Mackenzie-Helnwein 
        """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):
        # initialize a system model
        P = -10.0      # reference load on top nodes
        B = 6.0 * 12   # with of one bay in inches
        H = 8.0 * 12   # height of one bay in inches

        # material model parameters
        params = {'E': 10000., 'A': 3., 'nu': 0.0, 'fy': 1.e30}

        model = System()
        model.setSolver(NewtonRaphsonSolver())

        # create nodes
        nd0 = Node(0.0, 0.0)
        nd1 = Node(  B, 0.0)
        nd2 = Node(2*B, 0.0)
        nd3 = Node(3*B, 0.0)
        nd4 = Node(4*B, 0.0)
        nd5 = Node(0.5*B, H)
        nd6 = Node(1.5*B, H)
        nd7 = Node(2.5*B, H)
        nd8 = Node(3.5*B, H)

        model.addNode(nd0, nd1, nd2, nd3, nd4, nd5, nd6, nd7, nd8)

        # create elements
        model.addElement(Truss(nd0, nd1, FiberMaterial(params)))  # bottom 1
        model.addElement(Truss(nd1, nd2, FiberMaterial(params)))  # bottom 2
        model.addElement(Truss(nd2, nd3, FiberMaterial(params)))  # bottom 3
        model.addElement(Truss(nd3, nd4, FiberMaterial(params)))  # bottom 4

        model.addElement(Truss(nd5, nd6, FiberMaterial(params)))  # upper 1
        model.addElement(Truss(nd6, nd7, FiberMaterial(params)))  # upper 2
        model.addElement(Truss(nd7, nd8, FiberMaterial(params)))  # upper 3

        model.addElement(Truss(nd0, nd5, FiberMaterial(params)))  # up right diag 1
        model.addElement(Truss(nd1, nd6, FiberMaterial(params)))  # up right diag 2
        model.addElement(Truss(nd2, nd7, FiberMaterial(params)))  # up right diag 3
        model.addElement(Truss(nd3, nd8, FiberMaterial(params)))  # up right diag 4

        model.addElement(Truss(nd1, nd5, FiberMaterial(params)))  # up left diag 1
        model.addElement(Truss(nd2, nd6, FiberMaterial(params)))  # up left diag 2
        model.addElement(Truss(nd3, nd7, FiberMaterial(params)))  # up left diag 3
        model.addElement(Truss(nd4, nd8, FiberMaterial(params)))  # up left diag 4

        # define support(s)
        nd0.fixDOF('ux', 'uy')    # horizontal support left end
        nd4.fixDOF('uy')            # vertical support right end

        # add loads
        # .. load only the upper nodes
        nd5.setLoad((P,), ('uy',))
        nd6.setLoad((P,), ('uy',))
        nd7.setLoad((P,), ('uy',))
        nd8.setLoad((P,), ('uy',))

        model.setLoadFactor(0.0)
        model.plot(factor=1., filename="truss02_undeformed.png", title="Undeformed System", show_bc=1)

        # analyze the model
        model.setLoadFactor(1.0)
        model.solve()

        # write out report
        model.report()

        # create plots
        model.plot(factor=50.,  filename="truss02_deformed.png", show_loads=1, show_reactions=1)
        model.beamValuePlot('f',filename="truss02_forces.png")


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleTruss02()
    ex.run()

