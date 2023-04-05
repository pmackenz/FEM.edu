"""
======================================================
Three-span continuous beam under uniform load.
======================================================

The system is statically indeterminate but simple enough to validate
deformation, reactions and internal forces.

Author: Peter Mackenzie-Helnwein
"""
from femedu.examples.Example import *

from femedu.domain.System import *
from femedu.domain.Node import *
from femedu.elements.linear.Beam2D import *
from femedu.materials.ElasticSection import *

class ExampleBeam02(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = """
        Three-span continuous beam under uniform load. 
        
        The system is statically indeterminate but simple enough to validate 
        deformation, reactions and internal forces.
        
        Author: Peter Mackenzie-Helnwein 
        """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):
        # initialize a system model
        SpanLengths = [ 8.0 * 12, 10.0 * 12, 8.0 * 12 ]
        Nelems = 2     # number of elements
        params = {'E': 29000., 'A': 5, 'I':50}

        # define load
        w = -1.00

        model = System()

        # meshing parameters
        Xnode  = 0.0
        Ynode  = 0.0
        Offset = 0.0

        # create left node
        nd0 = Node(Xnode, Ynode)
        nd0.fixDOF('ux', 'uy')     # pin support left end
        model += nd0

        # initialization for node and element creation
        ndi = nd0

        for SpanLength in SpanLengths:

            Le = SpanLength / Nelems

            for e in range(Nelems):
                # create next node
                Xnode += Le
                ndj = Node(Xnode, Ynode)
                model += ndj

                # create elements
                elem = Beam2D(ndi, ndj, ElasticSection(params))
                model += elem

                # load the element with a uniform load
                elem.setDistLoad(w)

                # shift one node to the right
                ndi = ndj

            # define support(s)
            ndj.fixDOF('uy')           # roller support right end

            # move on to the next span
            Offset = Xnode

        # done building the model

        # analyze the model
        model.solve()

        # write out report
        model.report()

        # create plots
        model.plot(factor=100., filename="beam02_deformed.png", show_bc=1, show_reactions=1)

        model.beamValuePlot('V', filename="beam02_shear.png")
        model.beamValuePlot('M', filename="beam02_moment.png")


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleBeam02()
    ex.run()
