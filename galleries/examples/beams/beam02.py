"""
Example
"""
from ...examples.Example import *

from ...domain.System import *
from ...domain.Node import *
from ...elements.Beam2D import *
from ...materials.ElasticSection import *


class ExampleBeam02(Example):

    def docString(self):
        s = """
        Three-span continuous beam under uniform load. 
        
        The system is statically determined and allows for easy validation of 
        calculated deformation, reactions and internal forces.
        
        Author: Peter Mackenzie-Helnwein 
        """
        return s

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
        model.plot(factor=100., filename="beam02_deformed.png")

        model.beamValuePlot('V', filename="beam02_shear.png")
        model.beamValuePlot('M', filename="beam02_moment.png")

