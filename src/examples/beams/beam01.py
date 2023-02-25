"""
Example
"""
from examples.Example import *

from domain.System import *
from domain.Node import *
from elements.Beam2D import *
from materials.ElasticSection import *


class ExampleBeam01(Example):

    def docString(self):
        s = """
        Single span beam under uniform load. 
        
        The system is statically determined and allows for easy validation of 
        calculated deformation, reactions and internal forces.
        
        Author: Peter Mackenzie-Helnwein 
        """
        return s

    def problem(self):
        # initialize a system model
        SpanLength = 10.0 * 12
        Nelems = 8    # number of elements
        params = {'E': 29000., 'A': 4.7, 'I':103}

        # meshing parameters
        Le = SpanLength / Nelems
        Xnode = 0.0
        Ynode = 0.0

        model = System()

        # create left node
        nd0 = Node(Xnode, Ynode)
        model += nd0

        # initialization for node and element creation
        ndi = nd0

        for e in range(Nelems):
            # create next node
            Xnode += Le
            ndj = Node(Xnode, Ynode)
            model += ndj

            # remember center node for loading
            if Xnode <= SpanLength/2:
                ndP = ndj

            # create elements
            model += Beam2D(ndi, ndj, ElasticSection(params))

            # shift one node to the right
            ndi = ndj

        # define support(s)
        nd0.fixDOF('ux', 'uy')     # pin support left end
        ndj.fixDOF('uy')           # roller support right end

        # add point loads
        # .. load only the center node
        ndP.setLoad([0.0, -10.0], ('ux', 'uy'))

        # analyze the model
        model.solve()

        # write out report
        model.report()

        # create plots
        model.plot(factor=100.)

