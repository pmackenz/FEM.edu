"""
Example
"""
from ...examples.Example import *

from ...domain.System import *
from ...domain.Node import *
from ...elements.Beam2D import *
from ...materials.ElasticSection import *


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
        w =  -1.0   # distributed load (positive if acting in local y-direction
        P =   -40.0   # center point load (uses global system)

        Nelems = 4    # number of elements
        params = {'E': 29000., 'A': 4.7, 'I':103}

        model = System()

        # meshing parameters
        Le = SpanLength / Nelems
        Xnode = 0.0
        Ynode = 0.0

        # create left node
        nd0 = Node(Xnode, Ynode)
        model += nd0

        ndP = None

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
            elem = Beam2D(ndi, ndj, ElasticSection(params))
            model += elem

            # load the element
            elem.setDistLoad(w)

            # shift one node to the right
            ndi = ndj

        # define support(s)
        nd0.fixDOF('ux', 'uy')     # pin support left end
        ndj.fixDOF('uy')           # roller support right end

        # add point loads
        # .. load only the center node
        if ndP:
            ndP.setLoad([0.0, P], ('ux', 'uy'))

        # analyze the model
        model.solve()

        # write out report
        model.report()

        # create plots
        model.plot(factor=10., filename="beam01_deformed.png")

        model.beamValuePlot('V', filename="beam01_shear.png")
        model.beamValuePlot('M', filename="beam01_moment.png")

