"""
Example: frame01

Buckling of a beam with pin-pin support

modeled using a 2D frame element

.. code::

    x============x=============x  <---
    ^                          o

    x ..... node
    === ... frame element
    <-- ... applied force
    ^ ..... pin support
    o ..... roller support

    degrees of freedom:
    0 ... horizontal displacement, u
    1 ... vertical displacement, v
    2 ... rotation, theta


.. list-table:: setting given parameters

    * - N  = 2
      - number of elements
    * - L  = 100.0
      - column length
    * - EA = 2000000.0
      - axial stiffness
    * - EI = 21000.0
      - flexural stiffness
    * - w  = 0.1
      - applied lateral load


"""
from ...examples.Example import *

from ...domain.System import *
from ...solver.NewtonRaphsonSolver import *
from ...domain.Node import *
from ...elements.Frame2D import *
from ...materials.ElasticSection import *


class ExampleFrame01(Example):

    def docString(self):
        s = """
    Buckling of a beam with pin-pin support
    
    modeled using a 2D frame element
    
        x============x=============x  <---
        ^                          o
        
        x ..... node
        === ... frame element
        <-- ... applied force
        ^ ..... pin support
        o ..... roller support
        
    degrees of freedom:
        0 ... horizontal displacement, u
        1 ... vertical displacement, v
        2 ... rotation, theta
        
    parameters 
        
        N  = 2     # number of elements
        L  = 100.0
        EA = 2000000.0
        EI = 21000.0
        w  = 0.1
        
    Author: Peter Mackenzie-Helnwein 
    """
        return s

    def problem(self):
        # initialize a system model

        N  = 16     # number of elements
        L  = 100.0
        E  = 20000.
        EA = 2000000.0
        EI = 210000.0
        w  = -0.1

        params = {'E': E, 'A': EA/E, 'I': EI/E}

        model = System()
        model.setSolver(NewtonRaphsonSolver())

        # create nodes

        nd0 = Node(0.0, 0.0)
        model += nd0

        ndi = nd0
        for i in range(N):
            # nodes
            ndj = Node( (i+1)*L/N, 0.0 )
            model += ndj

            # elements
            elem = Frame2D(ndi, ndj, ElasticSection(params))
            elem.setDistLoad(w)
            model += elem

            ndi = ndj

        # define support(s)
        nd0.fixDOF('ux', 'uy')    # horizontal support left end
        ndi.fixDOF('uy')          # vertical support right end

        # add loads
        # .. load only the upper nodes
        Pcr = np.pi**2 * EI / L**2
        ndi.setLoad((-0.5*Pcr,), ('ux',))

        # show model information
        print(model)

        model.solve(verbose=True)

        model.report()

        model.plot(factor=10.0)

        model.beamValuePlot("F")
        model.beamValuePlot("M")
        model.beamValuePlot("V")

        return









        # ========== setting global parameters ==============

        target_load_level = 0.99
        max_steps = 50
        load_levels = np.linspace(0, target_load_level, max_steps)

        # ========= build your structural model =============




        print(model)

        # ===== applying the load in multiple smaller time steps ========

        # set up data recorder
        model.initRecorder()
        model.trackStability(True)

        # initialize the analysis:
        model.resetDisplacements()   # set U to all zeros
        model.setLoadFactor(0.0)     # define a known equilibrium solution

        model.startRecorder()

        # solve for all load_levels
        for loadfactor in load_levels:

            # define node X2 as the controled node; downward direction is prescribed:
            model.setLoadFactor(loadfactor)
            model.NewtonSolver()

            # report results
            print('+')
            model.report()

            print("\n=== next load level ===\n")

        model.plot([3*(N//2),3*(N//2)+1,3*(N//2)+2,3*N])

        model.plotSystem(1.0)

        model.plotBucklingMode(10.)
