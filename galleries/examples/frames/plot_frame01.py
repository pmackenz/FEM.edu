"""
======================================================
Buckling of a beam with pin-pin support
======================================================

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

Author: Peter Mackenzie-Helnwein
"""

from femedu.examples.Example import *

from femedu.domain import *
from femedu.solver.NewtonRaphsonSolver import *
from femedu.elements.finite.Frame2D import *
from femedu.materials.ElasticSection import *


class ExampleFrame01(Example):

    # sphinx_gallery_start_ignore
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

    # sphinx_gallery_end_ignore
    def problem(self):

        #
        # ==== Initialization ====
        #

        # ========== setting mesh parameters ==============

        N = 8         # number of elements in the mesh
        L = 100.0     # column free length


        # ========== setting material parameters ==============

        params = dict(
            E = 20000.,    # Young's modulus
            A = 100.0,     # cross section area
            I = 10.0       # cross section moment of inertia
        )

        # ========== setting load parameters ==============

        w   = -0.1         # uniform lateral load on the column
        Pcr = np.pi**2 * params['E'] * params['I'] / L**2    # Euler buckling load

        # ========== setting analysis parameters ==============

        target_load_level = 0.99      # 99% of Euler load
        max_steps = 10                # solve max_steps points on the primary path


        w   *= 0.01
        Pcr *= 0.01
        target_load_level = 99.      # 99% of Euler load



        # define a list of target load levels
        load_levels = np.linspace(0, target_load_level, max_steps)

        #
        # ==== Build the system model ====
        #

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
            model += elem

            # ** apply the element portion of the reference load
            elem.setDistLoad(w)

            ndi = ndj    # jump to next element: make current end-node the next start-node

        # define support(s)

        nd0.fixDOF('ux', 'uy')    # horizontal support left end
        ndi.fixDOF('uy')          # vertical support right end

        # ==== complete the reference load ====

        # these are only nodal forces as part of the reference load
        # .. load only the upper node
        ndi.setLoad((-Pcr,), ('ux',))

        # show model information
        print(model)

        #
        # ==== perform the analysis ===
        #

        print("\n==== perform the analysis ===\n")

        # * apply the load in multiple smaller load steps

        # set up data recorder
        model.initRecorder()
        model.trackStability(True)

        # initialize the analysis:
        model.resetDisplacements()   # set U to all zeros
        model.setLoadFactor(0.0)     # define a known equilibrium solution

        model.startRecorder()

        detKt   = []
        lambdas = []

        # solve for all load_levels
        for loadfactor in load_levels:

            # define node X2 as the controled node; downward direction is prescribed:
            model.setLoadFactor(loadfactor)
            model.solve(verbose=True)

            # stability check
            lambdas.append(model.loadfactor)
            detKt.append(model.solver.checkStability())

            # report results
            print('+')
            #model.report()

            print("\n=== next load level ===\n")


        #
        # ==== create some nice plots ===
        #

        model.report()

        model.plot(factor=1.0, filename="frame1_deformed.png")

        fig, ax = plt.subplots()

        ax.plot(lambdas,detKt,'--*r')
        ax.grid(True)
        ax.set_xlabel('Load factor, $ \lambda $')
        ax.set_ylabel("Stability index, $ {det}\: {\\bf K}_t $")

        fig.savefig("frame1_stability.png")
        fig.show()

        model.beamValuePlot("F", filename="frame1_force.png")
        model.beamValuePlot("V", filename="frame1_shear.png")
        model.beamValuePlot("M", filename="frame1_moment.png")

        model.plotBucklingMode(factor=10., filename="frame1_buckling_mode0.png")


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleFrame01()
    ex.run()

