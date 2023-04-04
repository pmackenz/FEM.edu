"""
======================================================
Buckling of a vertical beam with pin-pin support
======================================================

modeled using a 2D frame element

.. code::

    |
    | P
    v
    x o ... support
    I
    x   ... node
    I
    x   ... node
    ^   ... support

    x ..... node
    I ... frame element
    <-- ... applied force
    ^ ..... pin support
    o ..... roller support

    degrees of freedom:
    0 ... horizontal displacement
    1 ... vertical displacement
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


class ExampleFrame02(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = """
    Buckling of a beam with pin-pin support
    
    modeled using a 2D frame element
        
        |
        | P
        v
        x o ... support
        I
        x   ... node
        I
        x   ... node
        ^   ... support
    
        x ..... node
        I ... frame element
        <-- ... applied force
        ^ ..... pin support
        o ..... roller support

    degrees of freedom:
        ux ... horizontal displacement
        uy ... vertical displacement
        rz ... rotation, theta
        
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
        # initialize a system model

        N  = 8     # number of elements
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
            ndj = Node( 0.0, (i+1)*L/N )
            model += ndj

            # elements
            elem = Frame2D(ndi, ndj, ElasticSection(params))
            elem.setDistLoad(w)
            model += elem

            ndi = ndj

        # define support(s)
        nd0.fixDOF('ux', 'uy')    # horizontal and vertical support bottom end
        ndi.fixDOF('ux')          # horizontal support top end

        # add loads
        # .. load only the upper nodes
        Pcr = np.pi**2 * EI / L**2
        ndi.setLoad((-0.5*Pcr,), ('uy',))

        # show model information
        print(model)

        model.solve(verbose=True)

        model.report()

        model.plot(factor=10.0, filename="frame2_deformed.png")

        model.beamValuePlot("F", filename="frame2_force.png")
        model.beamValuePlot("V", filename="frame2_shear.png")
        model.beamValuePlot("M", filename="frame2_moment.png")

        return


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleFrame02()
    ex.run()


