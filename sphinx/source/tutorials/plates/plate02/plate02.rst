Example: plate02
==================

.. figure:: plate02_undeformed.png
    :scale: 75%
    :align: center

    Initial system and meshing for the patch test.

We build the model based a few parameters as follows.

.. literalinclude:: ../../../../../src/femedu/examples/plates/plate02.py
   :lineno-start: 1
   :lines: 73-81

All mesh creation is based solely on the above parameters to allow for easy
manipulation of the model.

The actual model is built by the block below.

.. literalinclude:: ../../../../../src/femedu/examples/plates/plate02.py
   :lineno-start: 10
   :lines: 83-102

Line 10 instantiates one model space.

Line 11 switches from the default linear solver to the :py:class:`NewtonRaphsonSolver`, needed for nonlinear problems.

Lines 13-16 create the nodes, and
lines 18 adds them to the model space.

Lines 18-20 apply the geometric boundary condition.  This examples implements symmetry conditions at :math:`x=0` and :math:`y=0`.

Lines 24-25 create the elements and line 27 adds them to the model space.
You only need to create variables for `Node` and `Element` objects, respectively,
if you need to either add or retrieve information from that object later.

Lines 29 adds a surface load to face 2 of element :py:obj:`elemB`. See :doc:`../../../implementation/Elements/Triangle_class`
for the definition of faces for this element.

We initialize the system by solving for load-level 0.00.
This is not necessary for most models, though some thermo-elastic or elast-plastic problems, e.g., may not be in equilibrium in an undeformed state.
The cost is minimal and adding this step to a nonlinear analysis is a good habit.
The system equations are solved by a single call to the solver:

.. literalinclude:: ../../../../../src/femedu/examples/plates/plate02.py
   :lineno-start: 30
   :lines: 104-109

Next increase the load to the next target load-level, here we are using 1.0 (100% of the reference load).
The system equations are once again solved by a single call to the solver:

.. literalinclude:: ../../../../../src/femedu/examples/plates/plate02.py
   :lineno-start: 36
   :lines: 111-113

The solver report the residual norm during each step. The geometric nonlinearity of the element
requires multiple iteration steps before achieving equilibrium.  The following convergence print-out
shows that the element generates a quadratic rate of convergence when paired with the :py:class:`NewtonRaphsonSolver`.

    .. code::

        norm of the out-of-balance force:   7.0711e+00
        norm of the out-of-balance force:   1.1554e+00
        norm of the out-of-balance force:   1.7307e-02
        norm of the out-of-balance force:   4.2408e-06
        norm of the out-of-balance force:   2.2995e-13


.. figure:: plate02_deformed.png
    :align: center

    Deformed system at load level 1.00


You can obtain a debug-style report on the state of the system:

.. literalinclude:: ../../../../../src/femedu/examples/plates/plate02.py
   :lineno-start: 39
   :lines: 115

The report looks like this:

    .. code::

            System Analysis Report
            =======================

            Nodes:
            ---------------------
              Node_0:
                  x:    [0. 0.]
                  fix:  ['ux', 'uy']
                  u:    [0. 0.]
              Node_1:
                  x:    [10.  0.]
                  fix:  ['uy']
                  u:    [0.88033915 0.        ]
              Node_2:
                  x:    [10. 10.]
                  u:    [ 0.88033915 -0.27963653]
              Node_3:
                  x:    [ 0. 10.]
                  fix:  ['ux']
                  u:    [ 0.         -0.27963653]

            Elements:
            ---------------------
              LinearTriangle: nodes ( Node_0 Node_1 Node_3 )
                  material: PlaneStress
                  strain: xx=9.191e-02 yy=-2.757e-02 xy=0.000e+00 zz=-1.930e-02
                  stress: xx=9.191e-01 yy=2.054e-15 xy=0.000e+00 zz=0.000e+00
              LinearTriangle: nodes ( Node_2 Node_3 Node_1 )
                  material: PlaneStress
                  strain: xx=9.191e-02 yy=-2.757e-02 xy=1.727e-16 zz=-1.930e-02
                  stress: xx=9.191e-01 yy=3.886e-15 xy=6.641e-16 zz=0.000e+00
                  element forces added to node:
                      Node_2: [5. 0.]
                      Node_3: [0. 0.]
                      Node_1: [5. 0.]


**Importing the example**

.. code:: python

    from femedu.examples.plates.plate02 import *

    # load the example
    ex = ExamplePlate02()

**More frame examples**: :doc:`../../plate_examples`
