Example: beam02
==================

We build the model based a few parameters as follows.

.. literalinclude:: ../../../../../src/femedu/examples/beams/beam02.py
   :lineno-start: 1
   :lines: 26-32

`SpanLengths` is a list of span lengths for a series of spans on a continuous beam.
Each span will be modeled by `Nelems` elements.
This list serves of a wrapper around input similar to :doc:../beam01/beam01.rst.
The loop is emphasized in the next code segment (yellow line).

All mesh creation is based solely on the above parameters to allow for easy
manipulation of the model.

The actual model is built by the block below.

.. literalinclude:: ../../../../../src/femedu/examples/beams/beam02.py
   :lineno-start: 8
   :lines: 34-75
   :emphasize-lines: 16

Line 8 instantiates one model space.

Lines 16, 28-30 create the nodes, and
lines 18 and 31 add them to the model space.

Lines 34-35 create the elements and  add them to the model space.
You only need to create variables for `Node` and `Element` objects, respectively,
if you need to either add or retrieve information from that object later.

A uniform load `w` (positive is **up**) is applied directly to the beam elements
in line 38.

The support conditions are defined by providing the respective information
directly to the supported nodes.
Line 17 adds a pin support on the first node and line 44 adds a roller on the last one.


The system equations are solved by a single call to the solver:

.. literalinclude:: ../../../../../src/femedu/examples/beams/beam02.py
   :lineno-start: 50
   :lines: 77-78

You can obtain a debug-style report on the state of the system:

.. literalinclude:: ../../../../../src/femedu/examples/beams/beam02.py
   :lineno-start: 52
   :lines: 80-81

Resulting in an output like (may change as the code evolves).

    .. code-block:: text

        System Analysis Report
        =======================

        Nodes:
        ---------------------
          Node 0: {'uy': 0, 'rz': 1}
                  x:[0. 0.], fix:['ux', 'uy'],
                  P:[0. 0.], u:[ 0.        -0.0123663]
          Node 1: {'uy': 0, 'rz': 1}
                  x:[48.  0.], fix:[],
                  P:[0. 0.], u:[-0.292646    0.00326429]
          Node 2: {'uy': 0, 'rz': 1}
                  x:[96.  0.], fix:['uy'],
                  P:[0. 0.], u:[ 0.         -0.00069085]
          Node 3: {'uy': 0, 'rz': 1}
                  x:[156.   0.], fix:[],
                  P:[0. 0.], u:[-3.93139430e-01 -4.24063967e-19]
          Node 4: {'uy': 0, 'rz': 1}
                  x:[216.   0.], fix:['uy'],
                  P:[0. 0.], u:[0.         0.00069085]
          Node 5: {'uy': 0, 'rz': 1}
                  x:[264.   0.], fix:[],
                  P:[0. 0.], u:[-0.292646   -0.00326429]
          Node 6: {'uy': 0, 'rz': 1}
                  x:[312.   0.], fix:['uy'],
                  P:[0. 0.], u:[0.        0.0123663]

        Elements:
        ---------------------
          Beam2D: node 0 to node 1:
             material ElasticSection properties: {'E': 29000.0, 'A': 5, 'I': 50, 'nu': 0.0, 'fy': 1e+30}  strain:{'axial': 0.0, 'flexure': 0.0}   stress:{'axial': 0.0, 'flexure': 0.0}
             nodal forces: Vi:11.673913043478258 Mi:-191.99999999999994 Vj:-11.673913043478258 Mj:752.3478260869563
          Beam2D: node 1 to node 2:
             material ElasticSection properties: {'E': 29000.0, 'A': 5, 'I': 50, 'nu': 0.0, 'fy': 1e+30}  strain:{'axial': 0.0, 'flexure': 0.0}   stress:{'axial': 0.0, 'flexure': 0.0}
             nodal forces: Vi:-36.32608695652172 Mi:-752.347826086956 Vj:36.32608695652172 Mj:-991.3043478260865
          Beam2D: node 2 to node 3:
             material ElasticSection properties: {'E': 29000.0, 'A': 5, 'I': 50, 'nu': 0.0, 'fy': 1e+30}  strain:{'axial': 0.0, 'flexure': 0.0}   stress:{'axial': 0.0, 'flexure': 0.0}
             nodal forces: Vi:29.999999999999993 Mi:883.3043478260867 Vj:-29.999999999999993 Mj:916.6956521739128
          Beam2D: node 3 to node 4:
             material ElasticSection properties: {'E': 29000.0, 'A': 5, 'I': 50, 'nu': 0.0, 'fy': 1e+30}  strain:{'axial': 0.0, 'flexure': 0.0}   stress:{'axial': 0.0, 'flexure': 0.0}
             nodal forces: Vi:-29.99999999999999 Mi:-916.6956521739127 Vj:29.99999999999999 Mj:-883.3043478260865
          Beam2D: node 4 to node 5:
             material ElasticSection properties: {'E': 29000.0, 'A': 5, 'I': 50, 'nu': 0.0, 'fy': 1e+30}  strain:{'axial': 0.0, 'flexure': 0.0}   stress:{'axial': 0.0, 'flexure': 0.0}
             nodal forces: Vi:36.32608695652173 Mi:991.3043478260865 Vj:-36.32608695652173 Mj:752.347826086956
          Beam2D: node 5 to node 6:
             material ElasticSection properties: {'E': 29000.0, 'A': 5, 'I': 50, 'nu': 0.0, 'fy': 1e+30}  strain:{'axial': 0.0, 'flexure': 0.0}   stress:{'axial': 0.0, 'flexure': 0.0}
             nodal forces: Vi:-11.673913043478258 Mi:-752.347826086956 Vj:11.673913043478258 Mj:192.0


An easier way to look at the simulation results are plots.  A deformed system plot is obtained
using the `model.plot()` directive.  If a `filename` is given, the plot will be saved
to the harddrive using that file name.
An internal force plot is created equally simple.

.. literalinclude:: ../../../../../src/femedu/examples/beams/beam02.py
   :lineno-start: 54
   :lines: 83-87

.. figure:: beam02_deformed.png
    :align: center

    Showing file *beam02_deformed.png*

.. figure:: beam02_shear.png
    :align: center

    Showing file *beam02_shear.png*

.. figure:: beam02_moment.png
    :align: center

    Showing file *beam02_moment.png*



**Importing the example**

.. code:: python

    from femedu.examples.beams.beam02 import *

    # load the example
    ex = Examplebeam02()

**More beam examples**: :doc:`../../beam_examples`
