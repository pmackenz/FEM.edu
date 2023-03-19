Frame Examples
=====================

.. toctree::
    :hidden:

    frames/frame01/frame01.rst
    frames/frame02/frame02.rst
    frames/frame03/frame03.rst
    frames/frame04/frame04.rst

.. list-table:: Available Examples
    :widths:  30 70
    :header-rows: 1

    * - Example
      - Description
    * - :doc:`frames/frame01/frame01`
      - Axially loaded beam with lateral load
    * - :doc:`frames/frame02/frame02`
      - Axially loaded column with lateral load. Rotated **frame01** to test invariance.
    * - :doc:`frames/frame03/frame03`
      - Single story Building frame
    * - :doc:`frames/frame04/frame04`
      - Multi story Building frame

**More**: :doc:`../tutorials` and :ref:`examples-index`


How to run a frame example from the distribution
----------------------------------------------------

All frame examples are packaged in :code:`examples.frames`.
To run a specific example use, e.g.:

.. code:: python

    from femedu import *
    from femedu.examples.frames.frame01 import *

    # load the example
    ex = ExampleFrame01()

    # print the doc-string for the current example
    print(ex)

    # run the actual problem
    ex.run()

