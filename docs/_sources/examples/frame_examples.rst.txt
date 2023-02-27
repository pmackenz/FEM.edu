Frame Examples
=====================

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


.. list-table:: Available Examples
    :widths:  30 70
    :header-rows: 1

    * - Example
      - Description
    * - **frame01**
      - Axially loaded beam with lateral load
    * - **frame02**
      - Axially loaded column with lateral load. Rotated **frame01** to test invariance.
    * - **frame03**
      - Building frame


