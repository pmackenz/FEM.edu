Beam Examples
=====================

All Plate examples are packaged in :code:`examples.beams`.
To run a specific example use, e.g.:

.. code:: python

    from femedu import *
    from femedu.examples.beams.beam01 import *

    # load the example
    ex = ExampleBeam01()

    # print the doc-string for the current example
    print(ex)

    # run the actual problem
    ex.run()


.. list-table:: Available Examples
    :widths:  30 70
    :header-rows: 1

    * - Example
      - Description
    * - **beam01**
      - Single span beam with point load
    * - **beam02**
      - Three-span continuous beam with point loads



