Beam Examples
=====================

.. toctree::
    :hidden:

    beams/beam01/beam01.rst
    beams/beam02/beam02.rst

.. list-table:: Available Examples
    :widths:  30 70
    :header-rows: 1

    * - Example
      - Description
    * - :doc:`beams/beam01/beam01`
      - Single span beam with point load
    * - :doc:`beams/beam02/beam02`
      - Three-span continuous beam with point loads

**More**: :doc:`../tutorials` and :ref:`examples-index`

How to run a beam example from the distribution
----------------------------------------------------

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


