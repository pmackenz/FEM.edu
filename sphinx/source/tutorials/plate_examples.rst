Plate Examples
=====================

.. toctree::
    :hidden:

    plates/plate01/plate01.rst
    plates/plate02/plate02.rst
    plates/plate03/plate03.rst
    plates/plate04/plate04.rst

.. list-table:: Available Examples
    :widths:  30 70
    :header-rows: 1

    * - Example
      - Description
    * - :doc:`plates/plate01/plate01`
      - Testing the internal force code
    * - :doc:`plates/plate02/plate02`
      - Testing the convergence behavior of a mini-problem
    * - :doc:`plates/plate03/plate03`
      - Patch test for in-plane loading.


**More**: :doc:`../tutorials` and :ref:`examples-index`

How to run a plate example from the distribution
----------------------------------------------------

All Plate examples are packaged in :code:`examples.plates`.
To run a specific example use, e.g.:

.. code:: python

    from femedu import *
    from femedu.examples.plates.plate01 import *

    # load the example
    ex = ExamplePlate01()

    # print the doc-string for the current example
    print(ex)

    # run the actual problem
    ex.run()


