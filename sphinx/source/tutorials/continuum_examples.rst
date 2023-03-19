Continuum Examples
=====================

.. toctree::
    :hidden:

    solids/solid01/solid01.rst

.. list-table:: Available Examples
    :widths:  30 70
    :header-rows: 1

    * - Example
      - Description
    * - **solid01**
      - WIP

**More**: :doc:`../tutorials` and :ref:`examples-index`

How to run a solid model example from the distribution
--------------------------------------------------------

All continuum examples are packaged in :code:`examples.solids`.
To run a specific example use, e.g.:

.. code:: python

    from femedu import *
    from femedu.examples.solids.solid01 import *

    # load the example
    ex = ExampleSolid01()

    # print the doc-string for the current example
    print(ex)

    # run the actual problem
    ex.run()
