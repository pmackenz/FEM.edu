Mixed Structures Examples
==========================

.. toctree::
    :hidden:

    mixed/mixed01/mixed01.rst
    mixed/mixed02/mixed02.rst

.. list-table:: Available Examples
    :widths:  30 70
    :header-rows: 1

    * - Example
      - Description
    * - **mixed01**
      - WIP
    * - **mixed02**
      - WIP

**More**: :doc:`../tutorials` and :ref:`examples-index`

How to run a mixed structures example from the distribution
--------------------------------------------------------------

All mixed structures examples are packaged in :code:`examples.mixed`.
To run a specific example use, e.g.:

.. code:: python

    from femedu import *
    from femedu.examples.mixed.mixed01 import *

    # load the example
    ex = ExampleMixed01()

    # print the doc-string for the current example
    print(ex)

    # run the actual problem
    ex.run()

