Truss Examples
=====================

.. toctree::
    :hidden:

    trusses/truss01/truss01.rst
    trusses/truss02/truss02.rst
    trusses/truss03/truss03.rst
    trusses/truss04/truss04.rst

.. list-table:: Available Examples
    :widths:  30 70
    :header-rows: 1

    * - Example
      - Description
    * - :doc:`trusses/truss01/truss01`
      - Basic truss triangle
    * - :doc:`trusses/truss02/truss02`
      - Simple statically determinate bridge
    * - :doc:`trusses/truss03/truss03`
      - Exploring alternative input format for **truss01**
    * - :doc:`trusses/truss04/truss04`
      - 3D truss example

How to run a truss example from the distribution
----------------------------------------------------

All Truss examples are packaged in :code:`examples.truss`.
To run a specific example use, e.g.:

.. code:: python

    from femedu import *
    from femedu.examples.trusses.truss01 import *

    # load the example
    ex = ExampleTruss01()

    # print the doc-string for the current example
    print(ex)

    # run the actual problem
    ex.run()


**More**: :doc:`../tutorials` and :ref:`examples-index`
