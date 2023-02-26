Truss Examples
=====================

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


.. list-table:: Available Examples
    :widths:  30 70
    :header-rows: 1

    * - Example
      - Description
    * - **truss01**
      - Basic truss triangle
    * - **truss02**
      - Simple statically determinate bridge
    * - **truss03**
      - Exploring alternative input format for **truss01**
    * - **truss04**
      - 3D truss example


