System
===========

Creates an instance of a truss model

.. list-table:: System class methods
   :widths: 25 25 25 50
   :header-rows: 1

   * - method
     - input
     - returns
     - description
   * - `__init__(...)`
     - 
     - 
     - constructor.
   * - `addNode(newNode)`
     - `Node(...)` object
     - 
     - add one `Node()` object to your list of elements (the model)
   * - `addElement(newElem)`
     - `Element(...)` object
     - 
     - add one `Element()` object to your list of elements (the model)
   * - `solve()`
     - 
     - 
     - assemble :math:`[K_t]` and :math:`\{P\}`, solve for :math:`\{u\} = [K_t]^{-1}\{P\}`,
       loop through nodes and update nodal displacement, compute unbalanced force :math:`\{R\}
       = \{P\} - \{F\}`
   * - `plot(factor=1.0)`
     - 
     - 
     - collect node info and send it to the plotter. Request the plot.
   * - `report()`
     - 
     - 
     - print a summary report: list of nodal position, load, displacement, unbalanced force.


.. list-table:: System class variables
   :widths: 25 25 50
   :header-rows: 1

   * - name
     - type
     - description
   * - `nodes`
     - List of `Node()` objects
     - holds all the nodes in the model
   * - `elements`
     - List of `Element()` objects
     - holds all the elements in the model
   * - `solver`
     - `Solver()`
     - pointer to `Solver()` object to handle plotting
   * - `plotter`
     - `Plotter()`
     - pointer to `Plotter()` object to handle plotting

