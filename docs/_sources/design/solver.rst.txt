Solver
===========

Creates an instance of a Solver

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
   * - `disp`
     - `np.array([...])`
     - system sized displacement vector
   * - `loads`
     - `np.array([...])`
     - system sized load vector


**Equations**

1. each element has `node0 = elem.nodes[0]` and `node1 = elem.nodes[1]`
#. node indices `i = node0.index` and `j = node1.index`
#. a local d.o.f. :math:`u\to k=0` or :math:`v\to k=1` of node :math:`i` belongs at global index :math:`K = 4*i + k`
#. a local d.o.f. :math:`u\to m=0` or :math:`v\to m=1` of node :math:`j` belongs at global index :math:`M = 4*j + m`

*Assembly*:

* :math:`{\bf F}`: element force from `elem.getForce()`
* :math:`{\bf K_t}`: element stiffness from `elem.getStiffness()`
* :math:`{\bf R}_{sys}`: system force
* :math:`{\bf K_t}_{sys}`: system stiffness

5. Loop over `nodes`: :math:`{\bf R}_{sys}[K] = nodes[i].getLoad()[k]`   (this should return 0 if no load at this node and d.o.f.)
#. Loop over `elements`: :math:`{\bf R}_{sys}[K] = {\bf R}_{sys}[K] -  {\bf F}[i][k]`
#. Loop over `elements`: :math:`{\bf K_t}_{sys}[K,M] = {\bf K_t}_{sys}[K,M] +  {\bf K_t}[i,j][k,m]`
#. Loop over `nodes`: if a d.o.f. at *K* is fixed, set :math:`{\bf R}_{sys}[K] = 0` and :math:`{\bf K_t}_{sys}[K,K] = 1.0e20`.
#. Solve system of equations:  :math:`{\bf U} = {\bf K_t}_{sys}^{-1}\,{\bf R}_{sys}`
#. Loop over `nodes`: :math:`nodes[i].setDisp(u,v)` using :math:`u = {\bf U}[2*i]` and :math:`v = {\bf U}[2*i+1]`
#. Recompute: :math:`{\bf R}_{sys}` as in steps 5 and 6 (**do not repeat steps 7-11**).  If
   everything was done correctly, fixed d.o.f.s will contain the support reactions and free
   d.o.f.s should hold numeric zeros.

