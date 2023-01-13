Node
=============

Each instance represents one node in the system

.. list-table:: Node class methods
   :widths: 25 25 25 50
   :header-rows: 1

   * - method
     - input
     - returns
     - description
   * - `__init__(x,y)`
     - coordinates of the point as two floats
     - 
     - constructor. Sets position and initializes displacement and force to zeros.
   * - `fixDOF(idx)`
     - `idx` of the degree of freedom (dof)
     - 
     - set internal flag for this dof accordingly.
   * - `isFixed(idx)`
     - `idx` of the degree of freedom (dof)
     - True|False
     - test function returning True if dof at `idx` is fixed, False otherwise.
   * - `setDisp(u,v)`
     - components of displacement
     - 
     - overwrited the displacements for this node.
   * - `getDisp()`
     - 
     - `np.array([u,v])`
     - returns displacement vector
   * - `getPos()`
     - 
     - `np.array([x,y])`
     - returns initial position
   * - `getDeformedPos(factor)`
     - 
     - `np.array([x+factor*u,y+factor*v])`
     - returns current position with displacement magnified by factor.  Would be good to have
       a default factor of 1.0 if none given.
   * - `addLoad(Px,Py)`
     - components of load
     - 
     - add this load to nodal load
   * - `setLoad(Px,Px)`
     - components of load
     - 
     - replace current load by provided load
   * - `getLoad()`
     - 
     - `np.array([Px,Py])`
     - returns current load


.. list-table:: Node class variables
   :widths: 25 25 50
   :header-rows: 1

   * - name
     - type
     - description
   * - `pos`
     - `np.array([x,y])`
     - holds `x` and `y` coordinates of the points
   * - `index`
     - `int`
     - index position of this `Node()` in a `System().nodes` list.  This needs to be set by
       `System().addNode(thisNode)`, so coordinate this with the `System()` team.
   * - `disp`
     - `np.array([u,v])`
     - holds `x` and `y` components of nodal displacement
   * - `fixity`
     - list of two `True|False`
     - `fixity[i]` is `True` if `i`-th degree of freedom is fixed, `False` otherwise.  Note:
       `i=0|1`
   * - `force`
     - `np.array([Px,Py])`
     - holds `x` and `y` components of the
       nodal load vector.


**Equations**

1. Deformed position node 0: :math:`{\bf x}_0 = {\bf X}_0 + (factor)*{\bf U}_0`
#. Deformed position node 1: :math:`{\bf x}_1 = {\bf X}_1 + (factor)*{\bf U}_1`

