Element
==============

Each instance represents one truss member

.. list-table:: Element class methods
   :widths: 25 25 25 50
   :header-rows: 1

   * - method
     - input
     - returns
     - description
   * - `__init__(nd0, nd1, material)`
     - two `Node()` objects, one `Material()` object.
     - 
     - constructor.
   * - `getForce()`
     - 
     - list of (1d) np.array objects
     - A list of nodal forces.  Each nodal force shall be represented as a 1d np.array with
       two components of the force.
   * - `getStiffness()`
     - 
     - 2-by-2 list of 2-by-2 np.array
     - A list of lists (matrix) containing nodal tangent matrices as `np.array([[.,.],[.,.]])`

**Note**: a `Node()` object may have changed its state between calls, so you need to
recompute every time!


.. list-table:: Element class variables
   :widths: 25 25 50
   :header-rows: 1

   * - name
     - type
     - description
   * - `nodes`
     - List of `Node()` instances
     - representing the two nodes at either end of a truss.
   * - `material`
     - `Material()`
     - pointer to an instance of `Material()`. Needed to compute stress and tangent modulus.
   * - `force`
     - 2-list of (1d) np.array objects.
     - holding nodal forces `P0` and `P1`
   * - `Kt`
     - 2-by-2 array of 2-by-2 np.array objects.
     - tangent stiffness matrix. Representing all nodal stiffness matrices.

**Equations**

1. :math:`{\bf L} = {\bf X}_1 - {\bf X}_0`
#. :math:`\ell = ||{\bf L}||`
#. :math:`{\bf n} = \frac{1}{\ell} \, {\bf L}`
#. Strain: :math:`\varepsilon = \frac{1}{\ell} \, {\bf n}\cdot( {\bf U}_1 - {\bf U}_0)`
#. Force: :math:`f = \sigma(\varepsilon) A` using `material.setStrain(eps)` and `material.getStress()`.
#. Nodal force vector: :math:`{\bf P}^e = f \, {\bf n}`
#. :math:`{\bf P}_0 = -{\bf P}^e ~~~~~~ {\bf P}_1 = {\bf P}^e`
#. Nodal stiffness matrix: :math:`{\bf k}^e = \frac{E_t(\varepsilon)\,A}{\ell}\, {\bf n}\otimes{\bf n}` using `material.getStiffness()` to find :math:`E_t`.
#. :math:`{\bf k}_{00} = {\bf k}_{11} = {\bf k}^e ~~~~~~~~ {\bf k}_{01} = {\bf k}_{10} = -{\bf k}^e`

