Material
==================

**This class is provided as a demonstration example.**


.. list-table:: Material class methods
   :widths: 25 25 25 50
   :header-rows: 1

   * - method
     - input
     - returns
     - description
   * - `__init__(...)`
     - parameters as `{'E':10.0}`
     - 
     - constructor. Sets parameters for this material and initializes all internal variables
   * - `getArea()`
     - 
     - :math:`A`
     - return cross section area from `parameters['A']`
   * - `getStress()`
     - 
     - :math:`\sigma`
     - request axial stress
   * - `getStiffness`
     - 
     - :math:`E_t`
     - request axial stiffness
   * - `setStrain(eps)`
     - strain :math:`\varepsilon`
     - 
     - update state for a user provided axial strain value

.. list-table:: Element class variables
   :widths: 25 25 50
   :header-rows: 1

   * - name
     - type
     - description
   * - `params`
     - dict
     - default parameters: `{'E':100., 'nu':0.0,  'fy:1.0e30}`
       Holds user provided parameters (MOE, Poisson's ratio, yield stress)
   * - `plastic_strain`
     - float
     - internal state variable.
   * - `sig`
     - float
     - holds current stress
   * - `Et`
     - float
     - holds current materil tangent modulus


**Equations**

1. Elastic trial state:

   1.1  :math:`\sigma = E * (\varepsilon - \varepsilon_P)`
   
   1.2  :math:`E_t = E`

2. Yield check: :math:`f = ||\sigma|| - f_y`


3. IF :math:`f \ge 0`:

   3.1. :math:`\Delta\varepsilon_P = \text{sign}(\sigma) * \frac{f}{E}`
   
   3.2. :math:`\sigma = \sigma -  E * \Delta\varepsilon_P`
   
   3.3. :math:`E_t = E_t - E`


