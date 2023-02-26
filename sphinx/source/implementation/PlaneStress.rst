PlaneStress material class
============================

.. math::
    \left\{ \begin{array}{c}
        S_{XX} \\ S_{YY} \\ S_{XY}
        \end{array} \right\}
    =
    \frac{E t}{1-\nu^2}
    \left[ \begin{array}{ccc}
          1  & \nu  &   0        \\
        \nu  &   1  &   0        \\
          0  &   0  & \frac{1-\nu}{2}
    \end{array} \right]
    \left(
    \left\{ \begin{array}{c}
        E_{XX} \\ E_{YY} \\ 2E_{XY}
    \end{array} \right\}
    -
    \left\{ \begin{array}{c}
        E_{XX}^p \\ E_{YY}^p \\ 2E_{XY}^p
    \end{array} \right\}
    \right)

where

.. list-table::

    * - :math:`E_{IJ}`
      - Component IJ of total strain
    * - :math:`E_{IJ}^p`
      - Component IJ of plastic strain



Parent class
---------------
* :doc:`Material_class`

Class doc
-------------

.. automodule:: femedu.materials.PlaneStress
  :members:
