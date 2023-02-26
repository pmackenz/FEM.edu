ElasticSection
==========================

The elastic section implements the simple linear elastic relations

.. math::

    f = EA \, \varepsilon
    \qquad    \qquad
    M = EI \, \phi

where

.. list-table::

    * - :math:`f`
      - resulting axial force (output)
    * - :math:`M`
      - resulting moment (output)
    * - :math:`\varepsilon`
      - axial strain (input)
    * - :math:`\phi`
      - curvature change (input)
    * - :math:`EA`
      - axial stiffness; constant (input parameters E and A)
    * - :math:`EI`
      - flexural stiffness; constant (input parameters E and I)


Parent class
---------------
* :doc:`Sections`


.. automodule:: femedu.materials.ElasticSection
  :members:

