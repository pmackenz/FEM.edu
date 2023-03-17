FiberSection
==========================

.. warning::

    This material is not yet available

The fiber section model implements a complex cross section
built from multiple :doc:`FiberMaterial` instances

.. math::

    \varepsilon(z) = \varepsilon_0 - z \,\phi

.. math::

    f = \int_A \sigma(z)\, dA
    \qquad    \qquad
    M = -\int_A z\sigma(z)\, dA

.. math::

    df = \int_A \mathcal{C}(z)\, dA \, d\varepsilon_0
        -\int_A z\mathcal{C}(z)\, dA \, d\phi
      =: \mathbb{D} \, d\varepsilon_0
        +\mathbb{C} d\phi

.. math::

    dM = -\int_A z\mathcal{C}(z)\, dA \, d\varepsilon_0
        +\int_A z^2\mathcal{C}(z)\, dA \, d\phi
      =: \mathbb{C} \, d\varepsilon_0
        + \mathbb{B} \, d\phi

where

.. list-table::

    * - :math:`\varepsilon_0`
      - axial strain (input)
    * - :math:`\phi`
      - curvature change (input)
    * - :math:`f`
      - resulting axial force (output)
    * - :math:`M`
      - resulting moment (output)
    * - :math:`\mathcal{C}(z)`
      - material tangent stiffness of a fiber-material at distance *z*
    * - :math:`\mathbb{D}`
      - axial stiffness
    * - :math:`\mathbb{B}`
      - flexural stiffness
    * - :math:`\mathbb{C}`
      - coupling stiffness


Parent class
---------------
* :doc:`Sections`


.. . automodule:: femedu.materials.ElasticSection
  :members:

