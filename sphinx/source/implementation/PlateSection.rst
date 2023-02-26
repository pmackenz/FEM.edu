PlateSection
==========================

.. warning::

    This material is not yet available

The plate section model implements a complex cross section
built from multiple :doc:`PlaneStress` instances

.. math::

    \varepsilon_{kl}(z) = \varepsilon^0_{kl} - z \,\phi_{kl}

.. math::

    n^{ij} = \int_{-h/2}^{h/2} \sigma^{ij}(z)\, dz
    \qquad    \qquad
    m^{ij} = -\int_{-h/2}^{h/2} z\sigma^{ij}(z)\, dz


Stress :math:`\sigma^{ij}(z)` is obtained from a :doc:`PlaneStress` material layer at distance *z*.
Through-the-thickness integration is performed numerically as a sum over layers.

The **PlateSection** material will further compute the tangent stiffness tensors defined as:

.. math::

    dn^{ij}
      =: \mathbb{D}^{ijkl} \, d\varepsilon^{0}_{kl}
        +\mathbb{C}^{ijkl} d\phi_{kl}

.. math::

    dm^{ij}
      =: \mathbb{C}^{ijkl} \, d\varepsilon^{0}_{kl}
        + \mathbb{B}^{ijkl} \, d\phi_{kl}

where

.. list-table::

    * - :math:`\varepsilon^0_{ij}`
      - membrane strain (input)
    * - :math:`\phi_{ij}`
      - curvature change (input)
    * - :math:`n^{ij}`
      - component of resulting membrane force (output)
    * - :math:`m^{ij}`
      - component of resulting moment (output)
    * - :math:`\mathcal{C}^{ijkl}(z)`
      - material tangent stiffness of a plane-stress material layer at distance *z*
    * - :math:`\mathbb{D}^{ijkl}`
      - axial stiffness
    * - :math:`\mathbb{B}^{ijkl}`
      - flexural stiffness
    * - :math:`\mathbb{C}^{ijkl}`
      - coupling stiffness

