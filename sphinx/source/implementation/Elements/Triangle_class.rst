.. _triangle_class:

Triangle class
==========================

Coordinate-free formulation for a finite deformation bi-linear triangle.

.. dropdown:: Linear Theory (small displacements)

    The element uses an  iso-parametric mapping for both the reference and the spatial configuration.
    In the reference configuration,

    .. math::
        :label: eq-triangle_lin1

        \mathbf{X}=\mathbf{X}_I \xi^I \quad I=0,1,2\equiv u,s,t \quad \xi^I \ge 0 ~\&~ \xi^0+\xi^1+\xi^2 = 1

    where :math:`\xi^I` are standard triangle coordinates.
    Alternatively, the mapping may be expressed using just two coordinates, :math:`\xi^1\equiv\xi^s` and
    :math:`\xi^2\equiv\xi^t`, as

    .. math::
        :label: eq-triangle_lin2

        \mathbf{X}=\mathbf{X}_0+\mathbf{G}_i \xi^i \quad i=1,2 \quad \xi^i \ge 0 ~\&~ \xi^1+\xi^2 \leq 1

    The covariant base vectors follow from :eq:`eq-triangle_lin1` and :eq:`eq-triangle_lin2` as

    .. math::
        :label: eq-triangle_lin3

        \mathbf{G}_i = \mathbf{X}_{,i} = \mathbf{X}_i - \mathbf{X}_0
        \quad i=1,2

    We can express differential :math:`d\xi^i` in terms of :math:`d\mathbf{X}` by means of the dual base as

    .. math::
        d\mathbf{X}=\mathbf{G}_i d \xi^i \rightarrow d \xi^i=\mathbf{G}^i \cdot d \mathbf{X} \quad w/ ~~ \mathbf{G}^i \mathbf{G}_j=\delta_j^i

    The displacement field is represented as

    .. math::
        \mathbf{u} = \mathbf{u}_I \xi^I
        \quad I=0,1,2

    Using :eq:`eq-triangle_lin3` we obtain the displacement gradient as

    .. math::
        d\mathbf{u}=\frac{\partial\mathbf{u}}{\xi^i} d \xi^i=\left(\mathbf{u}_I \otimes \mathbf{G}^I\right) \cdot d \mathbf{X}
        =: \nabla\mathbf{u} \cdot d\mathbf{\bar{X}} \Longrightarrow \nabla\mathbf{u}=\mathbf{u}_I \otimes \mathbf{G}^I

    where :math:`d\xi^0+d\xi^1+d\xi^2=0` and :math:`\mathbf{G}^0=-\mathbf{G}^1-\mathbf{G}^2` were used.


    The strain tensors follows as

    .. math::
        \symbf{\varepsilon} = \frac{1}{2} (\nabla\mathbf{u}^t + \nabla\mathbf{u})
                    = \frac{1}{2} (\mathbf{G}^I \otimes \mathbf{u}_I + \mathbf{u}_I \otimes \mathbf{G}^I)

    The variation of strain is needed for the *principle of virtual displacements*. It  follows as

    .. math::
        \delta \symbf{\varepsilon}
            = \frac{1}{2} (\mathbf{G}^I \otimes \delta \mathbf{u}_I + \delta \mathbf{u}_I \otimes \mathbf{G}^I)
                                                   \qquad I=0,1,2


    **Application of virtual work**

    The *principle of virtual displacements* (PVD) takes the form

    .. math::

        \int\!\!\!\int\limits_{\mathcal{A}} \delta \symbf{\varepsilon}:\symbf{\sigma}\: d\mathcal{A}
        = \int\!\!\!\int\limits_{\mathcal{A}} \delta\mathbf{u}\cdot\mathbf{b} \: d\mathcal{A}
        + \int\limits_{\partial\mathcal{A}_\sigma} \delta\mathbf{u}\cdot\bar{\mathbf{T}} \: d\mathcal{s}

    where :math:`\symbf{\sigma}=\mathbb{C}:d\symbf{\varepsilon}` is the Cauchy stress tensor,
    and :math:`\mathbb{C}` is the material stiffness tensor.

    The left hand term in the PVD yields

    .. math::

        \delta\symbf{\varepsilon}:\symbf{\sigma}
            = \frac{1}{2} (\mathbf{G}^I \otimes \delta \mathbf{u}_I + \delta \mathbf{u}_I \otimes \mathbf{G}^I):\symbf{\sigma}
            = \delta \mathbf{u}_I \cdot \symbf{\sigma} \cdot\mathbf{G}^I
            = \delta \mathbf{u}_I \cdot \mathbf{T}^I

    where :math:`\mathbf{T}^I= \symbf{\sigma} \cdot \mathbf{G}^I` is the traction on the edge opposing node *I*.

    Expressing the stress as :math:`\symbf{\sigma}=\mathbb{C}:\symbf{\varepsilon}` yields under consideration of
    symmetries in :math:`\mathbb{C}`

    .. math::

        \delta\symbf{\varepsilon}:\mathbb{C}:\symbf{\varepsilon}
            = \delta\mathbf{u}_I \left( \mathbb{C}^{IKLJ} \mathbf{G}_K\otimes\mathbf{G}_L \right) \mathbf{u}_J
            =: \delta\mathbf{u}_I \: \mathbb{K}^{IJ} \:\mathbf{u}_J

    leading to the stiffness matrix as

    .. math::

        \mathbf{K}^{IJ} = \int\!\!\!\int_\mathcal{A} \mathbb{K}^{IJ} \: d\mathcal{A}

    **Implementation notes**

    Writing the kinematic equation in cartesian component form, we obtain

    .. math::

        \left\{\begin{array}{c}
            \varepsilon_{xx} \\
            \varepsilon_{yy} \\
            \gamma_{xy}
        \end{array}\right\}
        = \underbrace{
          \left[\begin{array}{cc}
                G^I_x & 0 \\
                0 & G^I_y  \\
                G^I_y & G^I_x
          \end{array}\right]
            }_{=: [\mathbf{B}^I]}
        \left\{\begin{array}{c}
            u_{I,x} \\ u_{I,y}
        \end{array}\right\}

    leading to

    .. math::

        [ \mathbb{K}^{IJ} ] = [\mathbf{B}^I]^t [\mathbb{C}] [\mathbf{B}^J]


.. dropdown:: Finite deformation theory

    The element uses an  iso-parametric mapping for both the reference and the spatial configuration.
    In the reference configuration,

    .. math::
        :label: eq-triangle1

        \mathbf{X}=\mathbf{X}_I \xi^I \quad I=0,1,2\equiv u,s,t \quad \xi^I \ge 0 ~\&~ \xi^0+\xi^1+\xi^2 = 1

    where :math:`\xi^I` are standard triangle coordinates.
    Alternatively, the mapping may be expressed using just two coordinates, :math:`\xi^1\equiv\xi^s` and
    :math:`\xi^2\equiv\xi^t`, as

    .. math::
        :label: eq-triangle2

        \mathbf{X}=\mathbf{X}_0+\mathbf{G}_i \xi^i \quad i=1,2 \quad \xi^i \ge 0 ~\&~ \xi^1+\xi^2 \leq 1

    The covariant base vectors follow from :eq:`eq-triangle1` and :eq:`eq-triangle2` as

    .. math::
        :label: eq-triangle3

        \mathbf{G}_i = \mathbf{X}_{,i} = \mathbf{X}_i - \mathbf{X}_0
        \quad i=1,2

    We can express differential :math:`d\xi^i` in terms of :math:`d\mathbf{X}` by means of the dual base as

    .. math::
        d\mathbf{X}=\mathbf{G}_i d \xi^i \rightarrow d \xi^i=\mathbf{G}^i \cdot d \mathbf{X} \quad w/ ~~ \mathbf{G}^i \mathbf{G}_j=\delta_j^i

    Similar to the above, we map the spatial configuration as

    .. math::
        \mathbf{x} = \mathbf{x}_I \xi^I = \mathbf{x}_0+ \mathbf{g}_i \xi^i
        \quad I=0,1,2 \quad i=1,2

    Using :eq:`eq-triangle3` we obtain the deformation gradient as

    .. math::
        d\mathbf{x}=\mathbf{g}_i d \xi^i=\left(\mathbf{g}_i \otimes \mathbf{G}^i\right) \cdot d \mathbf{X}=\mathbf{F} \cdot d\mathbf{\bar{x}} \Longrightarrow \mathbf{F}=\mathbf{g}_i \otimes \mathbf{G}^i

    Various strain tensors can be computed from it.  We use the right Cauchy-Green deformation tensor,

    .. math::
        \mathbf{C} = \mathbf{F}^t \mathbf{F}=g_{i j} \mathbf{G}^i \otimes \mathbf{G}^j

    or the Green-Lagrange strain tensor

    .. math::
        \mathbf{E}=\frac{1}{2}\left(\mathbf{F}^t \mathbf{F}-\mathbf{1}\right)=\frac{1}{2}\left(g_{i j}-G_{i j}\right) \mathbf{G}^i \otimes \mathbf{G}^j

    Variation/linearization as needed for the *principle of virtual displacements*  follows as

    .. math::
        \delta \mathbf{F}=\delta\mathbf{g}_i \otimes \mathbf{G}^i
        =\left(\delta \mathbf{x}_i-\delta \mathbf{x}_0\right) \otimes \mathbf{G}^i
        =\delta \mathbf{x}_I \otimes \mathbf{G}^I
        \quad \underline{\text { note }}: \quad \begin{array}{l}
                                                    i=1,2 \\
                                                    I=0,1,2
                                                \end{array}

    where :math:`\mathbf{G}^0=-\mathbf{G}^1-\mathbf{G}^2` and

    .. math::

        \begin{align}
        \delta \mathbf{E}
        &= \frac{1}{2} \left( \delta \mathbf{F}^t \mathbf{F} + \mathbf{F}^t \delta \mathbf{F} \right) \\
        &= \frac{1}{2} \left( \mathbf{G}^I \otimes \delta \mathbf{x}_I \mathbf{g}_i \otimes \mathbf{G}^i
            + \mathbf{G}^i \otimes \mathbf{g}_i \delta \mathbf{x}_I \otimes \mathbf{G}^I \right) \\
        &= \left( \delta \mathbf{x}_I \mathbf{g}_i \right)
             \frac{1}{2} \left( \mathbf{G}^I \otimes \mathbf{G}^i + \mathbf{G}^i \otimes \mathbf{G}^I \right)
            \qquad\qquad I=0,1,2 \quad i=1,2 \\
        &= \left( \delta \mathbf{x}_I \mathbf{x}_i \right)
             \frac{1}{2} \left( \mathbf{G}^I \otimes \mathbf{G}^i + \mathbf{G}^i \otimes \mathbf{G}^I \right)
        - \left( \delta \mathbf{x}_I \mathbf{x}_0 \right)
             \frac{1}{2} \left( \mathbf{G}^I \otimes (\mathbf{G}^1 + \mathbf{G}^2)
            +  (\mathbf{G}^1 + \mathbf{G}^2) \otimes \mathbf{G}^I \right) \\
        &= \left( \delta \mathbf{x}_I \mathbf{x}_J \right)
             \frac{1}{2} \left( \mathbf{G}^I \otimes \mathbf{G}^J + \mathbf{G}^J \otimes \mathbf{G}^I \right)
        \qquad\qquad I,J=0,1,2
        \end{align}

    **Application of virtual work**

    The incremental form of the *principle of virtual displacements* (PVD) takes the form

    .. math::

        \int\!\!\!\int\limits_{\mathcal{A}}\left(
                                \delta\mathbf{E}:\mathbb{C}:d\mathbf{E}
                                + d\delta \mathbf{E}:\mathbf{S}
                                \right) d\mathcal{A}
        = \int\!\!\!\int\limits_{\mathcal{A}} \delta\mathbf{u}\cdot\mathbf{b} \: d\mathcal{A}
        + \int\limits_{\partial\mathcal{A}_\sigma} \delta\mathbf{u}\cdot\bar{\mathbf{T}} \: d\mathcal{s}

    where :math:`\mathbf{S}` is the :math:`2^{nd}` Piola-Kirchhoff stress tensor,
    :math:`\mathbb{C}=\partial\mathbf{S}/\partial\mathbf{E}` is the material stiffness tensor,
    and

    .. math::

        \begin{align}
        d\delta \mathbf{E}
        &= \left( \delta \mathbf{x}_I \:d\mathbf{x}_J \right)
             \frac{1}{2} \left( \mathbf{G}^I \otimes \mathbf{G}^J + \mathbf{G}^J \otimes \mathbf{G}^I \right)
        \qquad I,J=0,1,2
        \end{align}

    The left hand term in the PVD yields

    .. math::

        \delta\mathbf{E}:\mathbb{C}:d\mathbf{E}
            = \delta\mathbf{x}_I \left( \mathbb{C}^{IKLJ} \mathbf{x}_K\otimes\mathbf{x}_L \right) d\mathbf{x}_J
            =: \delta\mathbf{x}_I \: \mathbb{K}^{IJ}_{m} \:d\mathbf{x}_J

    and

    .. math::

        d\delta\mathbf{E}:\mathbf{S}
            = \delta\mathbf{x}_I \: S^{IJ} \mathbf{1} \: d\mathbf{x}_J

    leading to the tangent stiffness matrix as

    .. math::

        \mathbf{K}^{IJ} = \int\!\!\!\int_\mathcal{A} \left( \mathbb{K}^{IJ}_{m} + S^{IJ} \mathbf{1} \right) d\mathcal{A}

    **Implementation notes**

    .. math::

        \delta\mathbf{F}
        =\left(\mathbf{G}^J \otimes \delta\mathbf{x}_J\right)\left(\mathbf{g}_i \otimes \mathbf{G}^i\right)
        =\delta\mathbf{x}_J \cdot  \mathbf{g}_i \otimes \mathbf{G}^J \otimes \mathbf{G}^i
        =:\delta\mathbf{x}_J \cdot{\mathbf{B^{J}}^t}


Parent class
---------------
* :doc:`Element_class`

See also
------------
* :doc:`Triangle6_class`

Class doc
-------------

.. automodule:: femedu.elements.linear.Triangle
  :members:
