Triangle class
==========================

Coordinate-free formulation for a finite deformation bi-linear triangle.

.. math::

    \mathbf{X}=\mathbf{X}_0+\mathbf{G}_i \xi^i \quad i=1,2,3 \quad \xi^i \geqslant 0 ~\&~ \xi^1+\xi^2+\xi^3 \leq 1

.. math::

    d\mathbf{X}=\mathbf{G}_i d \xi^i \rightarrow d \xi^i=\mathbf{G}^i \cdot d \mathbf{X} \quad w/ ~~ \mathbf{G}^i \mathbf{G}_j=\delta_j^i

.. math::

    \mathbf{x}=\mathbf{x}_0+ \mathbf{g}_i \xi^i

.. math::

    d\mathbf{x}=\mathbf{g}_i d \xi^i=\left(\mathbf{g} \otimes \mathbf{G}^i\right) \cdot d \mathbf{X}=\mathbf{F} \cdot d\mathbf{\bar{x}} \Longrightarrow \mathbf{F}=\mathbf{g}_i \otimes \mathbf{G}^i

.. math::

    \mathbf{F}^t \mathbf{F}=g_{i j} \mathbf{G}^i \otimes \mathbf{G}^j \rightarrow \mathbf{E}=\frac{1}{2}\left(\mathbf{F}^t \mathbf{F}-\mathbf{1}\right)=\frac{1}{2}\left(g_{i j}-G_{i j}\right) \mathbf{G}^i \otimes \mathbf{G}^j

.. math::

    \delta \mathbf{E}=\delta\mathbf{g}_i \otimes \mathbf{G}^i=\left(\delta \mathbf{x}_i-\delta \mathbf{x}_0\right) \otimes \mathbf{G}^i=\delta \mathbf{x}_\alpha \otimes \mathbf{G}^\alpha \quad \underline{\text { note }}: \quad \begin{array}{l}
    i=1,2,3 \\
    \alpha=0,1,2,3
    \end{array}

where :math:`\mathbf{G}^0=-\mathbf{G}^1-\mathbf{G}^2-\mathbf{G}^3` and

.. math::

    \delta\mathbf{F}
    =\left(\mathbf{G}^\alpha \otimes \delta\mathbf{x}_\alpha\right)\left(\mathbf{g}_i \otimes \mathbf{G}^i\right)
    =\delta\mathbf{x}_\alpha \cdot \underbrace{\mathbf{g}_i \otimes \overbrace{\mathbf{G}^\alpha \otimes \mathbf{G}^i}^{\bar{\mathbf{G}}^a}}_{=: {\mathbf{B}^{\alpha}}^t}
    =\delta\mathbf{x}_\alpha \cdot{\mathbf{B^{\alpha}}^t}


Parent class
---------------
* :doc:`Element_class`

Class doc
-------------

.. automodule:: femedu.elements.linear.Triangle
  :members:
