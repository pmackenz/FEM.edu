Triangle class for Diffusion
==============================

Coordinate-free formulation for a bi-linear diffusion triangle.

Coordinate mapping:

.. math::

    \mathbf{X}=\mathbf{X}_0+\mathbf{G}_i \xi^i \quad i=1,2,3 \quad \xi^i \geqslant 0 ~\&~ \xi^1+\xi^2+\xi^3 \leq 1

Definition of covariant and dual base vectors:

.. math::

    d\mathbf{X}=\mathbf{G}_i d \xi^i \rightarrow d \xi^i=\mathbf{G}^i \cdot d \mathbf{X} \quad w/ ~~ \mathbf{G}^i \mathbf{G}_j=\delta_j^i

The dual base for the 0-direction of the triangle coordinates:

.. math::

    \mathbf{G}^u = -\mathbf{G}^s - \mathbf{G}^t

Leads to the following expression for the gradient of the scalar potential function:

.. math::

    \nabla \phi = \phi_0 \, \mathbf{G}^u + \phi_1 \, \mathbf{G}^s + \phi_2 \, \mathbf{G}^t

Nodal forces

.. math::

   \mathbf{F}_I = A t \lambda \, \mathbf{G}^I\cdot\nabla\phi
   \qquad\text{with}\quad I=u,s,t

Nodal tangent stiffness

.. math::

    \left({K}_{t}\right)^{IJ} = A t \lambda \, \mathbf{G}^I\cdot\mathbf{G}^J
   \qquad\text{with}\quad I,J=u,s,t

.. math::

    [\mathbf{K}_{t}] =
    \left[
    \begin{array}{ccc}
        \left({K}_{t}\right)^{uu} & \left({K}_{t}\right)^{us} & \left({K}_{t}\right)^{ut} \\
        \left({K}_{t}\right)^{su} & \left({K}_{t}\right)^{ss} & \left({K}_{t}\right)^{st} \\
        \left({K}_{t}\right)^{tu} & \left({K}_{t}\right)^{ts} & \left({K}_{t}\right)^{tt}
    \end{array}
    \right] =
    A t \lambda  \left[
    \begin{array}{ccc}
        \mathbf{G}^u\cdot\mathbf{G}^u & \mathbf{G}^u\cdot\mathbf{G}^s & \mathbf{G}^u\cdot\mathbf{G}^t \\
        \mathbf{G}^s\cdot\mathbf{G}^u & \mathbf{G}^s\cdot\mathbf{G}^s & \mathbf{G}^s\cdot\mathbf{G}^t \\
        \mathbf{G}^t\cdot\mathbf{G}^u & \mathbf{G}^t\cdot\mathbf{G}^s & \mathbf{G}^t\cdot\mathbf{G}^t
    \end{array}
    \right]

Parent class
---------------
* :doc:`../Element_class`

Class doc
-------------

.. automodule:: femedu.elements.diffusion
  :members:
