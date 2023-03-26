Shape Function classes
==========================

This family of classes computes common shape functions in 1D (Lines), 2D (Quads and Triangles), and
3D (Bricks and Tetrahedron elements).

**LineShapes** are defined over a unit-interval :math:`[0,1]`.

**QuadShapes** are defined over a bi-unit-square :math:`[-1,1]\times[-1,1]`.

**TriangleShapes** are defined over a domain (triangular coordinates) :math:`[0,1]\times[0,1]\times[0,1]`.

**BrickShapes** are defined over a tri-unit-cube :math:`[-1,1]\times[-1,1]\times[-1,1]`.

**TetraShapes** are defined over a domain (tetrahedral coordinates) :math:`[0,1]\times[0,1]\times[0,1]\times[0,1]`.

In its final form, they will support up to tri-quadratic shape functions, and up to 5th-order polynomials for lines.

**Usage**:

Integrating a bi-quadratic function g(s,t,u) over a quad-shaped domain

.. code:: python

    integrator = QuadIntegration(order=2)
    F = 0.0
    for xi, wi in integrator.parameters():
        F += g(xi[0], xi[1], xi[2]) * J(xi[0], xi[1], xi[2]) * wi

    print(f"Int_V g(s,t,u) dV = {F}")


Abstract Shape Function Class
--------------------------------

.. automodule:: femedu.utilities.ShapeFunctions
  :members:


Derived Classes
-------------------

.. automodule:: femedu.utilities.LineShapes
  :members:


.. automodule:: femedu.utilities.TriangleShapes
  :members:


.. automodule:: femedu.utilities.QuadShapes
  :members:


.. .. automodule:: femedu.utilities.TetraShapes
..   :members:


.. .. automodule:: femedu.utilities.BrickShapes
..   :members:

