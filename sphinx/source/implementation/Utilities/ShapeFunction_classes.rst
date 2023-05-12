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

Integrating a bi-quadratic function :code:`g(s,t)*shape_function(s,t)`
over a quad-shaped domain, assuming a function
:code:`Jacobian(s,t)` is defined and returns the jacobian for the current location.

.. code:: python

    integrator    = QuadIntegration(order=2)
    interpolation = QuadShapes()

    # initialization
    F = np.zeros((9,))  # 9-node quad

    # integration loop
    xis, wis = integrator.parameters()
    for xi, wi in zip(xis, wis):
        # computing the array of nodal shape functions
        shape = interpolation.shape(  # requesting spage function array
                    order=2,          # polynomial order per direction: quadratic
                    s=xi[0], t=xi[1], # local coordinates for current position
                    n=(0,0),          # n-th derivative with respect to (s,t)
                    serendipity=False # if serendipity: "8 nodes" else: "9 nodes"
                    )
    `   # adding to the integral
        F += g(xi[0], xi[1]) * shape * Jacobian(xi[0], xi[1]) * wi

    print(f"Int_V g(s,t) Phi(s,t) dV = {F}")



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

