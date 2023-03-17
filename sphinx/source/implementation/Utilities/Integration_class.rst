Integration classes
==========================

This class and its sub-classes provide support for numeric integration over
common element domains.

Implement most of the methods presented in `<http://me.rice.edu/\~akin/Elsevier/Chap\_10.pdf>`_.

**Usage**:

* Instantiate the Integration object most suitable for your element domain. Provide :py:data:`order`
  as the highest polynomial order of the function you want to integrate exactly.

* call :py:meth:`points()` or :py:meth:`weights()` or :py:meth:`parameters()`
  to obatin a list of integration point coordinates or weights or both, respectively.

Integrating a 4th-order function f(s,t) over a triangular domain

.. code:: python

    integrator = TrangleIntegration(order=4)
    F = 0.0
    for xi, wi in integrator.parameters():
        F += f(xi[0], xi[1]) * J(xi[0], xi[1]) * wi

    print(f"Int_A f(s,t) dA = {F}")


Integrating a tri-quadratic function g(s,t,u) over a brick-shaped domain

.. code:: python

    integrator = BrickIntegration(order=2)
    F = 0.0
    for xi, wi in integrator.parameters():
        F += g(xi[0], xi[1], xi[2]) * J(xi[0], xi[1], xi[2]) * wi

    print(f"Int_V g(s,t,u) dV = {F}")


Abstract Integration Class
--------------------------------

.. automodule:: femedu.utilities.Integration
  :members:


Derived Classes
-------------------

.. automodule:: femedu.utilities.LineIntegration
  :members:


.. automodule:: femedu.utilities.TriangleIntegration
  :members:


.. automodule:: femedu.utilities.QuadIntegration
  :members:


.. automodule:: femedu.utilities.TetraIntegration
  :members:


.. automodule:: femedu.utilities.BrickIntegration
  :members:

