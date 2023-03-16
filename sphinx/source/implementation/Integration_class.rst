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

