Element class
==========================

This class represents no particular element but rather defines common functionality for all elements.

.. note::

    Every element in **FEM.edu** is derived from this class.

The :code:`Element` class inherits methods from :code:`DrawElement`, which provides default plotting functionality
for common element types.

.. note::

    Every element needs to declare its generic type through :code:`self.element_type` in its constructor
    for the default plotting mechanism to work.  The default constructor defines :code:`self.element_type = self.UNKNOWN`,
    which means that this element **will not be plotted**.

    See `Inherited methods`_ below for more detail on the :code:`DrawElement` class.


.. automodule:: elements.Element
  :members:

Inherited methods
---------------------

.. automodule:: elements.DrawElement
  :members:


Derived Classes
-------------------

.. toctree::
    :maxdepth: 1

    Truss_class.rst
    Beam2D_class.rst
    Frame2D_class.rst
    LinearTriangle_class.rst
