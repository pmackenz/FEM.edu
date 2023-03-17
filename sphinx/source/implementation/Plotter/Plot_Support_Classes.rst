Plot Support classes
==========================
All plotters are based on an abstract plotter class that defines the fundamental
plot features.  Actual plotter classes are derived from that class and provide specific implementations
of those features.

In addition, all :doc:`../Elements/Element_class` inherit plot support methods from :code:`DrawElement`.


.. automodule:: femedu.plotter.AbstractPlotter
   :members:

Derived Classes
-------------------

.. toctree::
    :maxdepth: 1

    Plotter_class.rst
    Plotter3D_class.rst
    ElementPlotter_class.rst
    ElementPlotter3D_class.rst

