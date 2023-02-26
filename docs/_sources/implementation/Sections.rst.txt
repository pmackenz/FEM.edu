Section Material classes
==========================

**Section** materials represent either an entire cross-section of a beam or frame element,
or the material column perpendicular to the middle plane of a plate or shell.
All Section materials are derived from the following base class.
That class may be used directly to represent a linear elastic section without coupling between
axial and flexural deformation.

Parent class
---------------
* :doc:`Material_class`


.. automodule:: femedu.materials.SectionMaterial
  :members:


Derived Classes
-------------------

.. toctree::
    :maxdepth: 1

    ElasticSection.rst
    FiberSection.rst
    PlateSection.rst
