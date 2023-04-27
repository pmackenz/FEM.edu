Diffusion Material classes
==============================

**Diffusion** materials cover a wide variety of problems governed by the Laplace equation.
Stationary problems require the diffusivity constant for the (isotropic) Darcy's equation.

Derived classes provide interfaces for material properties commonly used in the context of one particular
diffusion type problem, translating them to the :py:class:`DiffusionGeneral` class' specific parameters,
and inherit all other methods from :py:class:`DiffusionGeneral`.

Parent class
---------------
* :doc:`Material_class`


.. automodule:: femedu.materials.DiffusionGeneral
  :members:


Derived Classes
-------------------

.. automodule:: femedu.materials.Thermal
  :members:

