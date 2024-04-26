.. _recorder_class:

Recorder class
================

A `Recorder` object takes care of recording time-history data during load-stepping, displacement-control,
arc-length, creep, or dynamic analyses.

That data can be exported as text, CSV, JSON, or Excel files for further processing outside of |PackageName|.

.. dropdown::  Abstract Recorder class

    .. automodule:: femedu.recorder.Recorder
      :members:

.. dropdown::  Derived classes

    .. toctree::
        :maxdepth: 1

        ModelRecorder_class.rst
        NodeRecorder_class.rst
        ElementRecorder_class.rst
        MaterialRecorder_class.rst


.. dropdown::  Helper classes

    .. toctree::
        :maxdepth: 1

        Record_class.rst


