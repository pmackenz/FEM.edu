import numpy as np
import pandas as pd

from ..domain import Node


class Recorder():
    """A time history recorder

    :param \**kwargs: additional parameters passed through to specialized recorders.
        See below for more details.

    :Keyword Arguments:
        * **variable** (``list`` of ``str`` or ``All``)
          -- Extra stuff
        * **nodes** (``list`` of ``Node`` ptrs or ``All``)
          -- nodes for which those variables shall be recorded. Use ``nodes=All`` to apply to all nodes.
        * **elements** (``list`` of ``Element`` ptrs or ``All``)
          -- elements for which those variables shall be recorded. Use ``elements=All`` to apply to all elements.
    """

    # hard-coded variable assignments
    system_types = (
        'lam',
        'time',
        'stability',
        'arc-length',
    )
    node_types = (
        'ux',
        'uy',
        'uz',
        'rx',
        'ry',
        'rz',
        'vx',
        'vy',
        'vz',
        'omx',
        'omy',
        'omz',
        'T',
        'stress',
        'strain'
    )
    element_types = (
        'stress',
        'strain',
    )

    def __init__(self, **kwargs):
        self.active = False
        self.data = {}
        self.nodes = []
        self.elements = []

        nodevars = []
        elemvars = []

        if 'variables' in kwargs:
            for var in kwargs['variables']:
                var_classified = False
                if var in self.system_types or 'type' in kwargs:
                    self.data[var] = []
                    var_classified = True
                if var in self.node_types:
                    nodevars.append(var)
                    var_classified = True
                if var in self.element_types:
                    elemvars.append(var)
                    var_classified = True
                if not var_classified:
                    # this is an unknown/new variable
                    # maybe from a new element type?
                    # just for safety, add it to all types
                    self.data[var] = []
                    nodevars.append(var)
                    elemvars.append(var)

        if 'nodes' in kwargs:
            self.nodes = kwargs['nodes']
            for node in self.nodes:
                node.setRecorder(Recorder(variables=nodevars, type='node'))

        if 'elements' in kwargs:
            self.elements = kwargs['elements']
            for elem in self.elements:
                elem.setRecorder(Recorder(variables=elemvars, type='element'))

    def fetchRecord(self, keys=None, source=None):
        """
        Request recorded time history data for the listed keys.

        :param keys: If a single key is given as a string, a single `np.array()` is returned.
                     If a list of keys is given, a `list of np.array()` is returned.
        :param source: If a single :py:class:`Node` is given, pick record from that node.
                       If a list of :py:class:`Node` objects is given, match keys and nodes
                       based on order in lists.  **source** and **keys** list must match in shape.
                       If a value should be picked from the model domain instead of a node, enter **None**
                       in the respective slot.
        :returns: time history data for the listed keys.
        """
        if not keys:
            return self.data

        if isinstance(keys, str):
            if source and isinstance(source,Node.Node):
                if source.recorder:
                    return source.recorder.fetchRecord(keys)
                else:
                    return {}
            elif source and (isinstance(source,list) or isinstance(source,tuple)) \
                    and isinstance(source[0],Node.Node):
                if source[0].recorder:
                    return source[0].recorder.fetchRecord(keys)
                else:
                    return {}
            else:
                if keys in self.data:
                    return {keys:self.data[keys]}
                else:
                    return {}

        elif isinstance(keys, list) or isinstance(keys, tuple):
            ans = {}

            if source and (isinstance(source,list) or isinstance(source,tuple)):
                for key, src in zip(keys, source):
                    if src and isinstance(src, Node.Node):
                        if src.recorder:
                            node_data = src.recorder.fetchRecord(key)
                            for field in node_data:
                                ans[field] = node_data[field]
                    else:
                        if key in self.data:
                            ans[key] = self.data[key]
                        else:
                            ans[key] = [0,]

            else:
                for key in keys:
                    if key in self.data:
                        ans[key] = self.data[key]
                    else:
                        ans[key] = [0,]
            return ans
        else:
            return {}

    def addData(self, dta):
        """
        *For internal use only.*

        :param dta: variable code as key and scalar value pairs.
        :type dta: dict
        """
        for var in self.data:
            if var in dta:
                self.data[var].append(dta[var])
            else:
                self.data[var].append(np.nan)
                print(f"Recorder.addData: '{var}' not not in data set: padding with nan")

        for var in dta:
            if var not in self.data:
                print(f"Recorder.addData: '{var}' not initialized by the recorder: ignored")

    def getVariables(self):
        return list(self.data.keys())
    def enable(self):
        """
        Enables the recorder.  You may suspend data collection without loss of previous data by calling ``disable()``

        **Note**: This method should not be called by the user.  Use ``System.startRecorder()`` instead.
        """
        self.active = True

    def disable(self):
        """
        Disables the recorder.  You may resume data collection without loss of previous data by calling ``enable()``

        **Note**: This method should not be called by the user.
        Use ``System.pauseRecorder()`` or  ``System.stopRecorder()`` instead.
        """
        self.active = False

    def isActive(self):
        """
        :returns: **True** if recording is on, False otherwise.
        """
        return self.active

    def reset(self):
        """
        This will reset the data recorder to a *disabled* state and wipe all previously collected data.
        """
        self.active = False
        for var in self.data:
            self.data[var] = []

    def export(self, filename='unknown.txt'):
        """
        :param filename: full path to file where recorded data shall be written to.
            The file type will be determined from the given extension.
            Recognized file types are

            .. list-table::

                * - .txt
                  - tab-separated text file
                * - .csv
                  - comma-separated text file
                * - .hdf
                  - HDF5 file
                * - .json, .jsn
                  - comma-separated text file
                * - .xlsx
                  - MicroSoft Excel file (see the `pandas` package for more detail)
        :type filename: ``str``
        """
        file_parts = filename.split['.']
        if file_parts and len(file_parts)>1:
            suffix = file_parts[-1].lower()
        else:
            suffix = 'txt'   # use text file as a default

        if suffix in ('xlsx',):
            self.export_excel(filename)
        elif suffix in ('csv',):
            self.export_csv(filename)
        elif suffix in ('hdf','hdf5'):
            self.export_hdf(filename)
        elif suffix in ('jsn','json'):
            self.export_json(filename)
        else:
            self.export_text(filename)

    def export_text(self, filename):
        """
        .. warning::

            This method is has not yet been implemented.

        """
        df = pd.DataFrame(self.data)
        df.to_csv(filename, sep='\t', header=True)

    def export_csv(self, filename):
        """
        .. warning::

            This method is has not yet been implemented.

        """
        df = pd.DataFrame(self.data)
        df.to_csv(filename, sep=',', header=True)

    def export_hdf(self, filename):
        """
        .. warning::

            This method is has not yet been implemented.

        """
        df = pd.DataFrame(self.data)
        df.to_hdf(filename)

    def export_json(self, filename):
        """
        .. warning::

            This method is has not yet been implemented.

        """
        df = pd.DataFrame(self.data)
        df.to_json(filename)

    def export_excel(self, filename):
        """
        .. warning::

            This method is has not yet been implemented.

        """
        df = pd.DataFrame(self.data)
        df.to_excel(filename)




