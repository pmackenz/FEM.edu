import numpy as np
import pandas as pd


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

    def __init__(self, **kwargs):
        self.active = False
        self.data = {}
        self.nodes = []
        self.elements = []

        if 'variables' in kwargs:
            for var in kwargs['variables']:
                self.data[var] = []

        if 'nodes' in kwargs:
            self.nodes = kwargs['nodes']
            for node in self.nodes:
                node.setRecorder(Recorder(variable=self.data.keys()))

        if 'elements' in kwargs:
            self.elements = kwargs['elements']
            for elem in self.elements:
                elem.setRecorder(Recorder(variable=self.data.keys()))


    def fetchRecord(self, keys=None):
        """
        Request recorded time history data for the listed keys.
        If a single key is given as a string, a single `np.array()` is returned.
        If a list of keys is given, a `list of np.array()` is returned.

        :returns: time history data for the listed keys.
        """
        if not keys:
            return self.data

        if isinstance(keys, str):
            if keys in self.data:
                return self.data[keys]
            else:
                return []
        elif isinstance(keys, list) or isinstance(keys, tuple):
            ans = []
            for key in keys:
                ans.append(self.data[key])
            else:
                ans.append([])
            return ans
        else:
            return []

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




