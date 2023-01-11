Plotter
=============

Creates undeformed and deformed plots of the system.


.. list-table:: Plotter class methods
   :widths: 25 25 25 50
   :header-rows: 1

   * - method
     - input
     - returns
     - description
   * - `__init__()`
     - 
     - 
     - constructor. Initialize the plotter object to sensible default settings, as needed.
   * - `setMesh(verts,lines)`
     - list of points, list of line indices
     - 
     - replace `self.vertices` and `self.lines` information.
   * - `setDisplacements(disp)`
     - list of displacement vectors
     - 
     - replace `self.disp` information.
   * - `setValues(vals)`
     - list of line (force) values.
     - 
     - replace `self.values` information.
   * - `displacementPlot(file=None)`
     - a string
     - 
     - creates a plot showing undeformed in black and deformed model in red lines. 
       If `file` is given, save a copy of the plot to a file
       of that name
   * - `valuePlot(deformed=False, file=None)`
     - a string
     - 
     - creates a plot showing the undeformed|deformed system (based on the user input) with
       lines colored based on `values`. Add a colormap/colorbar as legend.
       If `file` is given, save a copy of the plot to a file
       of that name

.. list-table:: Plotter class variables
   :widths: 25 25 50
   :header-rows: 1

   * - name
     - type
     - description
   * - `vertices`
     - List of `np.array([X,Y])`
     - list of coordinate pairs representing points (nodes in the model)
   * - `lines`
     - List of List
     - list of 2-element lists of indices.  The two lists shall contain the indices of the
       start and end point of a line in the `vertices list`, respectively.  
   * - `disp`
     - list of `np.array([u,v])`
     - list of point displacements for deformed plot.  This list must be of identical shape
       as the `vertices` list such that respective entries represent point position and
       displacement, respectively.
   * - `values`
     - `np.array([...])`
     - list containing the force values for each line (element).  This list must be of
       identical shape as the `lines` list.


