from .Faces import *

class Face3D(Faces):
    """
    Implementation representing one face of a 3D brick or tetrahedral element

    The face may be defined by

    * three nodes: liner interpolation of triangular face
    * four nodes:  linear interpolation of an 8-node brick.
    * six nodes:   quadratic interpolation of triangular face
    * eight nodes: quadratic interpolation of a 20-node brick.
    * nine nodes:  quadratic interpolation of a 27-node brick.

    .. note::

        The node order of 3D faces follows the node order of similar shaped 2D elements.

    """

    def __init__(self, id, *nds, **kwargs):
        super(Face3D, self).__init__(id, *nds, **kwargs)
        
    def initialize(self):
        """
        This method computes area and normals for :py:class:`Face3D`.
        
        It is called by the constructor of the :py:class:`Faces` class.
        """
        if self.dim < 3:
            msg = "{} requires 3D points - {} dimensions given".format(self.__class__.__name__, self.dim)
            raise TypeError(msg)

        super(Face3D, self).initialize()

