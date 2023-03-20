from .Faces import *

class Face2D(Faces):

    def __init__(self, id, *nds, **kwargs):
        super(Face2D, self).__init__(id, *nds, **kwargs)
        if self.dim < 2:
            msg = "{} requires 2D or 3D points - {} dimensions given".format(self.__class__.__name__, self.dim)
            raise TypeError(msg)
