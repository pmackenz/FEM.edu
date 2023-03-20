from .Faces import *

class Face3D(Faces):

    def __init__(self, id, *nds, **kwargs):
        super(Face3D, self).__init__(id, *nds, **kwargs)
        if self.dim < 3:
            msg = "{} requires 3D points - {} dimensions given".format(self.__class__.__name__, self.dim)
            raise TypeError(msg)

