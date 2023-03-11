from .Mesher import *

class PatchMesher(Mesher):
    """
    Mesher for plate and shell elements on a quadrilateral domain,
    including in-plane and out-of plane loaded plates.
    """

    def __init__(self):
        super(PatchMesher, self).__init__()

