from .Mesher import *

class TriPatchMesher(Mesher):
    """
    Mesher for plate and shell elements on a triangular domain,
    including in-plane and out-of plane loaded plates.
    """

    def __init__(self):
        super(TriPatchMesher, self).__init__()

