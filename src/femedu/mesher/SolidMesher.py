from .Mesher import *

class SolidMesher(Mesher):
    """
    Mesher for 8-node, 20-node, and 27-node brick elements.
    """

    def __init__(self):
        super(SolidMesher, self).__init__()
