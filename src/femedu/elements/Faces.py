import numpy as np


class Faces():

    def __init__(self, id, *nds, **kwargs):
        self.id = id
        self.nodes = nds
        self.num_nodes = len(nds)

        if self.num_nodes > 0:
            self.dim = nds[0].getPos().size
        else:
            self.dim = 0

    def __str__(self):
        s = "{}_{}".format(self.__class__.__name__, self.id)
        return s
