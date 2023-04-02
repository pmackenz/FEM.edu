
class Constraint():
    """
    The head parent for all constrint implementations
    """

    COUNT = 0

    def __init__(self, nodes_list, params, **kwargs):
        self.id = self.COUNT  # creating a unique ID
        self.COUNT += 1       # ensuring the next ID will be unique as well
        self.start = -1       # index in system matrix

        self.numnds = len(nodes_list)
        self.nodes  = nodes_list
        self.params = params

        self.neqns  = 0

    def __str__(self):
        s = "{}: ".format(self.__class__.__name__)
        for node in self.nodes:
            s += "{} ".format(node.getID())
        return s

    def __repr__(self):
        return str(self)

    def countConditions(self):
        """
        :returns: the count of constraint equations in this Constraint object
        """
        return self.neqns

    def on_converged(self):
        pass

    def revert(self):
        pass

