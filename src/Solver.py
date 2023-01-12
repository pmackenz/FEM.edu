
class Solver():
    """
    Abstract class for any solver implementation.

    This class describes the functions needed by any solver
    """

    def __init__(self):
        """
        Initialize a solver instance with empty elements and nodes lists
        """
        self.elements = []
        self.nodes = []

    def assemble(self):
        """

        :return:
        """
        pass

    def solve(self):
        """

        :return:
        """
        pass

    def initialize(self):
        """

        :return:
        """
        pass

    def reset(self):
        """

        :return:
        """
        pass


