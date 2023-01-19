from Solver import *

class LinearSolver(Solver):
    """
    A linear system solver

    This solver implies :math:`\{{\\bf P}\} = [{\\bf K}] \{{\\bf u}\}`
    """

    def __init__(self):
        super().__init__()

