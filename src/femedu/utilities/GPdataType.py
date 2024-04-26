class GPdataType():
    """
    Gauss Point data type

    Holds:
    :param base: base vectors of local coordinate system
    :param dual_base: dual base vectors of local coordinate system
    :param state: state dictionary at current point. Holds strain, stress, etc.
    :param J: Jacobian (area or volume ratio) at current point
    :param B: B-matrix (kinematic matrix) at current point
    :param C: C-matrix (material tangent) at current point
    :param material: pointer to material object at current point
    """

    def __init__(self):
        self.X         = None
        self.base      = None
        self.dual_base = None
        self.J         = 1.0
        self.B         = None
        self.state     = {}
        self.C         = 1.0
        self.material  = None

    def __str__(self):
        s  = f'GP @ X={self.base}\n'
        s += f'* base:\n\t{self.base}\n'
        s += f'* dual base:\n\t{self.dual_base}\n'
        s += f'* Jacobian (J):\n\t{self.J}\n'
        s += f'* state:\n\t{self.state}\n'
        s += f'* B:\n\t{self.B}\n'
        s += f'* C:\n\t{self.C}\n'
        s += f'* material:\n\t{repr(self.material)}'
        return s