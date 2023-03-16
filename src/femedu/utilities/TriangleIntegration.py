from .Integration import Integration

class TriangleIntegration(Integration):
    """
    Provides integration points and weights for a triangular domain (2d)
    """

    def __init__(self, order=2):
        super(TriangleIntegration, self).__init__(order=order, dimension=1)

        xi, w = self.dunavant( p=order )
        self.xi = xi
        self.w  =  w
