from .Integration import Integration

class TriangleIntegration(Integration):
    """
    Provides integration points and weights for 1D integration
    """

    def __init__(self, order=2):
        super(TriangleIntegration, self).__init__(order=order, dimension=1)

        xi, w = self.gauss1D( nGP=(order + 2)//2 )
        self.xi = xi
        self.w  =  w
