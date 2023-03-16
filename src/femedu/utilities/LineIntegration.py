from .Integration import Integration

class LineIntegration(Integration):
    """
    Provides integration points and weights for 1D integration
    """

    def __init__(self, order=2):
        super(LineIntegration, self).__init__(order=order, dimension=1)

        xi, w = self.gauss1D( nGP=(order + 2)//2 )
        self.xi = xi
        self.w  =  w
