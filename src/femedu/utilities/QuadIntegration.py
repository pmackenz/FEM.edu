from .Integration import Integration

class QuadIntegration(Integration):
    """
    Provides integration points and weights for quadrilateral domains (2D)
    """

    def __init__(self, order=2):
        super(QuadIntegration, self).__init__(order=order, dimension=1)
        self.init_points_and_weights(nGP=(order + 2)//2)

    def init_points_and_weights(self, nGP):
        xi, w = self.gauss1D( nGP, biunit=True )
        self.xi = [ (xi[i], xi[j]) for j in range(nGP) for i in range(nGP) ]
        self.w  = [  w[i]* w[j] for j in range(nGP) for i in range(nGP) ]

