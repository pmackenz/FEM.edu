from .Integration import Integration

class BrickIntegration(Integration):
    """
    Provides integration points and weights for brick-like domains (3d)
    """

    def __init__(self, order=2):
        super(BrickIntegration, self).__init__(order=order, dimension=1)
        self.init_points_and_weights(nGP=(order + 2)//2)

    def init_points_and_weights(self, nGP):
        xi, w = self.gauss1D( nGP )
        self.xi = [ xi[i]*xi[j]*xi[k] for k in range(nGP) for j in range(nGP) for i in range(nGP) ]
        self.w  = [  w[i]* w[j]* w[k] for k in range(nGP) for j in range(nGP) for i in range(nGP) ]
