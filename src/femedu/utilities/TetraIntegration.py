import numpy as np
from .Integration import Integration

class TetraIntegration(Integration):
    """
    Provides integration points and weights for tetrahedral domains
    """

    def __init__(self, order=2):
        super(TetraIntegration, self).__init__(order=order, dimension=1)
        self.init_points_and_weights()

    def init_points_and_weights(self):

        if self.order <= 1:
            xi = [ ( 0.25, 0.25, 0.25, 0.25 ) ]
            w  = [ 1/6. ]
        elif self.order == 2:
            a = ( 5. + 3.*np.sqrt( 5. )) / 20.
            b = ( 5. - np.sqrt( 5. )) / 20.
            wi = 1./24.
            xi = [ (a,b,b,b), (b,a,b,b), (b,b,a,b), (b,b,b,a) ]
            w  = [ wi, wi, wi, wi]
        elif self.order == 3:
            a = 0.500
            b = 1./6.
            wi = 9./120.
            xi = [ ( 0.25, 0.25, 0.25, 0.25 ), (a,b,b,b), (b,a,b,b), (b,b,a,b), (b,b,b,a) ]
            w  = [ -2./15., wi, wi, wi, wi ]
        else:
            a = 11./14.
            b = 1./14.
            wi = 343./45000.
            xi = [ ( 0.25, 0.25, 0.25, 0.25 ), (a,b,b,b), (b,a,b,b), (b,b,a,b), (b,b,b,a) ]
            w  = [ -74./5625., wi, wi, wi, wi ]
            a = 0.25 * ( 1. + np.sqrt( 5./14. ))
            b = 0.25 * ( 1. - np.sqrt( 5./14. ))
            wi = 56. / 2250.
            xi += [ (a,a,b,b), (a,b,a,b), (a,b,b,a), (b,a,a,b), (b,a,b,a), (b,b,a,a) ]
            w  += [ wi, wi, wi, wi, wi, wi ]

        self.xi = xi
        self.w  =  w
