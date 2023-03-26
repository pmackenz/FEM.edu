
class Integration():
    """
    Abstract interface definition for all Integration classes.
    """

    def __init__(self, order=2, dimension=1):
        self.order = order
        self.xi = [ 0.0 ]
        self.w  = [ 1.0 ]

    def weights(self):
        """
        :returns: list of weights for all integration points
        """
        return self.w

    def points(self):
        """
        :returns: list of integration point coordinates
        """
        return self.xi

    def parameters(self):
        """
        This is a short form of

        .. code::

            xi = points()
            w  = weights()

        :returns: list of integration point coordinates, list of weights for those points
        """
        return (self.xi, self.w)

    def gauss1D(self, nGP, biunit=False):
        """
        Gauss integration data for :py:data:`nGP` on the domain :math:`[-1,+1]`

        For internal use only.
        """

        if (nGP == 1):
            xi = [ 0.5 ]
            w  = [ 1.0 ]
        elif (nGP == 2):
            xi = [ 0.211324865405187, 0.788675134594813 ]
            w  = [ 0.500000000000000, 0.500000000000000 ]
        elif (nGP == 3):
            xi = [0.112701665379258, 0.50000000000000, 0.88729833462074]
            w  = [0.277777777777778, 0.444444444444444, 0.277777777777778]
        elif (nGP == 4):
            xi = [ 0.069431844202974, 0.33000947820757, 0.66999052179243, 0.93056815579703 ]
            w  = [ 0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727 ]
        else:
            xi = [ 0.046910077030668, 0.23076534494716, 0.50000000000000, 0.76923465505284, 0.95308992296933 ]
            w  = [ 0.118463442528095, 0.239314335249683, 0.284444444444444, 0.239314335249683, 0.118463442528095 ]

        if biunit:
            xi = [ 2*s-1. for s in xi ]
            w  = [ 2*s    for s in  w ]

        return (xi, w)

    def dunavant(self, p=0):
        """
        Triangular domain integration for polynomial of order <= p

        For internal use only.
        """

        if (p <= 1):
            xi = [ (1./3, 1./3., 1./3.) ]
            w  = [ 0.5 ]
        elif (p == 2):
            wi = 1./6.
            xi = [ (1./6., 1./6., 4./6.), (4./6., 1./6., 1./6.), (1./6., 4./6., 1./6.) ]
            w  = [ wi, wi, wi ]
        elif (p == 3):
            wi = 25./96.
            xi = [ (1./3, 1./3., 1./3.), (0.20, 0.20, 0.60), (0.60, 0.20, 0.20), (0.20, 0.60, 0.20) ]
            w  = [ -9./32., wi, wi, wi ]
        else:
            wi = 1./15.
            xi = [ (0.0, 0.0, 1.0), (0.5, 0.0, 0.5),
                   (1.0, 0.0, 0.0), (0.5, 0.5, 0.0),
                   (0.0, 1.0, 0.0), (0.0, 0.5, 0.5),
                   (1./3, 1./3., 1./3.) ]
            w  = [ 0.025, wi, 0.025, wi, 0.025, wi, 0.225]

        return (xi, w)

