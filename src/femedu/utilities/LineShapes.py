import sys
import numpy as np

from .ShapeFunctions import ShapeFunctions

class LineShapes(ShapeFunctions):
    """
    Shape functions or their n-th derivative
    for one-dimensional domains :math:`[0,+1]`.

    Derivatives are with respect to the normalized coordinate and need to be scaled accordingly.

    """

    def __init__(self):
        super(LineShapes, self).__init__(ShapeFunctions.BEAMS)

    def shape(self, order, xi, n=0, Le=1.):
        '''

        :param order: polynomial order of the shape functions. Options: 0,1,2,3,4,5
        :param xi:    normalized location on interval [0,1]
        :param n:     return n-th derivative of the requested shape function
        :param Le:    element length; needed to correct for rotational dof
        :return: PHI: a list of shape function values
        '''

        if (order == 0):
        #
        #   constant shape function
        #
            if (n == 0):
                PHI = np.array([1.0])
            else:
                PHI = np.array([0.0])

        elif (order == 1):
        #
        #   linear shape functions
        #
            if (n == 0):
                PHI = np.array([ 1-xi, xi])
            elif (n == 1):
                PHI = np.array([ -1.0, 1.0 ]) / Le
            else:
                PHI = np.array([0.0, 0.0])

        elif (order == 2):
        #
        #   quadratic shape functions
        #
            if (n == 0):
                PHI = np.array([ 1.-3.*xi+2.*xi*xi, -xi+2.*xi*xi, 4*xi-4.*xi*xi])
            elif (n == 1):
                PHI = np.array([ -3.+4.*xi, -1.+4.*xi, 4.-8.*xi]) / Le
            elif (n == 2):
                PHI = np.array([ 4., 4., -8.]) / (Le*Le)
            else:
                PHI = np.array([0.0, 0.0, 0.0])

        elif (order == 3):
        #
        #   cubic shape functions
        #
            if (n == 0):
                PHI = np.array([ (1 - 3 * xi**2 + 2 * xi**3) ,
                        Le * (xi - 2 * xi**2 + xi**3) ,
                        (3 * xi**2 - 2 * xi**3) ,
                        Le * (-xi**2 + xi**3) ])
            elif (n == 1):
                PHI = np.array([ (- 6. * xi + 6 * xi**2) ,
                        Le * (1. - 4. * xi + 3. * xi**2) ,
                        (6 * xi - 6 * xi**2) ,
                        Le * (-2. * xi + 3. * xi**2) ]) / Le
            elif (n == 2):
                PHI = np.array([ (- 6. + 12. * xi) , Le * (- 4. + 6. * xi) ,
                        (6 - 12 * xi) , Le * (-2. + 6. * xi) ]) / (Le*Le)
            elif (n == 3):
                PHI = np.array([ 12. , Le * 6. , - 12. , Le *  6. ]) / (Le*Le*Le)
            else:
                PHI = np.array([0.0, 0.0, 0.0, 0.0])

        elif (order == 4):
        #
        #   quartic shape functions
        #
            if (n == 0):
                PHI = np.array([ ( 1 - 11*xi**2 + 18*xi**3 - 8*xi**4 ) ,
                        Le * ( xi - 4*xi**2 + 5*xi**3 - 2*xi**4 ) ,
                        ( -5*xi**2 + 14*xi**3 - 8*xi**4 ) ,
                        Le * ( xi**2 - 3*xi**3 + 2*xi**4 ) ,
                        ( 16*xi**2 - 32*xi**3 + 16*xi**4 ) ])
            elif (n == 1):
                PHI = np.array([ ( - 22*xi + 54*xi**2 - 32*xi**3 ) ,
                        Le * ( 1 - 8*xi + 15*xi**2 - 8*xi**3 ) ,
                        ( -10*xi + 42*xi**2 - 32*xi**3 ) ,
                        Le * ( 2*xi - 9*xi**2 + 8*xi**3 ) ,
                        ( 32*xi - 96*xi**2 + 64*xi**3 ) ]) / (Le)
            elif (n == 2):
                PHI = np.array([ ( - 22 + 108*xi - 96*xi**2 ) ,
                        Le * ( - 8 + 30*xi - 24*xi**2 ) ,
                        ( -10 + 84*xi - 96*xi**2 ) ,
                        Le * ( 2 - 18*xi + 24*xi**2 ) ,
                        ( 32 - 192*xi + 192*xi**2 ) ]) / (Le*Le)
            elif (n == 3):
                PHI = np.array([ ( 108 - 192*xi ) ,
                        Le * ( 30 - 48*xi ) ,
                        ( 84 - 192*xi ) ,
                        Le * ( -18 + 48*xi ) ,
                        ( -192 + 384*xi ) ]) / (Le*Le*Le)
            elif (n == 4):
                PHI = np.array([ ( -192 ) ,
                        Le * ( -48 ) ,
                        ( -192 ) ,
                        Le * ( 48 ) ,
                        ( 384 ) ]) / (Le*Le*Le*Le)
            else:
                PHI = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

        elif (order == 5):
        #
        #   quintic shape functions
        #
            if (n == 0):
                PHI = np.array([ ( 1 - 23 * xi**2 + 66 * xi**3 - 68 * xi**4 + 24 * xi**5 ),
                        Le * ( xi - 6 * xi**2 + 13 * xi**3 - 12 * xi**4 + 4 * xi**5 ),
                        ( 7 * xi**2 - 34 * xi**3 + 52 * xi**4 - 24 * xi**5 ),
                        Le * ( -xi**2 + 5 * xi**3 - 8 * xi**4 + 4 * xi**5 ),
                        ( 16 * xi**2 - 32 * xi**3 + 16 * xi**4 ),
                        Le * ( -8 * xi**2 + 32 * xi**3 - 40 * xi**4 + 16 * xi**5 ) ])
            elif (n == 1):
                PHI = np.array([ -46 * xi + 198 * xi**2 - 272 * xi**3 + 120 * xi**4,
                   Le (1 - 12 * xi + 39 * xi**2 - 48 * xi**3 + 20 * xi**4),
                    14 * xi - 102 * xi**2 + 208 * xi**3 - 120 * xi**4,
                   Le (-2 * xi + 15 * xi**2 - 32 * xi**3 + 20 * xi**4),
                   32 * xi - 96 * xi**2 + 64 * xi**3,
                   Le (-16 * xi + 96 * xi**2 - 160 * xi**3 + 80 * xi**4) ]) / (Le)
            elif (n == 2):
                PHI = np.array([ -46 + 396 * xi - 816 * xi**2 + 480 * xi**3,
                   Le (-12 + 78 * xi - 144 * xi**2 + 80 * xi**3),
                   14 - 204 * xi + 624 * xi**2 - 480 * xi**3,
                   Le (-2 + 30 * xi - 96 * xi**2 + 80 * xi**3),
                   32 - 192 * xi + 192 * xi**2,
                   Le (-16 + 192 * xi - 480 * xi**2 + 320 * xi**3) ]) / (Le*Le)
            elif (n == 3):
                PHI = np.array([ 396 - 1632 * xi + 1440 * xi**2,
                   Le (78 - 288 * xi + 240 * xi**2),
                   -204 + 1248 * xi - 1440 * xi**2,
                   Le (30 - 192 * xi + 240 * xi**2),
                   -192 + 384 * xi,
                   Le (192 - 960 * xi + 960 * xi**2) ]) / (Le*Le*Le)
            elif (n == 4):
                PHI = np.array([ -1632 + 2880 * xi,
                   Le (-288 + 480 * xi),
                   1248 - 2880 * xi,
                   Le (-192 + 480 * xi),
                   384,
                   Le (-960 + 1920 * xi) ]) / (Le*Le*Le*Le)
            elif (n == 5):
                PHI = np.array([2880, 480*Le, -2880, 480*Le, 0, 1920*Le ]) / (Le*Le*Le*Le*Le)
            else:
                PHI = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        else:
            msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
            raise NotImplementedError(msg)

        return PHI

