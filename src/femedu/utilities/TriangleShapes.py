import sys
import numpy as np

from .ShapeFunctions import ShapeFunctions

class TriangleShapes(ShapeFunctions):
    """
    Shape functions or their `*n-th` mixed derivative
    for the unit-triangle (two-dimensional domain :math:`(s,t)\in[0,+1]\\times[0,+1]` with :math:`s+t\le +1`).

    Derivatives are with respect to the normalized coordinate and need to be scaled accordingly.

    """

    def __init__(self):
        super(TriangleShapes, self).__init__(ShapeFunctions.TRIANGLE)

    def shape(self, order, s, t, n=(0,0)):
        '''

        :param order: polynomial order of the shape functions. Options: 0,1,2,
        :param s:    first coordinate of point of interest on interval [-1,1]
        :type s:  float or np.array
        :param t:   second coordinate location on interval [-1,1]
        :type t:  float or np.array
        :param n:     return n-th derivative of the requested shape function.
        :type n:  tuple or list

            Use :code:`n=(1,0)` for first derivative with respect to the first variable.
            Use :code:`n=(2,1)` for second derivative with respect to the first variable and
            first with respect to the second variable, i.e., a mixed derivative.

        :return: PHI: a list of shape function values
        '''
        ns, nt = n

        u = 1. - s - t

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
            if (ns == 0 and nt == 0):
                PHI = np.array([ u, s, t ])
            elif (ns == 1 and nt == 0):
                PHI = np.array([ -1., 1., 0. ])
            elif (ns == 0 and nt == 1):
                PHI = np.array([ -1., 0., 1. ])
            else:
                PHI = np.array([0.0, 0.0, 0.0])

        elif (order == 2):
        #
        #   quadratic shape functions
        #
            if (ns == 0 and nt == 0):
                PHI = np.array([ u*(2*u-1), s*(2*s-1), t*(2*t-1), 4*u*s, 4*s*t, 4*t*u ])
            elif (ns == 1 and nt == 0):
                PHI = np.array([ (4*u+1), (4*s-1), 0.0, 4*(u-s), 4*t, -4*t ])
            elif (ns == 0 and nt == 1):
                PHI = np.array([ (4*u+1), 0.0, (4*t-1), -4*s, 4*s, 4*(u-t) ])
            elif (ns == 2 and nt == 0):
                PHI = np.array([ -4.0, 4.0, 0.0,-8.0, 0.0, 0.0 ])
            elif (ns == 1 and nt == 1):
                PHI = np.array([ -4.0, 0.0, 0.0, -4.0, 4.0, -4.0 ])
            elif (ns == 0 and nt == 2):
                PHI = np.array([ -4.0, 0.0, 4.0, 0.0, 0.0, -8.0 ])
            else:
                PHI = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        else:
            msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
            raise NotImplementedError(msg)

        return PHI
