import sys
import numpy as np

from .ShapeFunctions import ShapeFunctions
from .LineShapes import LineShapes

class QuadShapes(ShapeFunctions):
    """
    Shape functions or their :code:`*n-th` mixed derivative
    for the bi-unit-square (two-dimensional domain :math:`(s,t)\in[-1,+1]\\times[-1,+1]`.

    Derivatives are with respect to the normalized coordinate and need to be scaled accordingly.

    """

    def __init__(self):
        super(QuadShapes, self).__init__(ShapeFunctions.QUADS)

    def shape(self, order, s, t, n=(0,0), serendipity=False):
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
            first with respect to the second variable, i.e., a mixed derivative..

        :param serendipity: set to **True** for 8-node serendipity or **False** for bi-quadratic (9-node) quads.
        :type serendipity: bool

        :return: PHI: a list of shape function values
        '''
        ns, nt = n

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
                PHI = np.array([ 0.25 * (1-s) * (1-t),
                                 0.25 * (1+s) * (1-t),
                                 0.25 * (1+s) * (1+t),
                                 0.25 * (1-s) * (1+t) ])
            elif (ns == 1 and nt == 0):
                PHI = np.array([ -0.25 * (1-t),
                                  0.25 * (1-t),
                                  0.25 * (1+t),
                                 -0.25 * (1+t) ])
            elif (ns == 0 and nt == 1):
                PHI = np.array([ -0.25 * (1-s),
                                 -0.25 * (1+s),
                                  0.25 * (1+s),
                                  0.25 * (1-s) ])
            elif (ns == 1 and nt == 1):
                PHI = np.array([  0.25, -0.25, 0.25, -0.25 ])
            else:
                PHI = np.array([0.0, 0.0])

        elif (order == 2 and serendipity):
        #
        #   quadratic shape functions: 8-node quad
        #
            if (ns == 0 and nt == 0):
                PHI = np.array([ -0.25 * (1 - s) * (1 - t) * (1 + s + t),
                                 -0.25 * (1 + s) * (1 - t) * (1 - s + t),
                                 -0.25 * (1 + s) * (1 + t) * (1 - s - t),
                                 -0.25 * (1 - s) * (1 + t) * (1 + s - t),
                                  0.50 * (1 - s) * (1 + s) * (1 - t),
                                  0.50 * (1 + s) * (1 + t) * (1 - t),
                                  0.50 * (1 - s) * (1 + s) * (1 + t),
                                  0.50 * (1 - s) * (1 + t) * (1 - t) ])
            elif (ns == 1 and nt == 0):
                PHI = np.array([
                    -(0.25) * (-1 + t) * (2 * s + t), 
                    -(0.25) * (2 * s - t) * (-1 + t), 
                     0.25 * (1 + t) * (2 * s + t),
                     0.25 * (2 * s - t) * (1 + t),
                     s * (-1 + t),
                    -(0.50) * (-1 + t) * (1 + t),
                    -s * (1 + t),
                     0.50 * (-1 + t) * (1 + t)
                ])
            elif (ns == 0 and nt == 1):
                PHI = np.array([
                    -(0.25) * (-1 + s) * (s + 2 * t),
                    -(0.25) * (1 + s) * (s - 2 * t),
                     0.25 * (1 + s) * (s + 2 * t),
                     0.25 * (-1 + s) * (s - 2 * t),
                     0.50 * (-1 + s) * (1 + s),
                    -(1 + s) * t,
                    -(0.50) * (-1 + s) * (1 + s),
                    (-1 + s) * t
                ])
            elif (ns == 2 and nt == 0):
                PHI = np.array([
                    (1 - t)/2, (1 - t)/2, (1 + t)/2, (1 + t)/2, -1 + t, 0, -1 - t, 0
                ])
            elif (ns == 1 and nt == 1):
                PHI = np.array([
                    0.25 * (1 - 2 * s - 2 * t),
                    0.25 * (-1 - 2 * s + 2 * t),
                    0.25 * (1 + 2 * s + 2 * t),
                    0.25 * (-1 + 2 * s - 2 * t),
                    s,  -t,  -s,  t
                ])
            elif (ns == 0 and nt == 2):
                PHI = np.array([
                    (1 - s)/2., (1 + s)/2., (1 + s)/2., (1 - s)/2.,
                    0., -1. - s, 0., -1. + s
                ])
            elif (ns == 2 and nt == 1):
                PHI = np.array([
                    -(0.50), -(0.50), 0.50, 0.50, 1, 0, -1, 0
                ])
            elif (ns == 1 and nt == 2):
                PHI = np.array([
                    -(0.50), 0.50, 0.50, -(0.50), 0, -1, 0, 1
                ])
            else:
                PHI = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        elif (order == 2):
        #
        #   quadratic shape functions: 9-node quad
        #
            line = LineShapes()
            phi_s = line.shape(order=2, xi=s, n=n[0])
            phi_t = line.shape(order=2, xi=t, n=n[1])

            PHI = phi_s*phi_t  # this is an element-by element multiplication

        elif (order == 3):
        #
        #   cubic shape functions: using tensor multiplication of beam functions.
        #   these require rotations!
        #
            line = LineShapes()
            phi_s = line.shape(order=3, xi=s, n=n[0])
            phi_t = line.shape(order=3, xi=t, n=n[1])

            PHI = phi_s*phi_t  # this is an element-by element multiplication

        else:
            msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
            raise NotImplementedError(msg)

        return PHI
