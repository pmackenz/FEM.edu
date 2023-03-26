import numpy as np


class ShapeFunctions():
    """
    Abstract parent class for shape functions or their
    n-th (mixed) derivative for various interpolation functions.

    """

    UNKNOW    = 0x0000
    BEAMS     = 0x1000
    QUADS     = 0x2000
    TRIANGLE  = 0x4000
    BRICK     = 0x8000

    CONST     = 0x0001
    LINEAR    = 0x0002
    QUADRATIC = 0x0004
    CUBIC     = 0x0008
    OTHER     = 0x0010

    def __init__(self, type):
        self.type = type

    def shape(self, order, *argvs, **kwargs):
        """
        interface function to gets shape functions or their
        n-th derivative for various interpolation functions.

        .. note::

           This is a virtual function that needs to be implemented
           by any subclass

        """
        msg = "** WARNING ** {}.{} not implemented".format(self.__class__.__name__, sys._getframe().f_code.co_name)
        raise NotImplementedError(msg)


def selfTestShapes():
    for order in range(5):
        print(f"order: {order}")
        for derivative in range(order+2):
            res = np.array( [ LineShapes().shape(order,xi,1.0,derivative) for xi in np.linspace(0,1,5) ] )
            print(res.T)
        print()

if __name__ == "__main__":
    selfTestShapes()
