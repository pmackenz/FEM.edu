"""
============================
9-noded quadrilateral
============================
Full bi-quadratic interpolation

         1
      x     y
  x^2   x*y   y^2
    x^2*y x*y^2
      x^2*y^2

"""
from .Element import *

class Quad9(Element):
    """

    """

    def __init__(self, *pts, **kwargs):
        super(Quad9, self).__init__(nodes=pts, **kwargs)
        self.element_type = DrawElement.QUAD
        self.createFaces()
