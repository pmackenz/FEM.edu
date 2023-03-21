from .Faces import *

class Face2D(Faces):
    """
    Implementation representing one face of a 2D element

    The face may be defined by two nodes (liner interpolation)
    or three nodes (quadratic interpolation).
    The node order for quadratic elements is nodes at `[-1, 0, +1]`
    """


    def __init__(self, id, *nds, **kwargs):
        super(Face2D, self).__init__(id, *nds, **kwargs)

    def initialize(self):
        """
        This method computes area and normals for :py:class:`Face2D`.

        It is called by the constructor of the :py:class:`Faces` class.
        """
        if self.dim < 2:
            msg = "{} requires 2D or 3D points - {} dimensions given".format(self.__class__.__name__, self.dim)
            raise TypeError(msg)

        rot = np.array([[0.,-1.],[1.,0.]])

        if len(self.nodes) == 2:
            # this includes the weight factor for single point integration over [-1,+1]
            self.pos     = [ 0.5*(self.nodes[1].getPos() + self.nodes[0].getPos()) ]
            self.tangent = [ self.nodes[1].getPos() - self.nodes[0].getPos() ]
            self.area    = [ rot @ self.tangent ]
        elif len(self.nodes) == 3:
            # this includes the weight factor for two-point integration over [-1,+1]
            x0 = 0.4553418012614795 * self.nodes[0].getPos() \
                 + 0.6666666666666667 * self.nodes[1].getPos() \
                 - 0.1220084679281462 * self.nodes[2].getPos()
            x1 = -0.1220084679281462 * self.nodes[0].getPos() \
                 - 0.6666666666666667 * self.nodes[1].getPos() \
                 + 0.4553418012614795 * self.nodes[2].getPos()
            t0 = -1.077350269189626 * self.nodes[0].getPos() \
                 + 1.154700538379252 * self.nodes[1].getPos() \
                 - 0.07735026918962576 * self.nodes[2].getPos()
            t1 = 0.07735026918962576 * self.nodes[0].getPos() \
                 - 1.154700538379252 * self.nodes[1].getPos() \
                 + 1.077350269189626 * self.nodes[2].getPos()
            self.pos     = [ x0, x1 ]
            self.tangent = [ t0, t1 ]
            self.area    = [ rot @ t0, rot @ t1 ]
        else:
            msg = "{} requires 2 or 3 points - {} given".format(self.__class__.__name__, len(self.nodes))
            raise TypeError(msg)

    def computeNodalForces(self):
        """
        Implementation of nodal force calculation from surface loads.

        The surface is described  as

        .. math::

            {\\boldsymbol\\varphi}(s) = \\sum_k \\phi_k(s) {\\bf X}_k

        where :math:`\\phi_k(s)` is the k-th shape function on `[-1,1]`,
        and :math:`{\\bf X}_k` is the coordinate of the k-th node.

        The nodal forces are obtained as

        .. math::

            {\\bf P}_i = \\int_{-1}^{+1} \\phi_i(s)
            \\left(
                p_n \, {\\boldsymbol\\varphi}_{,s} \\times {\\bf k}
                + p_s \, {\\boldsymbol\\varphi}_{,s}
            \\right) ds

        """
        if len(self.nodes) == 2:
            pn = self.load[0]
            ps = self.load[1]
            forces = np.outer( (pn * self.area + ps * self.tangent), np.array([ 0.50000, 0.50000 ]))

        elif len(self.nodes) == 3:
            pn = self.load[0]
            ps = self.load[1]


            forces = np.array([[ 0.4553418012614795, 0.6666666666666667, -0.1220084679281462 ],
                               [ -0.1220084679281462, 0.6666666666666667, 0.4553418012614795 ]])

        else:
            msg = "{} requires 2 or 3 points - {} given".format(self.__class__.__name__, len(self.nodes))
            raise TypeError(msg)

        for i, node in enumerate(self.nodes):
            node.setLoad(forces[:,i])  # this may need a different approach to account for the global load factor


