from .Constraint import Constraint

class SPconstraint(Constraint):
    """
    wrapper for a single point constraint
    """

    def __init__(self, node, dofs, a, u0, u1=0, lam0=0, **kwargs):
        """
        :param node: one :py:class:`Node` object
        :param dofs: list of dof codes for which to apply those constraints
        :param a: list or array of factors associated with each listed **dof**
        """
        params = {}
        params['u0']   = u0
        params['u1']   = u1
        params['lam0'] = lam0
        params['dofs'] = dofs
        params['a']    = a

        self.neqns  = 0

        super(SPconstraint, self).__init__([node], params, **kwargs)
