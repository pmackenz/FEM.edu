from .Constraint import Constraint

class TieNodes(Constraint):
    """
    Use this to tie some or all dofs of two or more nodes.

    The program will select a **lead** node and assign all remaining nodes
    as **followers**.

    .. note::

        Most textbooks refer to such nodes as "master" and "slave" nodes,
        which we avoid on purpose due to their historic use in slavery
        and reappearance in modern day racism.

        Please view this is an apology for mistakes of the past
        and help us eradicate racism along the way.


    """

    def __init__(self, *nodes, dofs=('all',), **kwargs):
        """
        :param nodes: two or mode :py:class:`Node` objects
        :param dofs:  (optional) a list of dof-codes to be tied
                The default code :py:code:`"all"` will tie all dofs shared by the nodes.

        """
        params = {}
        params['dofs'] = dofs

        if 'all' in dofs:
            self.neqns = len(nodes[0].ndofs)  # needs fancier logic for multi-node ties
        else:
            self.neqns  = len(dofs)

        super(TieNodes, self).__init__(nodes, params, **kwargs)
