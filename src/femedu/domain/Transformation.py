

class Transformation():
    """
    Nodal transformation class

    A transformation may be attached to Nodes or to Elements.

    - When attached to a Node, the transformation applies to
      all fixities and all loads defined at that node.
    - When attached to an element, the user may specify a node for which that transformation
      shall be used.  If no dode is specified, the transformation will be applied at all
      attached nodes.

    """

    def __init__(self):
        pass

    def vectorToGlobal(self, U, dof_list=[]):
        return U

    def vectorToLocal(self, U, dof_list=[]):
        return U

    def matrixToGlobal(self, U, dof_list=[]):
        return U

    def matrixToLocal(self, U, dof_list=[]):
        return U
