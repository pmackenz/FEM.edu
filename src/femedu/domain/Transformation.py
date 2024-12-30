import numpy as np
import scipy as sc

class Transformation():
    r"""
    Nodal transformation class

    A transformation may be attached to Nodes or to Elements.

    - When attached to a Node, the transformation applies to
      all fixities and all loads defined at that node.
    - When attached to an element, the user may specify a node for which that transformation
      shall be used.  If no node is specified, the transformation will be applied at all
      attached nodes.

    """

    def __init__(self, dir1=None, dir2=None, axis=None):

        self.dir1  = dir1
        self.dir2  = dir2
        self.axis  = axis

        dir1_is_vector = dir1 and (isinstance(dir1, np.ndarray) or isinstance(dir1, list))
        dir2_is_vector = dir2 and (isinstance(dir2, np.ndarray) or isinstance(dir2, list))
        axis_is_vector = axis and (isinstance(axis, np.ndarray) or isinstance(axis, list))

        if dir1_is_vector and dir2_is_vector:
            g1 = np.array(dir1)
            g1 /= np.linalg.norm(g1)
            g2 = np.array(dir2)
            g2 -= (g1@g2) * g1
            g2 /= np.linalg.norm(g2)
            g3 = np.linalg.cross(g1,g2)
            self.T = np.vstack((g1,g2,g3)).T
        elif axis_is_vector:
            w1,w2,w3 = axis
            omega = np.array([[0., -w3, w2],[w3, 0., -w1],[ -w2, w1, 0.]])
            self.T = sc.linalg.expm(omega)
        else:
            self.T = np.identity(3)

    def __str__(self):
        s = "{}(".format(self.__class__.__name__)
        if np.any(self.dir1):
            s += f"dir1={self.dir1}, "
        if np.any(self.dir2):
            s += f"dir2={self.dir2}, "
        if np.any(self.axis):
            s += f"axis={self.axis}, "
        s += ")"
        return s

    def getT(self):
        return self.T

    def vectorToGlobal(self, U, dof_list=[]):
        r"""

        """
        return U @ self.T

    def vectorToLocal(self, U, dof_list=[]):
        r"""

        """
        return self.T @ U

    def matrixToGlobal(self, U, dof_list=[]):
        r"""

        """
        return self.T @ U @ self.T.T

    def matrixToLocal(self, U, dof_list=[]):
        r"""

        """
        return self.T.T @ U @ self.T


if __name__ == "__main__":

    vecs = []
    vecs.append(np.array([1.,0.,0.]))
    vecs.append(np.array([0.,1.,0.]))
    vecs.append(np.array([0.,0.,1.]))

    matrices = []
    matrices.append(np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]))
    matrices.append(np.array([[10.,0.,0.],[0.,20.,0.],[0.,0.,30.]]))
    matrices.append(np.array([[3.,2.,1.],[2.,3.,2.],[1.,2.,3.]]))

    transforms = []
    transforms.append(Transformation(dir1=[1.,0.,0.], dir2=[0.,1.,0.]))   # identity transformation
    transforms.append(Transformation(axis=[0.,0.,np.pi/360]))             # rotate CCW by 1.0 degrees
    transforms.append(Transformation(axis=[0.,0.,np.pi/2]))             # rotate CCW by 90.0 degrees
    transforms.append(Transformation(dir1=[1.,1.,0.], dir2=[0.,1.,0.]))   # rotate CCW by 45 degrees
    transforms.append(Transformation(axis=[0.,0.,np.pi/4]))               # rotate CCW by 45 degrees (alt input)

    for k,T in enumerate(transforms):
        print(T)
        print("transformation map:\n", T.getT())
        for i,vec in enumerate(vecs):
            print(f"vec{i} : {vec}")
            print(f"vec{i} g->l {T.vectorToLocal(vec)}")
            print(f"vec{i} l->g {T.vectorToGlobal(vec)}")
        for i,mat in enumerate(matrices):
            print(f"mat{i} : \n{mat}")
            print(f"mat{i} g->l \n{T.matrixToLocal(mat)}")
            print(f"mat{i} l->g \n{T.matrixToGlobal(mat)}")
        print("---")