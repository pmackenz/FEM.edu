import numpy as np
import scipy as sc
from copy import deepcopy

from .Node import Node
from ..elements import Element


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
        """

        :param dir1: raw first direction
        :param dir2: raw second direction
        :param axis: axial rotation vector
        """

        self.dir1  = dir1          # user input raw first direction
        self.dir2  = dir2          # user input raw second direction
        self.axis  = axis          # user input axial rotation vector

        self.known_vectors = []    # create a list of vectors that should be transformed
        self.clients = []          # holds pointers to client nodes/elements

        dir1_is_vector = (isinstance(dir1, np.ndarray) or isinstance(dir1, list))
        dir2_is_vector = (isinstance(dir2, np.ndarray) or isinstance(dir2, list))
        axis_is_vector = (isinstance(self.axis, np.ndarray) or isinstance(self.axis, list))

        have_dirs = dir1_is_vector and dir2_is_vector
        msg = ""
        model_dim = 0

        if have_dirs:
            if len(dir1) == len(dir2):
                model_dim = len(dir1)
            else:
                msg = "Cannot define transformation from vectors of mismatched length: len(dir1)={}  len(dir2)={}".format(len(dir1),len(dir2))
        elif axis_is_vector:
            model_dim = len(axis)
        else:
            msg = f"Cannot define transformation from axis={axis}. Need a vector but got {type(axis)}"

        if model_dim not in (2,3):
            msg = f"Transformation only defined for 2d and 3d but got model dimension {model_dim}"

        if msg:
            raise TypeError(msg)

        if have_dirs and model_dim == 2:
            self.init_from_dirs_2d(have_dirs)
        elif have_dirs and model_dim == 3:
            self.init_from_dirs_3d(have_dirs)

    def init_from_dirs_2d(self, have_dirs):

        if have_dirs:
            g1 = np.array(self.dir1)
            g1 /= np.linalg.norm(g1)
            g2 = np.array(self.dir2)
            g2 -= (g1@g2) * g1
            g2 /= np.linalg.norm(g2)
            self.T = np.vstack((g1,g2)).T
        elif np.isreal(self.axis):
            w3 = self.axis
            omega = np.array([[0., -w3],[w3, 0.]])
            self.T = sc.linalg.expm(omega)
        else:
            self.T = np.identity(2)

    def init_from_dirs_3d(self, have_dirs):

        axis_is_vector = (isinstance(self.axis, np.ndarray) or isinstance(self.axis, list))

        if have_dirs:
            g1 = np.array(self.dir1)
            g1 /= np.linalg.norm(g1)
            g2 = np.array(self.dir2)
            g2 -= (g1@g2) * g1
            g2 /= np.linalg.norm(g2)
            g3 = np.linalg.cross(g1,g2)
            self.T = np.vstack((g1,g2,g3)).T
        elif axis_is_vector:
            w1,w2,w3 = self.axis
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

    def registerClient(self, client):
        """
        Register a Node or Element object with this transformation.

        :param client: pointer to the node/element to be registered
        """
        if isinstance(client, (Node, Element)):
            if client not in self.clients:
                self.clients.append(client)
                self.refreshMaps()
        else:
            msg = "Only Nodes and Elements can interact with a Transformation: {} given".format(type(client).__name__)
            raise TypeError(msg)

    def addVector(self, vecdef):
        """
        Register a vector to be transformed by this Transformation object.

        This method is most commonly used internally, though may be used to add user-defined
        vector dofs to any (nodal) transformation.

        :param vecdef: an array-like definitions of a vector.

        **Examples for vecdef**

        * 2d solids use :code:`['ux,'uy']`
        * 3d solids use :code:`['ux','uy','uz']`
        * 3d Frames use :code:`['ux','uy','uz']` and :code:`['rx','ry','rz']`
        * your own vector dofs may use :code:`['myx','myy','myz']` if 'myx' represents the x-component of a :code:`my`-vector. These dofs must be requested by the constructor of your element.
        """
        newvec = list(vecdef)  # ensure object is a list
        if newvec not in self.known_vectors:
            self.known_vectors.append({'target':newvec, 'map':None})
            self.refreshMaps()

    def refreshMaps(self):
        r"""
        Reinitialize dof-maps for vector transformations.

        .. note::

            This method shall be called every time a vector is registered for transformation **and**
            if an element is added to a node since a new element addition may change the list of
            locally available dofs!

        The map structure is defined as

        .. code::

            self.known_vectors = [
                # one entry per vector that needs transformation
                {
                    'dofs': [ dofs_defining_this_vector ],
                    'map': {
                                # one entry per client
                                Node|Element: [ idx_for_dof[0], idx_for_dof[1], ... ],
                                ...
                           }
                },
                ...
            ]

        The map for any requested vector can be obtained using

        .. code::

            node = ...     # pointer to node of interest
            elem = ...     # pointer to elem of interest

            for vec in self.known_vectors:

                vector_dofs = vec['dofs']  # list of named dofs forming this vector

                if node in vec['map']:
                    node_map = vec['map'][node]

                if elem in vec['map']:
                    elem_map = vec['map'][elem]

        """
        for vec in self.known_vectors:
            map = {}
            for client in self.clients:
                if isinstance(client, Node):
                    map[client] = client.getIdx4DOFs(dofs=vec['target'], local=True)
                elif isinstance(client, Element):
                    map[client] = client.getDofs()
            vec['map'] = map

    def getT(self):
        r"""
        returns the transformation matrix that maps components from a local 3d vector to the global system

        **Usage**

        .. code::

            T = transform.getT()

            # vectors
            Uglobal = T @ Ulocal          # local to global for a vector U
            Vlocal = T.T @ Vglobal        # global to local for a vector V
            Vlocal = Vglobal @ T          # alternative global to local for a vector V

            # matrices
            Mglobal = T @ Mlocal @ T.T    # local to global
            Mlocal  = T.T @ Mglobal @ T   # global to local

        """
        return self.T

    #
    # acting on 2d|3d vectors|matrices
    #

    def vectorToGlobal(self, U, dof_list=[]):
        r"""
        :param U: 2d or 3d vector
        """
        return U @ self.T

    def vectorToLocal(self, U, dof_list=[]):
        r"""
        :param U: 2d or 3d vector
        """
        return self.T @ U

    def matrixToGlobal(self, M, dof_list=[]):
        r"""
        :param M: [2x2] or [3x3] matrix
        """
        return self.T @ M @ self.T.T

    def matrixToLocal(self, M, dof_list=[]):
        r"""
        :param M: [2x2] or [3x3] matrix
        """
        return self.T.T @ M @ self.T

    #
    # acting on the entire Node vector
    #

    def v2g(self, U, caller=None):
        r"""
        :param U: [ndof] vector w/ ndof = number of dofs at Node||Element
        """
        Utransformed = deepcopy(U)

        maps = self._caller_check(caller)

        for map in maps:
            local_vec  = Utransformed[map]
            global_vec = self.T @ local_vec
            Utransformed[map] = global_vec

        return Utransformed

    def v2l(self, U, caller=None):
        r"""
        :param U: [ndof] vector w/ ndof = number of dofs at Node||Element
        """
        Utransformed = deepcopy(U)

        maps = self._caller_check(caller)

        for map in maps:
            local_vec  = Utransformed[map]
            global_vec = self.T.T @ local_vec
            Utransformed[map] = global_vec

        return Utransformed

    def m2g(self, M, caller=None):
        r"""
        :param M: [ndof x ndof] matrix w/ ndof = number of dofs at Node||Element
        """
        maps = self._caller_check(caller)

        return self.T @ M @ self.T.T

    def m2l(self, M, caller=None):
        r"""
        :param M: [ndof x ndof] matrix w/ ndof = number of dofs at Node||Element
        """
        Mtransformed = deepcopy(M)

        maps = self._caller_check(caller)

        for map in maps:
            map = np.array(map, dtype=int)
            local_vec  = Mtransformed[map,:]
            global_vec = self.T.T @ local_vec
            Mtransformed[map,:] = global_vec

        return Mtransformed

        # return self.T.T @ M @ self.T

    def _caller_check(self, caller):
        if not caller:
            msg = "caller must be provided when requesting a vector|matrix transformation from {}".format(self.__class__.__name__)
            raise TypeError(msg)

        maps = []

        for vec in self.known_vectors:
            if caller in vec['map']:
                maps.append(vec['map'][caller])

        return maps


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