"""
==========================================================
Pulling a notched plate : using quads
==========================================================

Using PatchMesher to model a quarter of the plate

"""
import math
import sys

import matplotlib.pyplot as plt
import numpy as np

from femedu.examples import Example

from femedu.domain import System, Node
# from femedu.solver import NewtonRaphsonSolverSparse as NewtonRaphsonSolver
from femedu.solver import NewtonRaphsonSolver

from femedu.elements.linear import Quad
# from femedu.elements.linear import Quad9 as Quad
from femedu.materials import PlaneStress
from femedu.mesher import *


class ExampleFinal(Example):

    # sphinx_gallery_start_ignore
    # sphinx_gallery_thumbnail_number = 3
    def docString(self):
        s = """
    ## Pulling a notched plate : using quads

    The full plate has a dimension of (2*Lx)*(2*Ly) , Ly < Lx,
    and a circular hole with radius R < Ly.

    Using PatchMesher to model a quarter of the plate,
    applying symmetry conditions, and face loading on the
    free surface x=a.

    Author: Peter Mackenzie-Helnwein
    """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):
        # self.modelD(4., 1.00)
        # self.modelD(8., 2.00)
        # self.modelD(13., 1.00)
        # self.modelD(17., 3.00)
        # self.modelD(21., 2.00)
        # self.modelD(25., 1.00)
        self.modelD(25., 1.00)
        # self.modelD(25., 30.00)
        # self.modelD(28., 3.00)
        # self.modelD(35., 2.00)
        # self.modelD(56., 2.00)
        # self.modelD(120., 2.00)
        # self.modelD(120., 2.00)

    def modelD(self, D, h):
        """
        :params D: diameter of drill in mm
        :params h: estimated element size near the notch in mm
        """

        # ====== units ======
        mm  = 1.
        N   = 1.
        m   = 1000 * mm
        kN  = 1000 * N
        Pa  = N/m**2
        kPa = 1000 * Pa
        MPa = 1000 * kPa
        GPa = 1000 * MPa

        # ========== setting model dimensions ==============

        Lx = 120.0 * mm   # length of plate in the x-direction
        Ly =  60.0 * mm   # length of plate in the y-direction
        R  = D/2

        # ========== setting material parameters ==============

        params = dict(
            E  = 200. * GPa,    # Young's modulus
            nu = 0.30,          # Poisson's ratio
            t  = 1.00 * mm      # thickness of the plate
        )

        # ========== setting load parameters ==============

        px  = 50.0 * MPa        # uniform load normal to x=const

        # ========== setting mesh parameters ==============

        Nx = int(Lx/h/3) + 1     # number of elements in the mesh
        Ny = int(Ly/h) + 1        # number of elements in the mesh

        # ========== setting analysis parameters ==============

        target_load_level = 1.00     # reference load
        max_steps = 2                # number of load steps: 2 -> [0.0, 1.0]

        # define a list of target load levels
        load_levels = np.linspace(0, target_load_level, max_steps)

        #
        # ==== Build the system model ====
        #

        model = System()
        model.setSolver(NewtonRaphsonSolver())

        # create nodes

        #  5--------------6--------------7
        #  |             / \             |
        #  9     10     11  \            |
        #  |           /     \           |
        #  3----8----4        \          |
        #             \   12   \         |
        #              \        \        |
        #               \        \       |
        #                0---13---1------2

        xm = 0.0
        ym = -Ly/2 - R

        y0 = -Ly
        y1 = -Ly
        y2 = -Ly

        y3 = ym + R*np.cos(np.radians( 0.0))
        y8 = ym + R*np.cos(np.radians(22.5))
        y4 = ym + R*np.cos(np.radians(45.0))

        y5 = 0.0
        y6 = 0.0
        y7 = 0.0

        x2 = Lx

        x3 = xm + R*np.sin(np.radians( 0.0))
        x8 = xm + R*np.sin(np.radians(22.5))
        x4 = xm + R*np.sin(np.radians(45.0))

        x5 = 0.0
        x7 = Lx
        x6 = 0.25 * (x5 + x7)

        x0 = x4 + (y4 - y0)
        x1 = 0.45 * (x0 + x2)

        x9  = 0.7*x3 + 0.3*x5
        x10 = 0.7*x8 + 0.3*0.5*(x5+x6)
        x11 = 0.7*x4 + 0.3*x6
        x12 = 0.7*0.5*(x0+x4) + 0.3*0.5*(x1+x6)
        x13 = 0.7*x0 + 0.3*x1

        y9  = 0.7*y3 + 0.3*y5
        y10 = 0.7*y8 + 0.3*0.5*(y5+y6)
        y11 = 0.7*y4 + 0.3*y6
        y12 = 0.7*0.5*(y0+y4) + 0.3*0.5*(y1+y6)
        y13 = 0.7*y0 + 0.3*y1

        pts = (
            ( x0, y0),
            ( x1, y1),
            ( x2, y2),
            ( x3, y3),
            ( x4, y4),
            ( x5, y5),
            ( x6, y6),
            ( x7, y7),
            ( x8, y8),
            ( x9, y9),
            ( x10, y10),
            ( x11, y11),
            ( x12, y12),
            ( x13, y13)
        )

        mesher1 = PatchMesher(model, pts[3], pts[4], pts[6], pts[5], pts[8], pts[11], None, pts[9], pts[10])
        nodes1, elements1 = mesher1.quadMesh(Nx, Ny, Quad, PlaneStress(params))

        mesher2 = PatchMesher(model, pts[0], pts[1], pts[6], pts[4], pts[13], None, pts[11], None, pts[12])                                               # center node
        nodes2, elements2 = mesher2.quadMesh(Ny, 2*Nx, Quad, PlaneStress(params))

        mesher3 = PatchMesher(model, pts[1], pts[2], pts[7], pts[6])
        nodes3, elements3 = mesher3.quadMesh(Nx, 2*Nx, Quad, PlaneStress(params))

        mesher1.tie(mesher2)
        mesher2.tie(mesher3)

        nodes    = nodes1    + nodes2    + nodes3
        elements = elements1 + elements2 + elements3

        # define support(s)

        ## find nodes at y==0 and x==0

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[0], 0.0):
                node.fixDOF('ux')    # horizontal support left side
            if math.isclose(X[1], 0.0):
                node.fixDOF('uy')    # vertical support at y==0

        # ==== complete the reference load ====

        for elem, face in model.findFacesAlongLine((Lx, 0.0), (0.0, 1.0)):
            face.setLoad(px, 0.0)

        # model.plot(factor=0, title="undeformed system", show_bc=1, show_loads=1)

        model.solve(verbose=True)

        # model.plot(factor=100., show_bc=1, show_loads=1, show_reactions=1)
        #
        # model.valuePlot('sxx', show_mesh=0)
        # model.valuePlot('syy', show_mesh=0)
        # model.valuePlot('sxy', show_mesh=0)

        # create a stress section

        data = []
        for node, _ in model.findNodesAlongLine((0.,0.),(0.,-1.)):
            y = node.getPos()[1]
            sxx = node.getMappedValue('sxx')
            syy = node.getMappedValue('syy')
            sxy = node.getMappedValue('sxy')
            data.append([y,sxx,syy,sxy])
        data = np.array(data)

        maxStress = data[:,1].max()
        SCF = maxStress / (2.*px)

        print(f"Number of elements:   {len(model.elements):12d}\n",
              f"Number of nodes:      {len(model.nodes):12d}\n",
              f"Degrees of freedom:   {model.solver.sdof:12d}\n",
              f"Max normal stress:           {maxStress:12.2f} MPa\n",
              f"Stress Concentration Factor: {SCF:12.3f}")

        # fig, ax = plt.subplots()
        # ax.plot(data[:,0],data[:,1],'r-',label=r'$\sigma_{xx}$')
        # ax.plot(-data[:,0],data[:,1],'r-',label=None)
        # ax.plot(data[:,0],data[:,2],'g-',label=r'$\sigma_{yy}$')
        # ax.plot(-data[:,0],data[:,2],'g-',label=None)
        # ax.plot(data[:,0],data[:,3],'b-',label=r'$\sigma_{xy}$')
        # ax.plot(-data[:,0],-data[:,3],'b-',label=None)
        # ax.legend()
        # ax.grid(True)
        # ax.set_title("Stress distribution at the narrow section")
        # ax.set_xlabel(r'$y$ (mm)')
        # ax.set_ylabel(r'stress $\sigma_{ij}$ (MPa)')
        # plt.show()


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":

    from timeit import default_timer as timer

    start = timer()
    ex = ExampleFinal()
    ex.run()
    end = timer()
    print(end - start) # time in seconds


