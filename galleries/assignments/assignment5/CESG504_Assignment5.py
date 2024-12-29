"""
==========================================================
Bending a cantilever beam
==========================================================

* Using PatchMesher to model the beam
* representing a moment through nodal forces

"""
import math
import sys
import numpy as np
import matplotlib.pyplot as plt

from femedu.examples import Example

from femedu.domain import System, Node
from femedu.solver import NewtonRaphsonSolver
#from femedu.elements.linear import Quad
from femedu.elements.linear import ReducedIntegrationQuad as Quad
#from ReducedIntegrationQuad import ReducedIntegrationQuad as Quad
from femedu.materials import PlaneStress
from femedu.mesher import *

# %%
# Using a regular mesh
# ------------------------

class BeamModel(Example):

    # sphinx_gallery_start_ignore
    # sphinx_gallery_thumbnail_number = 2
    def docString(self):
        s = """
    ## Bending a cantilever beam

    Using PatchMesher to model the beam

    Author: Peter Mackenzie-Helnwein 
    """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):
        # ========== setting mesh parameters ==============

        Nx = 16         # number of elements in the mesh
        Ny =  8         # number of elements in the mesh
        Lx = 8000.0     # length of plate in the x-direction in cm
        Ly = 1000.0     # length of plate in the y-direction in cm
        Moment = -50000000. # moment in N.mm

        # ========== setting material parameters ==============

        params = dict(
            E  = 70000.,  # Young's modulus in MPa
            nu = 0.350,  # Poisson's ratio
            t  =   350.  # thickness of the plate
        )

        Iz = params['t']*Ly**3/12.

        # ========== setting load parameters ==============

        px  = 0.0  # uniform load normal to x=Lx
        py  = 0.0  # uniform load normal to y=Ly
        pxy = 0.0  # uniform shear load on x=L

        # ========== setting analysis parameters ==============

        target_load_level = 1.00  # reference load
        max_steps = 2  # number of load steps: 2 -> [0.0, 1.0]

        # define a list of target load levels
        load_levels = np.linspace(0, target_load_level, max_steps)

        #
        # ==== Build the system model ====
        #

        model = System()
        model.setSolver(NewtonRaphsonSolver())

        # create nodes

        mesher = PatchMesher(model, (0., -Ly/2), (Lx, -Ly/2), (Lx, Ly/2), (0., Ly/2), (0.3*Lx, -Ly/2), None, (0.4*Lx, Ly/2), None)
        nodes, elements = mesher.quadMesh(Nx, Ny, Quad, PlaneStress(params))

        # define support(s)

        ## find nodes at y==0 and x==0

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[0], 0.0):
                node.fixDOF('ux', 'uy')  # fix left side

        # ==== complete the reference load ====

        Xo = np.array([Lx, 0.0])
        No = np.array([1.0, 0.0])

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[0], Lx):
                # locate the node at the centerline
                if math.isclose(X[1], 0.0):
                    end_node = node
                # load the end faces
                for elem in node.elements:
                    for face in elem.faces:
                        for x, area in zip(face.pos, face.area):
                            if np.abs((x - Xo) @ No) < 1.0e-2 and No @ area / np.linalg.norm(area):
                                y = x[1]
                                px = -(y * Moment / Iz)*params['t']
                                face.setLoad(px, -pxy)

        Xo = np.array([0.0, Ly/2.])
        No = np.array([0.0, 1.0])

        for node in nodes:
            X = node.getPos()
            if math.isclose(X[1], Ly/2.):
                for elem in node.elements:
                    for face in elem.faces:
                        for x, area in zip(face.pos, face.area):
                            if np.abs((x - Xo) @ No) < 1.0e-2 and No @ area / np.linalg.norm(area):
                                face.setLoad(-py, 0.0)

        # model.report()

        # # set up a recorder
        # model.initRecorder(variables=['ux', 'uy'], nodes=[end_node])
        # model.startRecorder()

        model.plot(factor=0, title="undeformed system", show_bc=1, show_loads=1)

        for lf in np.linspace(0.0, target_load_level, max_steps):
            model.setLoadFactor(lf)
            model.solve(verbose=True, tolerance=1.e-3)

        # find elements that contain gauss points at x=256.4

        Xs   = []
        Ys   = []
        SigXXs = []
        SigYYs = []
        SigXYs = []

        for elem in elements:
            # find centroid of the element
            xo = np.array([0.0,0.0])
            for node in elem.nodes:
                X = node.getPos()
                xo += X
            xo /= len(elem.nodes)
            ##print(xo)

            # pick those around 250 < xo[0] < 280
            if xo[0] < Lx/2 and xo[0] > (Lx/2-Lx/Nx):
                #print(elem)
                # find the average stress
                sigXavg = np.zeros(3)
                for stress in elem.stress:
                    ##print(stress['xx'])
                    sigXavg[0] += stress['xx']
                    sigXavg[1] += stress['yy']
                    sigXavg[2] += stress['xy']
                sigXavg /= len(elem.stress)
                print(f'average stress at x = {xo[0]}, y = {xo[1]} is {sigXavg}')

                Xs.append(xo[0])
                Ys.append(xo[1])
                SigXXs.append(sigXavg[0])
                SigYYs.append(sigXavg[1])
                SigXYs.append(sigXavg[2])

        # adding the exact solution
        Area = Ly
        Iz = Ly**3./12.
        # Mz = -py*(Lx - Xs[0])**2/2.
        # Vy = -py*(Lx - Xs[0])
        Mz = Moment
        Vy = 0.0
        y = np.linspace(-Ly/2, Ly/2, 10)
        tauBeam = Vy/Area * 1.5 * (1. - (2.*y/Ly)**2)
        sigBeam = -y * Mz/Iz

        # plot section stress
        plt.plot(np.zeros_like(y),y,'-k',lw=2,label=None)
        plt.fill_betweenx(y,sigBeam,ls='-',color='#ffaaaa88',label='$\\sigma_{xx}$ beam solution')
        plt.fill_betweenx(y,tauBeam,ls='-',color='#aaffaa88',label='$\\sigma_{xy}$ beam solution')
        plt.plot(SigXXs,Ys,'or',label='$\\sigma_{xx}$ element average')
        plt.plot(SigYYs,Ys,'ob',label='$\\sigma_{yy}$ element average')
        plt.plot(SigXYs,Ys,'og',label='$\\sigma_{xy}$ element average')
        plt.grid(True)
        plt.legend()
        plt.title(f'Stress at x={Xs[0]} for {Nx}/{Ny} elements per length/height')
        plt.xlabel('stress $\\sigma_{ii}$')
        plt.ylabel('y')
        plt.savefig('comparison.png')
        plt.show()

        #model.report()

        model.plot(factor=1000., show_bc=1, show_loads=1, show_reactions=1)
        #
        # model.valuePlot('ux')
        # model.valuePlot('uy', show_mesh=True)

        # create a history plot for the end node

        # model.historyPlot('lam', ['ux', 'uy'], nodes=[end_node, end_node])
        # model.historyPlot(('ux', end_node), 'uy', node=end_node)

        print(end_node)



# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":

    ex = BeamModel()
    ex.run()


