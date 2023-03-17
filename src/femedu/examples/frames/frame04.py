"""
Example: frame04

Buckling of a building frame

modeled using a 2D frame element

.. list-table:: setting given parameters

    * - N  = 2
      - number of elements
    * - L  = 100.0
      - column length
    * - EA = 2000000.0
      - axial stiffness
    * - EI = 21000.0
      - flexural stiffness
    * - w  = 0.1
      - applied lateral load


"""
from ...examples.Example import *

from ...domain import *
from ...solver.NewtonRaphsonSolver import *
from ...elements.Frame2D import *
from ...materials.ElasticSection import *


class ExampleFrame04(Example):

    def docString(self):
        s = """
    Buckling of a building frame

    degrees of freedom:
        ux ... horizontal displacement
        uy ... vertical displacement
        rz ... rotation, theta
        
    Author: Peter Mackenzie-Helnwein 
    """
        return s

    def problem(self):
        # initialize a system model

        N  = 8     # number of elements

        B = 720.
        H = 720.

        E  = 29000.0
        A = 150.0
        I = 250.0

        w = 0.10
        load_at_nodes_only = False # set to True to apply equivalent nodal forces and moments

        Ph = 0.01      # additional horizontal load per floor
        Ph = 0.10      # additional horizontal load per floor
        Ph = 1.00      # additional horizontal load per floor
        Ph = 0.00      # additional horizontal load per floor

        # ========== setting global parameters ==============

        target_load_level = 33.
        max_steps = 10
        load_levels = np.linspace(0, target_load_level, max_steps)

        # ========= build your structural model =============

        model = System()
        model.setSolver(NewtonRaphsonSolver())

        x0 = 0.0
        x1 = B / 3
        x2 = 2 * B / 3
        x3 = B

        y0 = 0.0
        y1 = H / 4
        y2 = 2 * H / 4
        y3 = 3 * H / 4
        y4 = H

        X10 = Node(x0, y0)
        X11 = Node(x0, y1)
        X12 = Node(x0, y2)
        X13 = Node(x0, y3)
        X14 = Node(x0, y4)

        X20 = Node(x1, y0)
        X21 = Node(x1, y1)
        X22 = Node(x1, y2)
        X23 = Node(x1, y3)
        X24 = Node(x1, y4)

        X30 = Node(x2, y0)
        X31 = Node(x2, y1)
        X32 = Node(x2, y2)
        X33 = Node(x2, y3)
        X34 = Node(x2, y4)

        X40 = Node(x3, y0)
        X41 = Node(x3, y1)
        X42 = Node(x3, y2)
        X43 = Node(x3, y3)
        X44 = Node(x3, y4)

        model.addNode(X10,X11,X12,X13,X14)
        model.addNode(X20,X21,X22,X23,X24)
        model.addNode(X30,X31,X32,X33,X34)
        model.addNode(X40,X41,X42,X43,X44)

        # columns

        params = {'E': E, 'A': A, 'I': I}

        C11 = Frame2D(X10, X11, ElasticSection(params))
        C12 = Frame2D(X11, X12, ElasticSection(params))
        C13 = Frame2D(X12, X13, ElasticSection(params))
        C14 = Frame2D(X13, X14, ElasticSection(params))

        model.addElement(C11,C12,C13,C14)

        params = {'E': E, 'A': 2*A, 'I': 1.5*I}

        C21 = Frame2D(X20, X21, ElasticSection(params))
        C22 = Frame2D(X21, X22, ElasticSection(params))
        C23 = Frame2D(X22, X23, ElasticSection(params))
        C24 = Frame2D(X23, X24, ElasticSection(params))

        model.addElement(C21,C22,C23,C24)

        C31 = Frame2D(X30, X31, ElasticSection(params))
        C32 = Frame2D(X31, X32, ElasticSection(params))
        C33 = Frame2D(X32, X33, ElasticSection(params))
        C34 = Frame2D(X33, X34, ElasticSection(params))

        model.addElement(C31,C32,C33,C34)

        params = {'E': E, 'A': A, 'I': I}

        C41 = Frame2D(X40, X41, ElasticSection(params))
        C42 = Frame2D(X41, X42, ElasticSection(params))
        C43 = Frame2D(X42, X43, ElasticSection(params))
        C44 = Frame2D(X43, X44, ElasticSection(params))

        model.addElement(C41,C42,C43,C44)

        # floors

        params = {'E': E, 'A': A, 'I': 3*I}

        F11 = Frame2D(X11, X21, ElasticSection(params))
        F12 = Frame2D(X21, X31, ElasticSection(params))
        F13 = Frame2D(X31, X41, ElasticSection(params))

        model.addElement(F11,F12,F13)

        F21 = Frame2D(X12, X22, ElasticSection(params))
        F22 = Frame2D(X22, X32, ElasticSection(params))
        F23 = Frame2D(X32, X42, ElasticSection(params))

        model.addElement(F21,F22,F23)

        F31 = Frame2D(X13, X23, ElasticSection(params))
        F32 = Frame2D(X23, X33, ElasticSection(params))
        F33 = Frame2D(X33, X43, ElasticSection(params))

        model.addElement(F31,F32,F33)

        F41 = Frame2D(X14, X24, ElasticSection(params))
        F42 = Frame2D(X24, X34, ElasticSection(params))
        F43 = Frame2D(X34, X44, ElasticSection(params))

        model.addElement(F41,F42,F43)

        # fixities
        X10.fixDOF('ux','uy','rz')   # fixed
        X20.fixDOF('ux','uy','rz')   # fixed
        X30.fixDOF('ux','uy','rz')   # fixed
        X40.fixDOF('ux','uy','rz')   # fixed

        # reference load
        #Pcr = np.pi**2 * EI / L**2
        model.resetLoad()            # size load vector and initialize
        #model.addLoad(Xn, -Pcr, dof=0) # add a horizontal force (first dof only) ; remember C-style indexing: 0,1,...,(n-1)

        if load_at_nodes_only:

            # floor loading as nodal loads ...

            Pe = w * B/3
            Mi = w * (B/3)**2 /12

            X11.addLoad([-Pe/2., -Mi],['uy','rz'])
            X21.addLoad([-Pe/2.,  0.],['uy','rz'])
            X31.addLoad([-Pe/2.,  0.],['uy','rz'])
            X41.addLoad([-Pe/2.,  Mi],['uy','rz'])

            X12.addLoad([-Pe/2., -Mi],['uy','rz'])
            X22.addLoad([-Pe/2.,  0.],['uy','rz'])
            X32.addLoad([-Pe/2.,  0.],['uy','rz'])
            X42.addLoad([-Pe/2.,  Mi],['uy','rz'])

            X13.addLoad([-Pe/2., -Mi],['uy','rz'])
            X23.addLoad([-Pe/2.,  0.],['uy','rz'])
            X33.addLoad([-Pe/2.,  0.],['uy','rz'])
            X43.addLoad([-Pe/2.,  Mi],['uy','rz'])

            X14.addLoad([-Pe/2., -Mi],['uy','rz'])
            X24.addLoad([-Pe/2.,  0.],['uy','rz'])
            X34.addLoad([-Pe/2.,  0.],['uy','rz'])
            X44.addLoad([-Pe/2.,  Mi],['uy','rz'])

        else:

            # floor loading as distributed loads ...

            F11.setDistLoad(-w)
            F12.setDistLoad(-w)
            F13.setDistLoad(-w)

            F21.setDistLoad(-w)
            F22.setDistLoad(-w)
            F23.setDistLoad(-w)

            F31.setDistLoad(-w)
            F32.setDistLoad(-w)
            F33.setDistLoad(-w)

            F41.setDistLoad(-w)
            F42.setDistLoad(-w)
            F43.setDistLoad(-w)


        # wind load ...

        X11.addLoad([Ph],['ux'])   # horizontal load
        X12.addLoad([Ph],['ux'])   # horizontal load
        X13.addLoad([Ph],['ux'])   # horizontal load
        X14.addLoad([Ph/2],['ux']) # horizontal load


        # show model information
        print(model)

        print("\n==== perform the analysis ===\n")

        # * apply the load in multiple smaller load steps

        # set up data recorder
        model.initRecorder()
        model.trackStability(True)

        # initialize the analysis:
        model.resetDisplacements()   # set U to all zeros
        model.setLoadFactor(0.0)     # define a known equilibrium solution

        model.startRecorder()

        detKt   = []
        lambdas = []

        # solve for all load_levels
        for loadfactor in load_levels:

            # define node X2 as the controled node; downward direction is prescribed:
            model.setLoadFactor(loadfactor)
            model.solve(verbose=True)

            # stability check
            lambdas.append(model.loadfactor)
            detKt.append(model.solver.checkStability())

            # report results
            print('+')
            #model.report()

            print("\n=== next load level ===\n")


        #
        # ==== create some nice plots ===
        #

        model.report()

        print(lambdas)
        print(detKt)

        plt.plot(lambdas,detKt,'--*r')
        plt.grid(True)
        plt.xlabel('Load factor, $ \lambda $')
        plt.ylabel("Stability index, $ {det}\: {\\bf K}_t $")
        plt.show()

        model.plot(factor=10.0, filename="frame4_deformed.png")

        model.beamValuePlot("F", filename="frame4_force.png")
        model.beamValuePlot("V", filename="frame4_shear.png")
        model.beamValuePlot("M", filename="frame4_moment.png")

        model.plotBucklingMode(factor=100., mode=0, filename="frame4_buckling_mode0.png")
        model.plotBucklingMode(factor=100., mode=1, filename="frame4_buckling_mode1.png")
        model.plotBucklingMode(factor=100., mode=2, filename="frame4_buckling_mode2.png")
        model.plotBucklingMode(factor=100., mode=3, filename="frame4_buckling_mode3.png")
