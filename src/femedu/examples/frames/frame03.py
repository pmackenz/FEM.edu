"""
Example: frame03

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


class ExampleFrame03(Example):

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

        E  = 20000
        EA = 2000000.0
        EI = 21000.0

        w = 0.0075

        Ph = 0.01      # additional horizontal load per floor
        Ph = 0.00      # additional horizontal load per floor

        # ========== setting global parameters ==============

        target_load_level = 1.00
        max_steps = 30
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

        params = {'E': E, 'A': EA/E, 'I': EI/E}

        C11 = Frame2D(X10, X11, ElasticSection(params))
        C12 = Frame2D(X11, X12, ElasticSection(params))
        C13 = Frame2D(X12, X13, ElasticSection(params))
        C14 = Frame2D(X13, X14, ElasticSection(params))

        model.addElement(C11,C12,C13,C14)

        params = {'E': E, 'A': 2*EA/E, 'I': 1.5*EI/E}

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

        params = {'E': E, 'A': EA/E, 'I': EI/E}

        C41 = Frame2D(X40, X41, ElasticSection(params))
        C42 = Frame2D(X41, X42, ElasticSection(params))
        C43 = Frame2D(X42, X43, ElasticSection(params))
        C44 = Frame2D(X43, X44, ElasticSection(params))

        model.addElement(C41,C42,C43,C44)

        # floors

        params = {'E': E, 'A': EA/E, 'I': 3*EI/E}

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

        X11.addLoad([Ph],['ux'])   # horizontal load
        X12.addLoad([Ph],['ux'])   # horizontal load
        X13.addLoad([Ph],['ux'])   # horizontal load
        X14.addLoad([Ph/2],['ux']) # horizontal load


        # show model information
        print(model)

        model.solve(verbose=True)

        #model.report()

        model.plot(factor=25.0)

        model.beamValuePlot("F")
        model.beamValuePlot("M")
        model.beamValuePlot("V")

        return
