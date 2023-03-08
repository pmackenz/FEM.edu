from ...examples.Example import *

from ...domain import *
from ...solver.NewtonRaphsonSolver import *
from ...elements.Frame2D import *
from ...materials.ElasticSection import *


class ExampleFinal01(Example):

    def docString(self):
        s = """
    CESG 507 - 2023 - Final Exam Problem
    """
        return s

    def problem(self):

        # ========== setting mesh parameters ==============
        N = 5                  # number of elements in the mesh
        alpha = 0.5            #alpha value
        gamma = 2              # stiffening factor
        L = 100.0              # column free length

        # ========== setting material parameters ==============
        params = dict(
            E = 20000.,    # Young's modulus
            A = 100.0,     # cross section area
            I = 10.0       # cross section moment of inertia
        )

        stiffer_params = dict(
            E = 20000.,    # Young's modulus
            A = 100.0,     # cross section area
            I = gamma * 10.0       # cross section moment of inertia
        )

        # ========== setting load parameters ==============
        w   = 0.01        # uniform lateral load on the column
        Pcr = np.pi**2 * params['E'] * params['I'] / L**2    # Euler buckling load

        # ========== setting analysis parameters ==============
        target_load_level = 0.99      # 99% of Euler load
        max_steps = 10                # solve max_steps points on the primary path
        #w   *= 0.01
        #Pcr *= 0.01
        # define a list of target load levels
        load_levels = np.linspace(0, target_load_level, max_steps)

        #
        # ==== Build the system model ====
        #
        model = System()
        model.setSolver(NewtonRaphsonSolver())
        # create nodes
        y = 0.0
        nd0 = Node(0.0, y)
        model += nd0

        ndi = nd0
        #params['E']
        for i in range(N):
            # nodes
            y += (1-alpha)/2*L/N
            ndj = Node(0.0, y)
            model += ndj
            # elements
            elem = Frame2D(ndi, ndj, ElasticSection(params))
            model += elem
            # ** apply the element portion of the reference load
            elem.setDistLoad(w)
            ndi = ndj    # jump to next element: make current end-node the next start-node

        #params['E']
        for i in range(N):
            # nodes
            y += alpha*L/N
            ndj = Node(0.0, y)
            model += ndj
            # elements
            elem = Frame2D(ndi, ndj, ElasticSection(stiffer_params))
            model += elem
            # ** apply the element portion of the reference load
            elem.setDistLoad(w)
            ndi = ndj    # jump to next element: make current end-node the next start-node

        #params['E']
        for i in range(N):
            # nodes
            y += (1-alpha)/2*L/N
            ndj = Node(0.0, y)
            model += ndj
            # elements
            elem = Frame2D(ndi, ndj, ElasticSection(params))
            model += elem
            # ** apply the element portion of the reference load
            elem.setDistLoad(w)
            ndi = ndj    # jump to next element: make current end-node the next start-node

        # define support(s)
        nd0.fixDOF('ux', 'uy')      # horizontal support left end
        ndi.fixDOF('ux', )          # vertical support right end
        # ==== complete the reference load ====
        # these are only nodal forces as part of the reference load
        # .. load only the upper node
        ndi.setLoad((-Pcr,), ('uy',))

        print(model)

        model.initRecorder()     # sets variables to track; defaults are load_level and stability index

        model.startRecorder()    # this starts the recording of requested variables
        model.trackStability(True)

        for lam in np.linspace(0.0, 0.9, 5):
            model.setLoadFactor(lam)
            model.solve(verbose=True)
            model.recordThisStep()

        model.plot(factor=100.00)

        model.beamValuePlot('F')
        model.beamValuePlot('V')
        model.beamValuePlot('M')

        model.historyPlot('stability')


