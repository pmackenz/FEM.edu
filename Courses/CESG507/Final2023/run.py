import matplotlib.pyplot as plt
import numpy as np

from setup import *

from femedu.domain import *
from femedu.solver.NewtonRaphsonSolver import *
from femedu.elements.Frame2D import *
from femedu.materials.ElasticSection import *
from femedu.examples.Example import *


class ExampleFinal01(Example):

    def docString(self):
        s = """
        CESG 507 - 2023 - Final Exam Problem #1
        """
        return s

    def problem(self):

        # self.debugProblem()
        self.fullProblem()

    def debugProblem(self):

        model = self.model1(alpha=0.80, gamma=1.0, bottom='a', top='ii')
        self.doAnalysis(model, target_load_level=1.2, num_steps=3, verbose=True)
        model.historyPlot('stability')

    def fullProblem(self):

        support_options = (
            {'bottom':'a', 'top': 'ii', 'target':1.1},
            {'bottom':'a', 'top':'iii', 'target':2.5},
            {'bottom':'b', 'top':  'i', 'target':0.55},
            {'bottom':'b', 'top': 'ii', 'target':1.1},
            {'bottom':'b', 'top':'iii', 'target':4.5}
        )

        supp_map = {
            'a':'pin',
            'b':'fix',
            'i':'free',
            'ii':'pin',
            'iii':'fix',
        }

        alphas     = np.linspace(0.0,1.0,41)
        gammas     = [1.0, 1.5, 2.0]
        linestyles = [':','-','-.']

        for support in support_options:

            # initiate results plot
            plt.figure( figsize=(8.0, 6.0) )

            for gam, ls in zip(gammas, linestyles):

                data = {
                        'key':'{gamma:.1f}_{bottom}_{top}'.format(
                            gamma=gam,
                            bottom=supp_map[support['bottom']],
                            top=supp_map[support['top']]
                        ),
                        'model1':[],
                        'model2':[],
                        'model3':[]
                }

                for alph in alphas:

                    model1 = self.model1(alpha=alph, gamma=gam, **support)
                    model2 = self.model2(alpha=alph, gamma=gam, **support)
                    model3 = self.model3(alpha=alph, gamma=gam, **support)

                    llvl = support['target']

                    data['model1'].append( self.doAnalysis(model1, target_load_level=llvl, num_steps=3) )
                    data['model2'].append( self.doAnalysis(model2, target_load_level=llvl, num_steps=3) )
                    data['model3'].append( self.doAnalysis(model3, target_load_level=llvl, num_steps=3) )

                # create the stability plot
                txt = "system I_{}".format( data['key'] )
                plt.plot(alphas, data['model1'], linestyle=ls, color='r', label=txt)
                txt = "system II_{}".format( data['key'] )
                plt.plot(alphas, data['model2'], linestyle=ls, color='b', label=txt)
                txt = "system III_{}".format( data['key'] )
                plt.plot(alphas, data['model3'], linestyle=ls, color='g', label=txt)

            plt.xlabel('$\\alpha$')
            plt.ylabel('$P_{cr}~/~(\pi^2 EI / L^2)$')
            plt.grid(True)
            plt.legend()
            plt.savefig("study_{}.png".format(data['key']))
            #plt.show()

    def model1(self, alpha=0.33333333333, gamma=1.0, bottom='a', top='ii', **kwargs):

        # ========== setting mesh parameters ==============
        N = 10                 # number of elements in the mesh
        L = 100.0              # column free length

        # ========== setting material parameters ==============
        params_A = dict(
            E = 20000.,       # Young's modulus
            A = 100.0,        # cross section area
            I = 10.0          # cross section moment of inertia
        )

        params_B = dict(
            E = 20000.,       # Young's modulus
            A = 100.0,        # cross section area
            I = gamma * 10.0  # cross section moment of inertia
        )

        # ========== setting load parameters ==============
        w   = 0.0        # uniform lateral load on the column
        Pcr = np.pi**2 * params_A['E'] * params_A['I'] / L**2    # Euler buckling load

        # ========== setting analysis parameters ==============
        target_load_level = 0.99      # 99% of Euler load
        max_steps = 10                # solve max_steps points on the primary path

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
        #bottom section
        Le = 0.5*(1.0 - alpha) * L/N
        if Le > 0.0:
            for i in range(N):
                # nodes
                y += Le
                ndj = Node(0.0, y)
                model += ndj
                # elements
                elem = Frame2D(ndi, ndj, ElasticSection(params_A))
                model += elem
                # ** apply the element portion of the reference load
                elem.setDistLoad(w)
                ndi = ndj    # jump to next element: make current end-node the next start-node

        # center section
        Le = alpha * L/N
        if Le > 0.0:
            for i in range(N):
                # nodes
                y += Le
                ndj = Node(0.0, y)
                model += ndj
                # elements
                elem = Frame2D(ndi, ndj, ElasticSection(params_B))
                model += elem
                # ** apply the element portion of the reference load
                elem.setDistLoad(w)
                ndi = ndj    # jump to next element: make current end-node the next start-node

        # top section
        Le = 0.5*(1.0 - alpha) * L/N
        if Le > 0.0:
            for i in range(N):
                # nodes
                y += Le
                ndj = Node(0.0, y)
                model += ndj
                # elements
                elem = Frame2D(ndi, ndj, ElasticSection(params_A))
                model += elem
                # ** apply the element portion of the reference load
                elem.setDistLoad(w)
                ndi = ndj    # jump to next element: make current end-node the next start-node

        # define support(s)
        if bottom == 'a':
            nd0.fixDOF('ux', 'uy')        # pin support at bottom
        else:
            nd0.fixDOF('ux', 'uy', 'rz')  # fixed support at bottom

        if top == 'ii':
            ndi.fixDOF('ux', )            # horizontal support right end
        elif top == 'iii':
            ndi.fixDOF('ux', 'rz')        # horizontal support right end, fix rotation

        # ==== complete the reference load ====

        # these are only nodal forces as part of the reference load
        # .. load only the upper node
        ndi.setLoad((-Pcr,), ('uy',))

        return model

    def model2(self, alpha=0.33333333333, gamma=1.0, bottom='a', top='ii', **kwargs):

        # ========== setting mesh parameters ==============
        N = 10                 # number of elements in the mesh
        L = 100.0              # column free length

        # ========== setting material parameters ==============
        params_A = dict(
            E = 20000.,       # Young's modulus
            A = 100.0,        # cross section area
            I = gamma * 10.0  # cross section moment of inertia
        )

        params_B = dict(
            E = 20000.,       # Young's modulus
            A = 100.0,        # cross section area
            I = 10.0          # cross section moment of inertia
        )

        # ========== setting load parameters ==============
        w   = 0.0        # uniform lateral load on the column
        Pcr = np.pi**2 * params_B['E'] * params_B['I'] / L**2    # Euler buckling load

        # ========== setting analysis parameters ==============
        target_load_level = 0.99      # 99% of Euler load
        max_steps = 10                # solve max_steps points on the primary path

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
        #bottom section
        Le = 0.5 * alpha * L/N
        if Le > 0.0:
            for i in range(N):
                # nodes
                y += Le
                ndj = Node(0.0, y)
                model += ndj
                # elements
                elem = Frame2D(ndi, ndj, ElasticSection(params_A))
                model += elem
                # ** apply the element portion of the reference load
                elem.setDistLoad(w)
                ndi = ndj    # jump to next element: make current end-node the next start-node

        # center section
        Le = (1.0 - alpha) * L/N
        if Le > 0.0:
            for i in range(N):
                # nodes
                y += Le
                ndj = Node(0.0, y)
                model += ndj
                # elements
                elem = Frame2D(ndi, ndj, ElasticSection(params_B))
                model += elem
                # ** apply the element portion of the reference load
                elem.setDistLoad(w)
                ndi = ndj    # jump to next element: make current end-node the next start-node

        # top section
        Le = 0.5 * alpha * L/N
        if Le > 0.0:
            for i in range(N):
                # nodes
                y += Le
                ndj = Node(0.0, y)
                model += ndj
                # elements
                elem = Frame2D(ndi, ndj, ElasticSection(params_A))
                model += elem
                # ** apply the element portion of the reference load
                elem.setDistLoad(w)
                ndi = ndj    # jump to next element: make current end-node the next start-node

        # define support(s)
        if bottom == 'a':
            nd0.fixDOF('ux', 'uy')        # pin support at bottom
        else:
            nd0.fixDOF('ux', 'uy', 'rz')  # fixed support at bottom

        if top == 'ii':
            ndi.fixDOF('ux', )            # horizontal support right end
        elif top == 'iii':
            ndi.fixDOF('ux', 'rz')        # horizontal support right end, fix rotation

        # ==== complete the reference load ====

        # these are only nodal forces as part of the reference load
        # .. load only the upper node
        ndi.setLoad((-Pcr,), ('uy',))

        return model

    def model3(self, alpha=0.33333333333, gamma=1.0, bottom='a', top='ii', **kwargs):

        # ========== setting mesh parameters ==============
        N = 10                 # number of elements in the mesh
        L = 100.0              # column free length

        # ========== setting material parameters ==============
        params_A = dict(
            E = 20000.,       # Young's modulus
            A = 100.0,        # cross section area
            I = gamma * 10.0  # cross section moment of inertia
        )

        params_B = dict(
            E = 20000.,       # Young's modulus
            A = 100.0,        # cross section area
            I = 10.0          # cross section moment of inertia
        )

        # ========== setting load parameters ==============
        w   = 0.0        # uniform lateral load on the column
        Pcr = np.pi**2 * params_B['E'] * params_B['I'] / L**2    # Euler buckling load

        # ========== setting analysis parameters ==============
        target_load_level = 0.99      # 99% of Euler load
        max_steps = 10                # solve max_steps points on the primary path
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
        #bottom section
        Le = alpha * L/N
        if Le > 0.0:
            for i in range(N):
                # nodes
                y += Le
                ndj = Node(0.0, y)
                model += ndj
                # elements
                elem = Frame2D(ndi, ndj, ElasticSection(params_A))
                model += elem
                # ** apply the element portion of the reference load
                elem.setDistLoad(w)
                ndi = ndj    # jump to next element: make current end-node the next start-node

        # top section
        Le = (1.0 - alpha) * L/N
        if Le > 0.0:
            for i in range(N):
                # nodes
                y += Le
                ndj = Node(0.0, y)
                model += ndj
                # elements
                elem = Frame2D(ndi, ndj, ElasticSection(params_B))
                model += elem
                # ** apply the element portion of the reference load
                elem.setDistLoad(w)
                ndi = ndj    # jump to next element: make current end-node the next start-node

        # define support(s)
        if bottom == 'a':
            nd0.fixDOF('ux', 'uy')        # pin support at bottom
        else:
            nd0.fixDOF('ux', 'uy', 'rz')  # fixed support at bottom

        if top == 'ii':
            ndi.fixDOF('ux', )            # horizontal support right end
        elif top == 'iii':
            ndi.fixDOF('ux', 'rz')        # horizontal support right end, fix rotation

        # ==== complete the reference load ====

        # these are only nodal forces as part of the reference load
        # .. load only the upper node
        ndi.setLoad((-Pcr,), ('uy',))

        return model

    def doAnalysis(self, model, target_load_level=1.0, num_steps=3, verbose=False):

        #print(model)

        model.initRecorder()     # sets variables to track; defaults are load_level and stability index

        model.startRecorder()    # this starts the recording of requested variables
        model.trackStability(True)

        for lam in np.linspace(0.0, target_load_level, num_steps):
            model.setLoadFactor(lam)
            model.solve(verbose=False)
            model.solve(verbose=True)
            if verbose:
                print(model.checkStability(num_eigen=5))
            model.recordThisStep()

        #model.plot(factor=100.00)

        #model.beamValuePlot('F')
        #model.beamValuePlot('V')
        #model.beamValuePlot('M')

        #model.historyPlot('stability')
        #model.plotBucklingMode(factor=25.0)

        # estimate critical load from recorded data
        if 'lam' in model.recorder.data and len(model.recorder.data['lam']) >= 2:
            lam   = model.recorder.data['lam'][-2:]
            detKt = model.recorder.data['stability'][-2:]

            dl =   lam[1]    -   lam[0]
            dk = detKt[1][0] - detKt[0][0]

            dkdl = dk/dl

            # k = detKt[0] + dkdl * (ll - lam[0])
            ll = lam[0] - detKt[0] / dkdl
            return ll
        else:
            return 0.0


if __name__ == "__main__":

    ex = ExampleFinal01()
    ex.run()
