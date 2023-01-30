
import sys, os
sys.path.insert(0, os.path.abspath("."))

from domain.System import *
from domain.Node import *
from elements.Truss import *
from elements.LinearTriangle import *
from materials.FiberMaterial import *
from materials.PlaneStress import *

#from domain import *


def problem1():
    # initialize a system model
    B = 6.0 * 12
    H = 3.0 * 12
    params = {'E': 10., 'A': 1., 'nu': 0.0, 'fy': 1.e30}

    model = System()

    # create nodes
    nd0 = Node(0.0, 0.0)
    nd1 = Node(  B, 0.0)
    nd2 = Node(0.5*B, H)

    model.addNode(nd0)
    model.addNode(nd1)
    model.addNode(nd2)

    # create elements
    model.addElement(Truss(nd0, nd1, FiberMaterial(params)))  # bottom 1
    model.addElement(Truss(nd0, nd2, FiberMaterial(params)))  # up right diag 1
    model.addElement(Truss(nd1, nd2, FiberMaterial(params)))  # up left diag 1

    # define support(s)
    nd0.fixDOF(0)    # horizontal support left end
    nd0.fixDOF(1)    # vertical support left end
    nd1.fixDOF(1)    # vertical support right end

    # add loads
    # .. load only the upper nodes
    nd2.setLoad(0.0, -1.0)

    # analyze the model
    model.solve()

    # write out report
    model.report()

    # create plots
    model.plot(factor=1.)

def problem2():
    # initialize a system model
    B = 6.0 * 12
    H = 8.0 * 12
    params = {'E': 10000., 'A': 3., 'nu': 0.0, 'fy': 1.e30}

    model = System()

    # create nodes
    nd0 = Node(0.0, 0.0)
    nd1 = Node(  B, 0.0)
    nd2 = Node(2*B, 0.0)
    nd3 = Node(3*B, 0.0)
    nd4 = Node(4*B, 0.0)
    nd5 = Node(0.5*B, H)
    nd6 = Node(1.5*B, H)
    nd7 = Node(2.5*B, H)
    nd8 = Node(3.5*B, H)

    model.addNode(nd0)
    model.addNode(nd1)
    model.addNode(nd2)
    model.addNode(nd3)
    model.addNode(nd4)
    model.addNode(nd5)
    model.addNode(nd6)
    model.addNode(nd7)
    model.addNode(nd8)

    # create elements
    model.addElement(Truss(nd0, nd1, FiberMaterial(params)))  # bottom 1
    model.addElement(Truss(nd1, nd2, FiberMaterial(params)))  # bottom 2
    model.addElement(Truss(nd2, nd3, FiberMaterial(params)))  # bottom 3
    model.addElement(Truss(nd3, nd4, FiberMaterial(params)))  # bottom 4

    model.addElement(Truss(nd5, nd6, FiberMaterial(params)))  # upper 1
    model.addElement(Truss(nd6, nd7, FiberMaterial(params)))  # upper 2
    model.addElement(Truss(nd7, nd8, FiberMaterial(params)))  # upper 3

    model.addElement(Truss(nd0, nd5, FiberMaterial(params)))  # up right diag 1
    model.addElement(Truss(nd1, nd6, FiberMaterial(params)))  # up right diag 2
    model.addElement(Truss(nd2, nd7, FiberMaterial(params)))  # up right diag 3
    model.addElement(Truss(nd3, nd8, FiberMaterial(params)))  # up right diag 4

    model.addElement(Truss(nd1, nd5, FiberMaterial(params)))  # up left diag 1
    model.addElement(Truss(nd2, nd6, FiberMaterial(params)))  # up left diag 2
    model.addElement(Truss(nd3, nd7, FiberMaterial(params)))  # up left diag 3
    model.addElement(Truss(nd4, nd8, FiberMaterial(params)))  # up left diag 4

    # define support(s)
    nd0.fixDOF(0)    # horizontal support left end
    nd0.fixDOF(1)    # vertical support left end
    nd4.fixDOF(1)    # vertical support right end

    # add loads
    # .. load only the upper nodes
    nd5.setLoad(0.0, -1.0)
    nd6.setLoad(0.0, -1.0)
    nd7.setLoad(0.0, -1.0)
    nd8.setLoad(0.0, -1.0)

    # analyze the model
    model.solve()

    # write out report
    model.report()

    # create plots
    model.plot(factor=250.)

def problem3():
    # initialize a system model
    B = 6.0 * 12
    H = 3.0 * 12
    params = {'E': 10., 'A': 1., 'nu': 0.0, 'fy': 1.e30}

    model = System()

    # create nodes
    nd0 = Node(0.0, 0.0)
    nd1 = Node(  B, 0.0)
    nd2 = Node(0.5*B, H)

    model += nd0
    model += nd1
    model += nd2

    # create elements
    model += Truss(nd0, nd1, FiberMaterial(params))  # bottom 1
    model += Truss(nd0, nd2, FiberMaterial(params))  # up right diag 1
    model += Truss(nd1, nd2, FiberMaterial(params))  # up left diag 1

    # define support(s)
    #nd0.fixDOF(0)    # horizontal support left end
    nd0 //= 0
    nd0.fixDOF(1)    # vertical support left end
    nd1.fixDOF(1)    # vertical support right end

    # add loads
    # .. load only the upper nodes
    nd2.setLoad(0.0, -1.0)

    # analyze the model
    model.solve()

    # write out report
    model.report()

    # create plots
    model.plot(factor=1.)

    # fix horizontal motion of node 1
    nd1.fixDOF(0)

    # add loads: same load -- nothing to do

    # RE-analyze the model
    model.resetDisp()
    model.solve()

    # skip the report
    model.report()

    # create plots
    model.plot(factor=2.)

def problem4():

    params = {'E': 10., 'A': 1., 'nu': 0.0, 'fy': 1.e30}

    nd0 = Node(  0.0,  0.0)
    nd1 = Node( 10.0, 10.0)
    nd2 = Node(  0.0, 20.0)

    mat = PlaneStress(params)

    elem1 = LinearTriangle(nd0, nd1, nd2, mat)

    nd0.setDisp( 0.0, 0.0)
    nd1.setDisp( 5.0, 5.0)
    nd2.setDisp( 0.0,-5.0)

    elem1.updateState()

    print(elem1)


if __name__ == "__main__":

    #problem1()
    problem2()
    #problem3()
    problem4()

