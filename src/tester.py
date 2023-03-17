from femedu.domain import *
from femedu.elements.Frame2D import *
from femedu.materials.ElasticSection import *
from femedu.mesher.CurveMesher import *

if __name__ == "__main__":

    model = System()
    mesher = CurveMesher(model, (0,0),(1.5,.25),(2,1),(3.,1.5))
    nodes, elements = mesher.mesh(10, Frame2D, ElasticSection())

    model.report()

    #mesher._create_doc()
