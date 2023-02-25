
import sys, os
sys.path.insert(0, os.path.abspath("."))

from examples.trusses.truss01 import *
from examples.trusses.truss02 import *
from examples.trusses.truss03 import *
from examples.trusses.truss04 import *
from examples.beams.beam01 import *


if __name__ == "__main__":

    #ex = ExampleTruss01()
    #ex = ExampleTruss02()
    #ex = ExampleTruss03()
    #ex = ExampleTruss04()

    ex = ExampleBeam01()

    # print the doc-string for the current example
    print(ex)

    # run the actual problem
    ex.run()

