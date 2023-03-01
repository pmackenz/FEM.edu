
import sys, os
sys.path.insert(0, os.path.abspath("."))

from femedu.examples.trusses.truss01 import *
from femedu.examples.trusses.truss02 import *
from femedu.examples.trusses.truss03 import *
from femedu.examples.trusses.truss04 import *
from femedu.examples.beams.beam01 import *
from femedu.examples.beams.beam02 import *
from femedu.examples.frames.frame01 import *
from femedu.examples.frames.frame02 import *
from femedu.examples.frames.frame03 import *
from femedu.examples.frames.frame04 import *


if __name__ == "__main__":

    #ex = ExampleTruss01()
    #ex = ExampleTruss02()
    #ex = ExampleTruss03()
    #ex = ExampleTruss04()

    #ex = ExampleBeam01()
    #ex = ExampleBeam02()

    #ex = ExampleFrame01()
    #ex = ExampleFrame02()
    ex = ExampleFrame03()
    #ex = ExampleFrame04()

    # print the doc-string for the current example
    print(ex)

    # run the actual problem
    ex.run()

