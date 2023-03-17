
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
from femedu.examples.plates.plate01 import *
from femedu.examples.plates.plate02 import *
from femedu.examples.plates.plate03 import *
from femedu.examples.plates.plate04 import *
from femedu.examples.course_help.final01 import *


if __name__ == "__main__":

    examples = (
        ExampleTruss01(),
        ExampleTruss02(),
        ExampleTruss03(),
        ExampleTruss04(),
        #
        ExampleBeam01(),
        ExampleBeam02(),
        #
        ExampleFrame01(),
        ExampleFrame02(),
        ExampleFrame03(),
        ExampleFrame04(),
        #
        ExamplePlate01(),
        ExamplePlate02(),
        ExamplePlate03(),
        ExamplePlate04()
    )




    #ex = ExampleFinal01()

    for ex in examples:

        # print the doc-string for the current example
        print(ex)

        # run the actual problem
        ex.run()

