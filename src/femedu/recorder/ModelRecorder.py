from .Recorder import Recorder
from .NodeRecorder import *
from .ElementRecorder import *

class ModelRecorder(Recorder):
    """

    """

    def __init__(self,**kwargs):
        super(ModelRecorder, self).__init__(**kwargs)

        for var in self.sysvars:
            self.data[var] = Record(key=var, label=f"sys:{var}")

        self.nodes    = []
        self.elements = []

        if 'nodes' in kwargs:
            self.nodes = kwargs['nodes']
            for node in self.nodes:
                node.setRecorder(NodeRecorder(variables=self.nodevars, label=node.getID(), type='node'))

        if 'elements' in kwargs:
            self.elements = kwargs['elements']
            for elem in self.elements:
                elem.setRecorder(ElementRecorder(variables=self.elemvars, type='element'))

