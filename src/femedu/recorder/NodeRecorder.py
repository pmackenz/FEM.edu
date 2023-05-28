from .Recorder import *
from .Record import *

class NodeRecorder(Recorder):
    """

    """

    def __init__(self, **kwargs):
        super(NodeRecorder, self).__init__(**kwargs)

        if 'label' in kwargs:
            template = kwargs['label'] + ":{}"
        else:
            template = "node:{}"

        for var in self.nodevars:
            lbl = template.format(var)
            self.data[var] = Record(key=var, label=lbl)

