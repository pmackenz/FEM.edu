from .Recorder import *
from .Record import *

class MaterialRecorder(Recorder):
    """
    Recorder class for state variables at integration points

    .. warning::

        Use this sparely and only for very few elements since is can generate
        **lots of mostly unnecessary data**.

    """

    def __init__(self, **kwargs):
        super(MaterialRecorder, self).__init__(**kwargs)

        if 'label' in kwargs:
            template = kwargs['label'] + ":{}"
        else:
            template = "material:{}"

        for var in self.materialvars:
            lbl = template.format(var)
            self.data[var] = Record(key=var, label=lbl)

