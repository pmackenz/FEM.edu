from .Recorder import *
from .Record import *

class ElementRecorder(Recorder):
    """

    """

    def __init__(self, **kwargs):
        super(ElementRecorder, self).__init__(**kwargs)

        if 'label' in kwargs:
            template = kwargs['label'] + ":{}"
        else:
            template = "elem:{}"

        for var in self.elemvars:
            lbl = template.format(var)
            if var == 'stress' or var == 'strain':
                for key in ('xx','yy','zz','xy','yz','zx',):
                    lbl += f":{key}"
                    reference = f"{var}:{key}"
                    self.data[reference] = Record(key=reference, label=lbl)
            else:
                self.data[var] = Record(key=var, label=lbl)

            print(lbl)

