from .Material import *


class SectionMaterial(Material):
    """
    class: representing an integrated Section Material

    """

    def __init__(self, params={'E':1.0, 'A':1.0, 'I':1.0, 'nu':0.0, 'fy':1.0e30}):
        super().__init__(params = params)

        self._type = self.SECTION1D

        # make sure all necessary parameters exist
        if 'E' not in self.parameters:
            self.parameters['E']  = 1.0
        if 'A' not in self.parameters:
            self.parameters['A']  = 1.0
        if 'I' not in self.parameters:
            self.parameters['I']  = 1.0
        if 'nu' not in self.parameters:
            self.parameters['nu'] = 0.0
        if 'fy' not in self.parameters:
            self.parameters['fy'] = 1.0e30

        # initialize strain
        self.plastic_strain = 0.0
        self.setStrain({'axial':0.0, 'flexure':0.0})
        self.updateState()

    def updateState(self):
        EA = self.parameters['E']*self.parameters['A']
        EI = self.parameters['E']*self.parameters['I']

        if 'axial' in self.strain:
            force = EA * self.strain['axial']
        else:
            force = 0.0

        if 'flexure' in self.strain:
            moment = EI * self.strain['flexure']
        else:
            moment = 0.0

        self.stress = {'axial':force, 'flexure':moment}
        self.Et = {'ax-ax':EA, 'flx-flx':EI, 'ax-flx':0.0, 'flx-ax':0.0}

