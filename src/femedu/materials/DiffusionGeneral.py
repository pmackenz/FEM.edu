from .Material import Material
import numpy as np

class DiffusionGeneral(Material):
    """
    Generic diffusion model for 2D and 3D diffusion problems.

    Use subclasses :py:class:`Thermal`, :py:class:`Seapage`, etc., for actual analyses.
    """
    
    def __init__(self, params={'diffusivity':1., 'capacity':1., 'density':1., 'thickness':1.}):
        super(DiffusionGeneral, self).__init__(params)

        self._type = self.DIFFUSION

        if 'diffusivity' not in self.parameters:
            self.parameters['diffusivity'] = 1.0
        if 'capacity' not in self.parameters:
            self.parameters['capacity'] = 1.0
        if 'density' not in self.parameters:
            self.parameters['density'] = 1.0
        if 'thickness' not in self.parameters:
            self.parameters['thickness'] = 1.0

        self.gradPhi = None
        self.flux    = None

    def setGrad(self, gradPhi):
        """
        :param gradPhi: gradient :math:`\\nabla \\phi` of the scalar potential :math:`\\phi`
        :type gradPhi: numpy.ndarray
        """
        self.gradPhi = gradPhi
        self.updateState()

    def getGrad(self):
        """
        :returns gradPhi: gradient :math:`\\nabla \\phi` of the scalar potential :math:`\\phi`
        :type gradPhi: numpy.ndarray
        """
        if self.gradPhi:
            return self.gradPhi
        else:
            return np.zeros(3)

    def getFlux(self):
        """
        :return: flux (numpy.ndarray)
        """
        if self.flux:
            return self.flux
        else:
            return np.zeros(3)

    def getThickness(self):
        """
        Returns thickness of a 2D plate, or 1.0 if none given or 3D problem
        """
        return self.parameters['thickness']

    def getCapacity(self):
        """
        returns :math:`C := \\varrho c`
        """
        return self.parameters['capacity']*self.parameters['density']

    def getDiffusivity(self):
        """
        Returns diffusivity of a 2D plate, or 1.0 if none given
        """
        return self.parameters['diffusivity']

    def updateState(self):
        """
        Update the material state after :math:`\\nabla \\phi`
        """
        if self.gradPhi:
            self.flux = -self.parameters['diffusivity'] * self.gradPhi

