import numpy as np
import matplotlib.pyplot as plt


class Material():
    """
    abstract class: representing a generic material

    """

    """
    defining class parameters
    
    Add a new one for every general type of material.
    This will be used by elements to test for compatible materials
    """
    UNKNOWN     = 0x000000
    DEFAULT     = 0x000001
    FIBER       = 0x000002
    SECTION1D   = 0x000004  # integrated section (beams and frames)
    SECTION2D   = 0x000008  # through the thickness integration (plates and shells)
    PLANESTRESS = 0x000010
    PLANESTRAIN = 0x000020
    CONTINUUM   = 0x000040

    ELASTIC     = 0x010000
    PLASTIC     = 0x020000
    CREEP       = 0x040000
    HARDENING   = 0x080000


    def __init__(self, params={'E':1.0, 'A':1.0, 'nu':0.0, 'fy':1.0e30}):
        """

        :param params:
        """
        self._type = self.UNKNOWN

        self.parameters = params

        # make sure all necessary parameters exist
        if 'E' not in self.parameters:
            self.parameters['E']  = 1.0
        if 'A' not in self.parameters:
            self.parameters['A']  = 1.0
        if 'nu' not in self.parameters:
            self.parameters['nu'] = 0.0
        if 'fy' not in self.parameters:
            self.parameters['fy'] = 1.0e30

        self.sig = 0.0
        self.Et  = self.parameters['E']

    def __str__(self):
        s = "{}(Material)({})".format(self.__class__.__name__, self.parameters)
        return s

    def __repr__(self):
        return str(self)

    def materialType(self):
        return self._type

    def getStress(self):
        """
        request axial stress

        :return: sigma
        """
        return self.stress

    def getStiffness(self):
        """
        request axial stiffness

        :return: Et ... tangent material stiffness
        """
        return self.Et

    def setStrain(self, eps):
        """
        update state for a user provided axial strain value

        :param eps:  strain or strain tensor
        """
        # update stress state
        self.strain = eps
        self.updateState()

    def updateState(self):
        raise NotImplementedError(self.__class__.__name__ + '.updateState() needs to be overloaded')

    def getStrain(self):
        return self.strain

    def converged(self):
        raise NotImplementedError(self.__class__.__name__ + '.converged() needs to be overloaded')


if __name__ == "__main__":
    # testing the Node class
    mat = Material(params={'E':100.0, 'nu':0.0, 'fy':1.0})

    eps = 0.02 * np.sin( np.linspace(0, 2.*np.pi, 100) )

    sig = []
    Et = []
    for strain in eps:
        # update material strain state
        mat.setStrain(strain)
        #collect stress response
        sig.append(mat.getStress())
        #collect material tangent stiffness response
        Et.append(mat.getStiffness())

    fig, (ax1,ax2) = plt.subplots(2,1)

    ax1.plot(eps, sig, '-r', label='stress')
    ax1.grid(True)
    ax1.set_xlabel('strain $\\varepsilon$')
    ax1.set_ylabel('stress $\sigma$')
    ax1.legend()

    ax2.plot(eps, Et, '-.b', label='tangent modulus')
    ax2.grid(True)
    ax2.set_xlabel('strain $\\varepsilon$')
    ax2.set_ylabel('tangent modulus $E_t$')
    ax2.legend()

    plt.show()
