import numpy as np

from .Material import *

class VonMises(Material):
    """
    class: representing a 3d von Mises Material

    Assumptions are isotropic linear elastic-hardening plastic behavior.

    """

    def __init__(self, params={'E':1.0, 'nu':0.0, 'fy':1.0e30}):
        super().__init__(params = params)

        self._type = self.PLANESTRESS | self.PLASTIC | self.HARDENING

        # isotropic hardening parameter
        if 'H' not in self.parameters:
            self.parameters['H'] = 0.0

        # kinematic hardening parameter
        if 'K' not in self.parameters:
            self.parameters['K'] = 0.0

        # initialize strain
        self.plastic_strain = np.zeros(6)
        self.alpha          = 0.0
        self.beta           = np.zeros(6)

        self.setStrain({'xx':0.0, 'yy':0.0, 'zz':0.0, 'xy':0.0, 'yz':0.0, 'zx':0.0})

    def getThickness(self):
        raise NotImplementedError('Thickness has no meaning in a 3d model')

    def updateState(self, eps):
        """
        update state for a user provided axial strain value

        :param eps:  strain or strain tensor
        :return: n/a
        """

        # update stress state
        E  = self.parameters['E']
        nu = self.parameters['nu']
        fy = self.parameters['fy']
        H  = self.parameters['H']   # kinematic hardening parameter
        K  = self.parameters['K']   # isotropic hardening parameter
        twoG  = E / (1 + nu)

        Idev = np.array([ [2./3., -1./3., -1./3.,  0.,  0.,  0.],
                          [2./3., -1./3., -1./3.,  0.,  0.,  0.],
                          [2./3., -1./3., -1./3.,  0.,  0.,  0.],
                          [   0.,     0.,     0., 0.5,  0.,  0.],
                          [   0.,     0.,     0.,  0., 0.5,  0.],
                          [   0.,     0.,     0.,  0.,  0., 0.5] ])
        Ce_inv = (1 / E) * np.array([ [1., -nu, -nu, 0., 0., 0.],
                                    [-nu, 1., -nu, 0., 0., 0.],
                                    [-nu, -nu, 1., 0., 0., 0.],
                                    [ 0., 0., 0.,  2. * (1. + nu), 0., 0.],
                                    [ 0., 0., 0.,  0., 2. * (1. + nu), 0.],
                                    [ 0., 0., 0.,  0., 0., 2. * (1. + nu)] ])
        Ce = np.linalg.inv( Ce_inv )

        # epsilon_{n+1}
        eps1 = np.array([self.strain['xx'],self.strain['yy'],self.strain['zz'],
                        self.strain['xy'],self.strain['yz'],self.strain['zx']])

        # set trial state
        self.plastic_strain1 = self.plastic_strain.copy()
        self.alpha1          = self.alpha
        self.beta1           = self.beta.copy()

        self.Et = Ce

        # elastic predictor (trial state)
        self.sig = self.Et @ ( eps1 - self.plastic_strain )
        q_iso = - np.sqrt(2./3.) * K * self.alpha1
        q_kin = - H * self.beta1

        # check yield condition
        one = np.array( [1.,1.,1.,0.,0.,0.] )
        I1 = self.sig @ one
        eta = self.sig - (I1/3.) * one + q_kin
        norm_eta = np.sqrt( eta @ eta )

        f = norm_eta + np.sqrt(2./3.) * (q_iso - fy)

        # plastic corrector as needed
        if f >= 0.0:
            # loop
            r = eta / norm_eta
            Cxr = twoG * r  # Ce . r = 2G r   since tr(r) == 0
            gamma = f / (twoG + K + (2/3)*H)

            # update state variables
            self.alpha1 += np.sqrt(2/3) * gamma
            self.beta1  += gamma * r

            rstrain = r.copy()
            rstrain[3] *= 2  # gamma_xy^p = 2 eps_xy^p
            rstrain[4] *= 2  # gamma_yz^p = 2 eps_yz^p
            rstrain[5] *= 2  # gamma_zx^p = 2 eps_zx^p
            self.plastic_strain1 += gamma * rstrain

            self.sig = Ce @ ( self.strain - self.plastic_strain1 )

            Cep  = Ce - twoG*twoG*gamma/norm_eta * Idev
            Cep -= twoG*twoG/(twoG + K + (2/3)*H) * (1 - f/norm_eta) * np.outer(r,r)

            self.Et  = Cep

    def getStrain(self):
        return self.sig / self.parameters['E'] + self.plastic_strain

    def converged(self):
        # update state now that the global analysis has converged
        self.plastic_strain = self.plastic_strain1
        self.alpha          = self.alpha1
        self.beta           = self.beta1


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    # testing the Node class
    mat = VonMises(params={'E':100.0, 'nu':0.0, 'fy':1.0})

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
