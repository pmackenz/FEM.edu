from .Material import *

class PlaneStress(Material):
    """
    class: representing a 2d Plane Stress Material

    """

    def __init__(self, params={'E':1.0, 't':1.0, 'nu':0.0, 'fy':1.0e30}):
        super().__init__(params = params)

        self._type = self.PLANESTRESS | self.PLASTIC

        if 't' not in self.parameters:
            self.parameters['t'] = 1.0

        # initialize strain
        self.total_strain   = np.zeros(3)
        self.plastic_strain = np.zeros(3)
        self.setStrain({'xx':0.0, 'yy':0.0, 'zz':0.0, 'xy':0.0, 'yz':0.0, 'zx':0.0})

    def getThickness(self):
        return self.parameters['t']

    def updateState(self):
        """
        update state for a user provided axial strain value

        :param eps:  strain or strain tensor
        :return: n/a
        """
        eps = np.array([self.strain['xx'], self.strain['yy'], self.strain['xy']])

        # update stress state
        E  = self.parameters['E']
        t  = self.parameters['t']
        nu = self.parameters['nu']
        fy = self.parameters['fy']

        # default consistency parameter
        gamma = 0.0

        self.Et = E*t/(1. - nu*nu) * np.array([[1.,nu,0.],[nu,1.,0.],[0.,0.,(1.-nu)/2.]])

        Cinv = 1/(E*t) * np.array([[1.,-nu,0.],[-nu,1.,0.],[0.,0.,2.*(1.+nu)]])
        Phi  = np.array([[2.,-1.,0.],[-1.,2.,0.],[0.,0.,6.]])

        # elastic predictor
        stress = self.Et @ ( eps - self.plastic_strain )

        self.strain['zz'] = -self.parameters['nu'] * (eps[0] + eps[1])  # elastic thickness strain

        # check yield condition
        (sxx, syy, sxy) = stress
        f = stress @ Phi @ stress / 2. - (t*fy)**2

        gamma = 0.0
        Xi = self.Et

        # plastic corrector as needed
        if f >= 0.0:
            # loop
            r = Phi @ stress
            Xixr = Xi @ r
            gamma += f / (r @ Xixr)
            Xi = np.linalg.inv(Cinv + gamma * Phi)
            stress = Xi @ ( eps - self.plastic_strain )
            self.Et  = Xi

            self.strain['zz']  = -0.5 * (self.plastic_strain[0] + self.plastic_strain[1])
            self.strain['zz'] += -nu * (stress[0] + stress[1]) / E

            print("material entering plastic state")

        #print(4*'{:12.8e}  '.format(eps, f, self.plastic_strain, self.sig ))
        self.sig    = stress
        self.stress = {'xx':stress[0], 'yy':stress[1], 'zz':0.0, 'xy':stress[2], 'yz':0.0, 'zx':0.0}

    def getStrain(self):
        return self.strain

    def converged(self):
        # update state now that the global analysis has converged
        E  = self.parameters['E']
        t  = self.parameters['t']
        nu = self.parameters['nu']

        Cinv = 1 / (E * t) * np.array([[1., -nu, 0.], [-nu, 1., 0.], [0., 0., 2. * (1. + nu)]])
        self.plastic_strain = self.total_strain - Cinv @ self.sig
