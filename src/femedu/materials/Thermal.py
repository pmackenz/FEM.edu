from .DiffusionGeneral import *

class Thermal(DiffusionGeneral):
    """
    Thermal diffusion material description

    .. list-table:: parameters inside the **params** dictionary
       :header-rows: 1

       * - key
         - default value
         - description
       * - **conductivity**
         - 1.0
         - thermal conductivity
       * - **specifict heat**
         - 1.0
         - specific heat capacity
       * - **density**
         - 1.0
         - mass density

    :param params: isotropic thermal properties
    :type params: dict
    """
    
    def __init__(self, params={'conductivity':1., 'specific_heat':1., 'density':1., 'thickness':1.}):

        valid_params = ('conductivity','specific_heat','density','thickness')

        for p in valid_params:
            if p not in params:
                params[p] = 1.0

        for p in params:
            if p not in valid_params:
                msg = f"Warning: unknown parameter {p} found in parameter list - ignored"
                print(msg)

        general_params = dict(
            diffusivity = params['conductivity'],
            capacity    = params['specific_heat'],
            density     = params['density'],
            thickness   = params['thickness']
        )
        super(Thermal, self).__init__(params=general_params)
