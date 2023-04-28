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
    
    def __init__(self, params):

        if 'conductivity' not in params:
            params['conductivity'] = 1.0
        if 'specific heat' not in params:
            params['specific heat'] = 1.0
        if 'density' not in params:
            params['density'] = 1.0

        general_params = dict(
            diffusivity = params['conductivity'],
            capacity    = params['specific heat'],
            density     = params['density']
        )
        super(Thermal, self).__init__(params=general_params)
