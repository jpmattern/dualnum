# -*- coding: utf-8 -*-

import numpy as np

#
# The parameters used in the functions below.
#

standard_parameters = {
    # baseline irradiance parameter
    'irr0':5.0,
    # maximum rate in Michaelis Menten formulation
    'Vmax':10.0,
    # nutrient half saturation in Michaelis Menten formulation
    'nuthalfsat':0.5,
    # multiplicative grazing parameter
    'grazphy':0.25,
    # grazing parameter used in exponential functions
    'grazlambda':0.5,
    # maximum grazing rate
    'grazmax':0.25,
    # phytoplankton mortality rate
    'mort_phy':0.2,
    # zooplankton mortality rate
    'mort_zoo':0.1,
    }

#
# A selection of light response functions. Compare Table 1 in Franks (2002).
#

def lightresponse_linear(irr, parameters):
    return irr/parameters['irr0']

def lightresponse_saturating(irr, parameters):
    return irr/(parameters['irr0']+irr)

def lightresponse_exp(irr, parameters):
    return 1.0 - np.exp(-irr/parameters['irr0'])

def lightresponse_tanh(irr, parameters):
    return np.tanh(-irr/parameters['irr0'])

def lightresponse_inhibit(irr, parameters):
    irr_norm = irr/parameters['irr0']
    return irr_norm * np.exp(1.0-irr_norm)

#
# A selection of nutrient uptake functions. Compare Table 2 in Franks (2002).
#

def nutrientuptake_michaelismenten(nut, parameters):
    return parameters['Vmax']/(parameters['nuthalfsat']+nut)

#
# A selection of zooplankton grazing functions. Compare Table 3 in Franks (2002).
#

def grazing_linear(phy, parameters):
    return parameters['grazphy']*phy

def grazing_bilinear(phy, parameters):
    return np.min(parameters['grazphy']*phy,parameters['grazmax'])

def grazing_ivlev(phy, parameters):
    return parameters['grazmax']*(1.0 - np.exp(-parameters['grazlambda']*phy))

#
# A selection of phytoplankton loss functions. Compare Table 4 in Franks (2002).
#

def phytoplanktonloss_linear(phy, parameters):
    return parameters['mort_phy']

def phytoplanktonloss_quadratic(phy, parameters):
    return parameters['mort_phy']*phy

#
# A selection of zooplankton loss functions. Compare Table 4 in Franks (2002).
#

def zooplanktonloss_linear(zoo, parameters):
    return parameters['mort_zoo']

def zooplanktonloss_quadratic(zoo, parameters):
    return parameters['mort_zoo']*zoo

#
# A generic function that can be used in place of any of the above in order to
# "switch off" a given segment. Using generic_nomod as the zooplankton grazing
# function, for example, will turn zooplankton grazing to zero.
#

def generic_nomod(*args, **kwargs):
    return 0.0

