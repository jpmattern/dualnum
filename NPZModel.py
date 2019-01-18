# -*- coding: utf-8 -*-

import modelcomponents
from dualnum import DualNumber
import numpy as np

i_irr = 0
i_nut = 1
i_phy = 2
i_zoo = 3

class NPZModel(object):
    ''' A simple NPZ model formulation containing tangent linear and adjoint code.
    
Methods
-------
run_nl: 
    Run the nonlinear model.
run_tl: 
    Run the tangent linear model.
run_ad: 
    Run the adjoint model.
    '''
    def __init__(self, delta_t=0.01, lightparameters=None,
                 lightresponse=None, nutrientuptake=None, grazing=None, 
                 phytoplanktonloss=None, zooplanktonloss=None, **parameters):
        '''Initialize an NPZ model.

Parameters
----------
delta_t: float, optional
    The length of a model time step.
    
lightparameters: dict, optional 
    A dict containing the parameters for a light parametrization.
    Valid options: 
    constant light: {'type':'const', 'value':<constant light value>}
    sinusiodal: {'type':'sinusoidal', 'min':<minimum light value>,
                 'max':<maximum light value>,
                 'wavelength':<wavelength of the sine function>}
lightresponse: function, optional
    A light response function (see examples in modelcomponents.py).
nutrientuptake: function, optional
    A nutrient uptake function (see examples in modelcomponents.py).
grazing: function, optional
    A zooplankton grazing on phytoplankton function (see examples in 
    modelcomponents.py).
phytoplanktonloss: function, optional
    A phytoplankton loss function (see examples in modelcomponents.py).
zooplanktonloss: function, optional
    A zooplankton loss function (see examples in modelcomponents.py).
**parameters: float, optional
    One or multiple parameter values.
    '''
        
        # set all functions, if they are not specified by the user,
        # select defaults from modelcomponents
        if lightresponse is None:
            lightresponse = modelcomponents.lightresponse_linear
        if nutrientuptake is None:
            nutrientuptake = modelcomponents.nutrientuptake_michaelismenten
        if grazing is None:
            grazing = modelcomponents.grazing_linear
        if phytoplanktonloss is None:
            phytoplanktonloss = modelcomponents.phytoplanktonloss_linear
        if zooplanktonloss is None:
            zooplanktonloss = modelcomponents.zooplanktonloss_quadratic
        
        self.lightresponse = lightresponse
        self.nutrientuptake = nutrientuptake
        self.grazing = grazing
        self.phytoplanktonloss = phytoplanktonloss
        self.zooplanktonloss = zooplanktonloss
        
        # load standard parameters and update them with user-supplied values
        self.parameters = modelcomponents.standard_parameters.copy()
        self.parameters.update(parameters)
        
        self.delta_t = delta_t
        
        # set lightparameters
        if lightparameters is None:
            self.lightparameters = {'type':'sinusoidal',
                                    'min':0.5*self.parameters['irr0'],
                                    'max':1.5*self.parameters['irr0'],
                                    'wavelength':250.0*self.delta_t,
                                    }
        else:
            self.lightparameters = lightparameters
        # test lightparameters briefly
        self.generate_light(num_t=5)
        
        # the stepsize h for the dual number-based tangent linear and adjoint models
        self.h_tl = 1e-6
        self.h_ad = 1e-6
    
    def __str__(self):
        return '''NPZModel
     lightresponse: "{}"        
    nutrientuptake: "{}"
           grazing: "{}"
 phytoplanktonloss: "{}"
   zooplanktonloss: "{}"
'''.format(self.lightresponse.__name__,
           self.nutrientuptake.__name__,
           self.grazing.__name__,
           self.phytoplanktonloss.__name__,
           self.zooplanktonloss.__name__)

    def __repr__(self):
        return '''NPZModel(delta_t={},
         lightresponse={},        
         nutrientuptake={},
         grazing={},
         phytoplanktonloss={},
         zooplanktonloss={})
'''.format(self.delta_t,
           self.lightresponse.__name__,
           self.nutrientuptake.__name__,
           self.grazing.__name__,
           self.phytoplanktonloss.__name__,
           self.zooplanktonloss.__name__)
        
    def generate_light(self, num_t):
        if self.lightparameters['type'] == 'const':
            light = np.full(shape=num_t, fill_value=self.lightparameters['value'])
        elif self.lightparameters['type'] == 'sinusoidal':
            light = np.sin(np.arange(num_t)*(self.delta_t*2.*np.pi/self.lightparameters['wavelength']))
            light = (light+1.0)*(0.5*(self.lightparameters['max']-self.lightparameters['min']))+self.lightparameters['min']
        else:
            raise ValueError('Invalid light parametrization.')
        if np.any(light < 0.):
            raise ValueError('Light parametrization creates negative values.')
        return light
    
    #
    # the nonlinear (nl) code
    #
        
    def p_growth(self, x_nl):
        #
        # a phytoplankton response to light
        # and nutrient uptake formulation
        #
        growth  = self.lightresponse(irr=x_nl[i_irr], parameters=self.parameters)
        growth *= self.nutrientuptake(nut=x_nl[i_nut], parameters=self.parameters)
        growth *= x_nl[i_phy]
        
        # apply growth to nutrients and phytoplankton
        x_nl[i_nut] -= growth*self.delta_t
        x_nl[i_phy] += growth*self.delta_t
        
        return x_nl

    def z_grazing(self, x_nl):
        #
        # a zooplankton grazing on
        # phytoplankton formulation
        #
        grazing  = self.grazing(phy=x_nl[i_phy], parameters=self.parameters)
        grazing *= x_nl[i_zoo]
        
        # apply grazing to phytoplankton and zooplankton
        x_nl[i_phy] -= grazing*self.delta_t
        x_nl[i_zoo] += grazing*self.delta_t
        
        return x_nl
    
    def p_loss(self, x_nl):
        #
        # a phytoplankton loss formulation
        #
        
        loss  = self.phytoplanktonloss(phy=x_nl[i_phy], parameters=self.parameters)
        loss *= x_nl[i_phy]
        
        # apply loss to phytoplankton and nutrients
        x_nl[i_phy] -= loss*self.delta_t
        x_nl[i_nut] += loss*self.delta_t
        
        return x_nl
    
    def z_loss(self, x_nl):
        #
        # a zooplankton loss formulation
        #
        loss  = self.zooplanktonloss(zoo=x_nl[i_zoo], parameters=self.parameters)
        loss *= x_nl[i_zoo]
        
        # apply loss to zooplankton and nutrients
        x_nl[i_zoo] -= loss*self.delta_t
        x_nl[i_nut] += loss*self.delta_t
        
        return x_nl
    
    def run_nl(self, npz_ini, num_t):
        '''Run the nonlinear model.

Parameters
----------
npz_ini: array-like
    The initial conditions for nutrients, phytoplankton, and zooplankton
num_t: int
    The number of time steps to run the model for.

Returns
-------
x_nl_history: array
    An `num_t`-by-4 array containing the nonlinear model state at each time 
    step.
        '''
        # generate light
        light = self.generate_light(num_t+1)
        
        # create history
        x_nl_history = np.full(shape=(num_t+1,4), fill_value=np.nan)
        x_nl_history[0,i_irr] = light[0]
        x_nl_history[0,1:] = npz_ini
        
        # current conditions
        x_nl = x_nl_history[0,:].copy()
        
        for t in range(num_t):
            
            # call each model segment in sequence
            self.p_growth(x_nl)
            self.z_grazing(x_nl)
            self.p_loss(x_nl)
            self.z_loss(x_nl)
   
            # obtain light
            x_nl[i_irr] = light[t+1]
            
            # record history
            x_nl_history[t+1,:] = x_nl
            
        return x_nl_history

    #
    # the dual number-based tangent linear (tl) code
    #
    
    # the dual number-based tangent linear phytoplankton growth function is written out 
    # explicitly, all others use the _generic_tl function below
    def p_growth_tl(self, x_nl, x_tl):
        
        # create a dual number array
        # using the entries of x_nl as the real parts and
        # entries of x_tl as dual parts
        x_dtl = np.array([DualNumber(x,x_d) for x,x_d in zip(x_nl,x_tl)])
        
        # call p_growth with dual number array
        x_dtl = self.p_growth(x_dtl)
        
        # copy dual parts back into x_tl
        x_tl = np.array([dn.x_d for dn in x_dtl])
        
        return x_tl
    
    # A function providing the dual number-based tangent linear code for 
    # any nonlinear model segment function (generic_nl).
    # For usage see definitions of z_grazing_tl, p_loss_tl, etc. below.
    def _generic_tl(self, generic_nl, x_nl, x_tl):
        
        # create a dual number array
        # using the entries of x_nl as the real parts and
        # entries of x_tl as dual parts
        x_dtl = np.array([DualNumber(x,x_d) for x,x_d in zip(x_nl,x_tl)])
        
        # call generic_nl with dual number array
        x_dtl = generic_nl(x_dtl)
        
        # copy dual parts back into x_tl
        x_tl = np.array([dn.x_d for dn in x_dtl])
        return x_tl
    
    def z_grazing_tl(self, x_nl, x_tl):
        return self._generic_tl(self.z_grazing, x_nl, x_tl)
        
    def p_loss_tl(self, x_nl, x_tl):
        return self._generic_tl(self.p_loss, x_nl, x_tl)
        
    def z_loss_tl(self, x_nl, x_tl):
        return self._generic_tl(self.z_loss, x_nl, x_tl)
        
    def run_tl(self, npz_ini, x_tl_ini, num_t):        
        '''Run the dual number-based tangent linear model.

Parameters
----------
npz_ini: array-like
    The initial conditions for nutrients, phytoplankton, and zooplankton.
x_tl_ini: array-like
    The initial conditions for the tangent linear model (including light, 
    nutrients, phytoplankton, and zooplankton).
num_t: int
    The number of time steps to run the model for.
    
Returns
-------
x_nl_history: array
    An `num_t`-by-4 array containing the nonlinear model state at each time 
    step.
x_tl_history: array
    An `num_t`-by-4 array containing the tangent linear model state at each 
    time step.
        '''
        
        # generate light
        light = self.generate_light(num_t+1)
        
        # create histories
        x_nl_history = np.full(shape=(num_t+1,4), fill_value=np.nan)
        x_nl_history[0,1:] = npz_ini
        x_nl_history[0,i_irr] = light[0]
        
        x_tl_history = np.full(shape=(num_t+1,4), fill_value=np.nan)
        x_tl_history[0,:] = x_tl_ini
        
        # current conditions
        x_nl = x_nl_history[0,:].copy()
        x_tl = x_tl_history[0,:].copy()
        
        for t in range(num_t):
            # for each of the 4 segments: 
            # step (1) update the tangent linear 
            # state x_tl using the correct 
            # value of the nonlinear state x_nl
            # step (2) propagate the nonlinear 
            # state forward
            
            # segment: P growth
            x_tl = self.p_growth_tl(x_nl, x_tl)
            x_nl = self.p_growth(x_nl)
            
            # segment: Z grazing
            x_tl = self.z_grazing_tl(x_nl, x_tl)
            x_nl = self.z_grazing(x_nl)
            
            # segment: P loss
            x_tl = self.p_loss_tl(x_nl, x_tl)
            x_nl = self.p_loss(x_nl)
            
            # segment: Z loss
            x_tl = self.z_loss_tl(x_nl, x_tl)
            x_nl = self.z_loss(x_nl)
            
            # obtain light
            x_nl[i_irr] = light[t+1]
            
            # record history
            x_nl_history[t+1,:] = x_nl
            x_tl_history[t+1,:] = x_tl
            
        return x_nl_history, x_tl_history
    
    #
    # the dual number-based adjoint (ad) code
    #
    
    # A function providing the dual number-based adjoint code for any 
    # nonlinear model segment function (generic_nl).
    # For usage see definitions of p_growth_ad, z_grazing_ad etc. below.
    def _generic_ad(self, generic_nl, ivars_in, x_nl, x_ad):
        
        # create copy of adjoint state 
        x_ad_orig = x_ad.copy()
        
        # n is the number of input variables in ivars_in
        n = len(ivars_in)

        # Create a dual number with n dual parts.
        # For now, all dual parts are zero
        # (np.zeros(n) creates a NumPy array with n zeros).
        # For NumPy arrays, operations like +,-,*,/, etc. are
        # performed element-wise, so no extension to the
        # DualNumber class is necessary to permit more than
        # one dual part.
        x_dad = np.array([DualNumber(x,np.zeros(n)) for x in x_nl])

        # Set the i-th dual part to one for the dual number 
        # associated with the i-th input variable.
        # example:
        # For ivars_in = (1,3) and hence n = 2 
        # x_dad is an array of length 4 with the following entries:
        #   x_dad[0] = DualNumber(x_nl[0], np.array([0.0, 0.0]))
        #   x_dad[1] = DualNumber(x_nl[1], np.array([1.0, 0.0]))
        #   x_dad[2] = DualNumber(x_nl[2], np.array([0.0, 0.0]))
        #   x_dad[3] = DualNumber(x_nl[3], np.array([0.0, 1.0]))
        for i, iv_in in enumerate(ivars_in):
            x_dad[iv_in].x_d[i] = 1.0

        # call nonlinear function generic_nl with dual number array
        x_dad = generic_nl(x_dad)
        
        # contruct x_ad from the x_dad:
        # The new entry x_ad[iv_in], associated with the i-th input 
        # variable, is the scalar product of 
        # x_dad[:].x_d[i] (the i-th dual part of all 4 dual numbers) 
        # and the original values in x_ad
        # (variables that are not input variables, are not modified
        # in the adjoint model).
        for i, iv_in in enumerate(ivars_in):
            delta_ad = 0.0
            for j in range(4):
                delta_ad += x_dad[j].x_d[i] * x_ad_orig[j]
            x_ad[iv_in] = delta_ad

        return x_ad
    
    # based on self._generic_ad, it is  easy to create the adjoint code 
    # for each model segment

    # for p_growth: light, N, and P act as input
    def p_growth_ad(self, x_nl, x_ad):
        return self._generic_ad(self.p_growth, (i_irr,i_nut,i_phy), x_nl, x_ad)
    
    # for z_grazing: P and Z act as input
    def z_grazing_ad(self, x_nl, x_ad):
        return self._generic_ad(self.z_grazing, (i_phy,i_zoo), x_nl, x_ad)
    
    # for p_loss: P acts as input
    def p_loss_ad(self, x_nl, x_ad):
        return self._generic_ad(self.p_loss, (i_phy,), x_nl, x_ad)
    
    # for z_loss: Z acts as input
    def z_loss_ad(self, x_nl, x_ad):
        return self._generic_ad(self.z_loss, (i_zoo,), x_nl, x_ad)
    
    def run_ad(self, npz_ini, x_ad_ini, num_t=None, x_nl_history=None):
        '''Run the dual number-based adjoint model.

Parameters
----------
npz_ini: array-like
    The initial conditions for nutrients, phytoplankton, and zooplankton.
x_ad_ini: array-like
    The initial conditions for the adjoint model (including light, 
    nutrients, phytoplankton, and zooplankton).
num_t: int, optional
    The number of time steps to run the model for (optional if `x_nl_history`
    is supplied).
x_nl_history: array
    An `num_t`-by-4 array containing the nonlinear model state at each time 
    step. If `x_nl_history` is provided, it is used as the nonlinear model 
    state in the computation of the adjoint. Otherwise `x_nl_history` is 
    computed using `npz_ini` and `num_t`.
    
Returns
-------
x_nl_history: array
    An `num_t`-by-4 array containing the nonlinear model state at each time 
    step.
x_ad_history: array
    An `num_t`-by-4 array containing the adjoint model state at each time 
    step.
        '''
        
        # run nonlinear 
        if x_nl_history is None:
            x_nl_history = self.run_nl(npz_ini=npz_ini, num_t=num_t)
        
        #print(x_nl_history)
        
        x_ad_history = np.full(shape=(num_t+1,4), fill_value=np.nan)
        x_ad_history[-1,:] = x_ad_ini
        
        x_ad = x_ad_history[-1,:].copy()
       
        for t in reversed(range(num_t)):
            # load appropriate nonlinear state
            x_nl = x_nl_history[t,:].copy()
            
            # First, propagate the nonlinear model 
            # solution and save it at every stage.
            
            # segment: P growth
            x_nl_p_growth = x_nl.copy()
            x_nl = self.p_growth(x_nl)
            
            # segment: Z grazing
            x_nl_z_grazing = x_nl.copy()
            x_nl = self.z_grazing(x_nl)
            
            # segment: P loss
            x_nl_p_loss = x_nl.copy()
            x_nl = self.p_loss(x_nl)
            
            # segment: Z loss
            x_nl_z_loss = x_nl.copy()
            x_nl = self.z_loss(x_nl)
            
            # Second, propagate adjoint model 
            # solution in reverse order using 
            # the appropriate saved nonlinear 
            # model solution.
            # Use saved nonlinear model solution 
            # for the next segment as the 
            # reference solution for the current 
            # segment (3rd input argument).
            
            # segment: Z loss
            x_ad = self.z_loss_ad(x_nl_z_loss, x_ad)
            
            # segment: P loss
            x_ad = self.p_loss_ad(x_nl_p_loss, x_ad)
            
            # segment: Z grazing
            x_ad = self.z_grazing_ad(x_nl_z_grazing, x_ad)
            
            # segment: P growth
            x_ad = self.p_growth_ad(x_nl_p_growth, x_ad)
            
            # record history
            x_ad_history[t,:] = x_ad
        return x_nl_history, x_ad_history
    
if __name__=='__main__':
    m = NPZModel()
    
    x_nl = m.run_nl(npz_ini=(10.,5.,5.), num_t=500)
    
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots()
    ax.plot(x_nl)
    ax.legend(('I','N','P','Z'))
    ax.grid(True)
    ax.set(title='nonlinear model simulation', xlabel='time step', 
           ylabel='variable value')
    
    # plot tl
    
    x_nl, x_tl = m.run_tl(npz_ini=(10.,5.,5.), x_tl_ini=(1.0,1.0,0.2,0.2), num_t=500)
        
    colors = ('C0','C1','C2','C3')
    
    fig, ax = plt.subplots()
    
    for i,name in enumerate('INPZ'):
        ax.plot(x_nl[:,i], color=colors[i], label=name)
        ax.plot(x_nl[:,i]+x_tl[:,i], color=colors[i], alpha=0.5)
        
    ax.legend()
    ax.grid(True)
    ax.set(title='nonlinear and tangent linear model simulation', 
           xlabel='time step', ylabel='variable value')
    
    plt.show()
    
