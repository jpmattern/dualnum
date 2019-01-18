#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Perform dot-product test, to check consistency between tangent linear and 
adjoint models.

The dot-product is a test for consistency based on a simple relationship, true 
for any suitable vectors x_tl and x_ad and suitable matrix, including the 
tangent linear matrix M:
    
    x_ad^T * M * x_tl = (M^T * x_ad)^T * x_tl
    
It follows that

    x_ad^T * M_tl(x_tl) = M_ad(x_ad)^T * x_tl
    
where M_tl and M_ad are the tangent linear and adjoint models, respectively.
The difference between the left and right hand side of the above equation acts 
as an indicator of the consistency between the two models.

In the code below, the dot-product test is applied to the FDA-based tangent 
linear and adjoint models.
'''

import numpy as np
from NPZModel import NPZModel

# create model and initial conditions
m = NPZModel()
npz_ini = (8.,5.,3.)

# perform one test each for different num_t (number of time steps)
for num_t in (1,3,10,30,100,300,1000): 

    # randomly generate x_tl and x_ad
    x_tl = np.random.uniform(size=4, low=-1.0, high=1.0)
    x_ad = np.random.uniform(size=4, low=-1.0, high=1.0)
    
    # compute the results of tangent linear and adjoint models
    res_tl = m.run_tl(npz_ini=npz_ini, x_tl_ini=x_tl, num_t=num_t)[1][-1,:]
    res_ad = m.run_ad(npz_ini=npz_ini, x_ad_ini=x_ad, num_t=num_t)[1][0,:]
    
    '''
    # the following could be used to test an individual model segment
    
    x_nl = np.array((10.,8.,5.,3.))
    res_tl = m.p_growth_tl(x_nl=x_nl.copy(), x_tl=x_tl.copy())
    res_ad = m.p_growth_ad(x_nl=x_nl.copy(), x_ad=x_ad.copy())
    '''
    
    # compute the left and right hand side of the equations above
    res0 = x_ad.dot(res_tl)
    res1 = x_tl.dot(res_ad)

    # ... and compute the absolute difference
    print('num_t = {:4d}, abs. difference: {:.3e}'.format(num_t,abs(res0-res1)))
