#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class DualNumber(object):
    def __init__(self,x=0.0,x_d=0.0):
        self.x = x
        self.x_d = x_d
        
    def __repr__(self):
        return 'DualNumber({},{})'.format(self.x,self.x_d)
    
    # overload operators:
    
    # self + other
    def __add__(self,other):
        if isinstance(other,DualNumber):
            return DualNumber(self.x + other.x, self.x_d + other.x_d)
        return DualNumber(self.x + other, self.x_d)
      
    # other + self
    def __radd__(self,other):
        return self.__add__(other)
    
    # self - other
    def __sub__(self,other):
        if isinstance(other,DualNumber):
            return DualNumber(self.x - other.x, self.x_d - other.x_d)
        return DualNumber(self.x - other, self.x_d)
      
    # other - self
    def __rsub__(self,other):
        if isinstance(other,DualNumber):
            return DualNumber(other.x - self.x, other.x_d - self.x_d)
        return DualNumber(other - self.x, -self.x_d)
    
    # self * other
    def __mul__(self,other):
        if isinstance(other,DualNumber):
            return DualNumber(self.x * other.x, self.x*other.x_d + self.x_d*other.x)
        return DualNumber(self.x * other, self.x_d * other)
    
    # other * self
    def __rmul__(self,other):
        return self.__mul__(other)
    
    # self / other
    def __truediv__(self,other):
        if isinstance(other,DualNumber):
            xsq = other.x*other.x
            return DualNumber( self.x / other.x, (self.x_d*other.x - self.x*other.x_d)/xsq)
        return DualNumber(self.x / other, self.x_d / other)
    
    # other / self
    def __rtruediv__(self,other):
        if isinstance(other,DualNumber):
             xsq = self.x*self.x
             return DualNumber( other.x / self.x, (other.x_d*self.x - other.x*self.x_d)/xsq)
        xsq = self.x*self.x
        return DualNumber(other/self.x , -other*self.x_d/xsq)
    
    # self ** other
    def __pow__(self,other):
        if isinstance(other,int) or isinstance(other,float):
            return DualNumber(self.x ** other, other * self.x**(other-1) * self.x_d )
        raise TypeError('unsupported operand type(s) for ** or pow(): \'DualNumber\' and \'{}\''.format(type(other)))
    
    # self > other
    def __gt__(self,other):
        if isinstance(other,DualNumber):
            return self.x > other.x
        return self.x > other
    
    # self < other
    def __lt__(self,other):
        if isinstance(other,DualNumber):
            return self.x < other.x
        return self.x < other
   
    # overload some mathematical functions

    def exp(self):
        tmp = np.exp(self.x)
        return DualNumber(tmp, tmp*self.x_d)
    
if __name__=='__main__':
    d = DualNumber(3.,4.)
    
    def f(x):
        return 2*(x-1)**2 + 3
       
    x = DualNumber(4,1)
    print('f({}) = {}'.format(x,f(x)))
    
