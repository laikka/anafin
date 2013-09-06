#! /usr/bin/env python
"""Analytical solver for fin temperature. 
    
    Modeled after the PhD thesis of Antti Lehtienen "Analytical Treatment of
    Heat Sinks Cooled by Forced Convection", Tampere University of Technology,
    2005.
    
    Written by Antti Mikkonen, a.mikkonen@iki.fi, fall 2013
    
    This program was written for the fun of it by a curious PhD student. 
    The program is supposed to be functional but no tedious testing was 
    performed. The code can be used for any purpose as long as the source  
    is visibly cited. 
                
"""
from __future__ import division

import numpy as np
from math import pi, sqrt, exp

class Fin(object):
    """Class for fin temperature solving.
    
    The equation and table references refer to the thesis by Lentinen. Same
    notations as in the thesis are used when convenient.
                    
    """


    def __init__(self, in_values):
        self.__dict__.update(in_values)
        
        self.b = self.b2 / 2.
        self.t = self.t2 / 2.
    
        self.Re = 4 * self.U * self.b / self.nu
        self.Pr = self.rho * self.c_p * self.nu / self.k_a
        
        # Dimensionless fin lenght in flow direction eqref{2.9}
        self.L_plus = self.L / (2 * self.b * self.Re * self.Pr)
    
    def Graetz_eigenvalue(self,n):
        """ tableref{2.1}
        
        """
        if   n == 0:    lambda_n = 3.885
        elif n == 1:    lambda_n = 13.09
        elif n == 2:    lambda_n = 22.32
        else:           lambda_n = 16 * sqrt(n/3) + 20/3*sqrt(1/3)
        return lambda_n    
            
    def Graetz_coefficient(self, n):
        """ tableref{2.1}
        
        """
        if   n == 0:    G_n = 1.717
        elif n == 1:    G_n = 1.139
        elif n == 2:    G_n = 0.952
        else:           G_n = 2.68 * self.Graetz_eigenvalue(n)**(-1/3)
        return G_n
    
    def delta(self, i):
        """eqref{A:.}
        """
        if i == 0:  return 1
        else:       return 0
    
    def solve_fin_xy_temperature_for_laminar(self,x,y):
        self.calculateM(x,y)
    
    def calculateM(self,x,y):
        """ eqref{5.14}
        
        """
        # init
        M2 = np.zeros([self.N,self.N])
    
        # diagonal elemets, eqref{5.14a}
        for i in range(self.N):
            
            alpha_i = i * pi / self.L 

            sum_term = 0            
            for n in range(self.N):
                
                G_n = self.Graetz_coefficient(n)
                lambda_n = self.Graetz_eigenvalue(n)
                delta_i = self.delta(i)
                
                sum_term += G_n \
                            * ( (i * pi)**2 / \
                                    ( (lambda_n**2 * self.L_plus)**2 \
                                      + (i * pi)**2) \
                                + 2 * (lambda_n**2 * self.L_plus)**3 \
                                  / (1 + delta_i) \
                                  * (1 - (-1)**i \
                                     * exp(-lambda_n**2 * self.L_plus)) \
                                  / ((lambda_n**2 * self.L_plus)**2 \
                                     + (i * pi)**2)**2
                                 )
             
             
             
            M2[i,i] = alpha_i**2 \
                        + self.k_a / ( self.k_f * self.b * self.t)  * sum_term
      
        print M2
    
    
    
    
    
    
    
    
    
    
    
    
    
##################
# Tester
##################    

def test():
    in_values = {
                  # Geometry
                  "l" : 0.05,
                  "L" : 0.10,
                  "t2" : 0.001,
                  "b2" : 0.001,
                  
                  # Fluid props
                  "k_a" : 0.0289,
                  "rho" : 1,
                  "c_p" : 1006,
                  "nu"  : 1.9e-05,
                  "U"   : 5,
                                   
                  # Fin props
                  "k_f" : 210,
        
                  # Calc props
                  "N" : 5,
                            
                  }
    
    fin = Fin(in_values)
    fin.solve_fin_xy_temperature_for_laminar(0.05, 0.025)

    
    
# *****************************************************************************
if __name__ == "__main__":
    test()
    print "DONE"
