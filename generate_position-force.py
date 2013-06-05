#!/usr/bin/env python

# Script to create a position-force file from
# a chosen analytic potential V(r) and its
# derived force function = V'(r)

# NOTE: THE main() FUNCTION MUST BE EDITED SO TO USE
#       THE DESIRED FORCE FUNCTION!!! (feel free to
#       add force functions to the code)

#=== GLOBAL VARIABLES ===#
global r0
global rN
global N
r0 = 3.200
rN = 12.4125
N = 1000000

import os, sys
import Numeric, math
import random

  #========================#
  # CREATE FORCE FUNCTIONS #
  #========================

def F_hard_sphere(r):
    """ Based on: V(r) = 1000/(r - 1.0)^12 """
    F = -12000/(r - 1.0)**13
    return F

def F_quadratic(r):
    """ Based on: V(r) = r^2 """
    F = 2*r
    return F    

def F_cos(r):
    """ Based on: V(r) = sin(r) """
    F = math.cos(r)
    return F

  #==============#
  # RUN THE CODE #
  #==============#
  
def main():
    """
    Choose an array of certain length and stepsize
    to sufficiently represent the desired interatomic
    distances.
    """
#    r = Numeric.arange(r0,rN+dr,dr,typecode=Numeric.Float)
    r = [(random.random()*(rN-r0)+r0) for i in range(N)]
    F = []
    out = open('position-force.dat','w')
    
    for i in range(len(r)):
        
        # CHOOSE THE FORCE FUNCTION YOU WANT
        F.append( F_hard_sphere(r[i]) )
#        F.append( F_quadratic(r[i]) )
#        F.append( F_cos(r[i]) )        

        out.write(str(i+1)+'\n')
        out.write('0.0  0.0  0.0  '+str(F[i])+'  0.0  0.0\n')
        out.write(str(r[i])+'  0.0  0.0  '+str(-F[i])+'  0.0  0.0\n')

    out.write(str(i+2))
    
    out.close()

if __name__ == '__main__':
    main()
