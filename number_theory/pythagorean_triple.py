#!/usr/bin/env python
"""
Finds all Pythagorean triples ( (x,y,z) satisfying x^2 + y^2 = z^2
where x, y, & z have no factor in common ), given a certain range of values
"""
import math
from Numeric import arange, Float

global MIN
global MAX
MIN = 10
MAX = 30

def pythagorean_triples(top, bottom=2): 

    triples = []
    x = arange(bottom, top+1, typecode=Float)
    y = arange(bottom, top+1, typecode=Float)

    for i in range(len(x)):

        for j in range(len(y)):
            z = ( x[i]**2 + y[j]**2 )**.5
            sharing = factor_check(x[i],y[j],z)

            if z == math.floor(z) and sharing == False:
                triple = (x[i],y[j],z)

                if (y[j],x[i],z) not in triples:
                    triples.append(triple)

    return triples

def factor_check(x,y,z):

    x_factors = []
    y_factors = []
    z_factors = []
    
    check_x = arange(3, x+1,typecode=Float)
    check_y = arange(3, y+1,typecode=Float)
    check_z = arange(3, z+1,typecode=Float)

    for i in range(len(check_x)):
        if x % check_x[i] == 0:
            x_factors.append(check_x[i])

    for j in range(len(check_y)):
        if y % check_y[j] == 0:
            y_factors.append(check_y[j])

    for k in range(len(check_z)):
        if z % check_z[k] == 0:
            z_factors.append(check_z[k])

    for factor in x_factors:
        if factor in y_factors:
            x_share_y = True
            break
        else: x_share_y = False

    for factor in x_factors:
        if factor in z_factors:
            x_share_z = True
            break
        else: x_share_z = False

    for factor in y_factors:
        if factor in z_factors:
            y_share_z = True
            break
        else: y_share_z = False

    if y_share_z or x_share_z or x_share_y:
        sharing = True
    else: sharing = False

    return sharing
        
def main():

    """ Print out the Pythagorean triples for the defined range """

    triples = pythagorean_triples(MAX,bottom=MIN)

    print '\nPythagorean triples for '+str(MIN)+' <= x, y <= '+str(MAX)+':'
    for triple in triples:
        print triple
    print ''
    
if __name__ == '__main__':
    main()

