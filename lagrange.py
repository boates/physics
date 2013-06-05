#!/usr/bin/env python
import Numeric
import pylab

from Numeric import pi

def poly(x,y,xnew):
    """Lagrange polynomial interpolation.

    x: The x values.
    y: The corresponding y values.
    xnew: The x values to interpolate at.

    Returns: The interpolated y values.
    """
   

    # Force everything to be a Numeric array
    xnew = Numeric.array(xnew,typecode=Numeric.Float)
    x = Numeric.array(x,typecode=Numeric.Float)
    y = Numeric.array(y,typecode=Numeric.Float)
    
    ynew=[] # New y values corresponding to new x values

    # Do the calculation for each xnew value
    for x1 in xnew:
        l = []
        
        for j,xj in enumerate(x):

            # Numerator
            tmp = x1-x
            tmp[j] = 1.
            tmp1 = Numeric.product(tmp)

            # Denominator
            tmp = xj-x
            tmp[j] = 1.
            tmp2 = Numeric.product(tmp)
            
            l.append(tmp1/tmp2)

        ynew.append(Numeric.sum(y*Numeric.array(l)))

    return Numeric.array(ynew)


# Test program
if __name__== '__main__':

    # Create some data to fit
    x = Numeric.arange(0.,2*pi+pi/4,pi/4)
    y = Numeric.sin(x)
    y[-1]=0.8
    y[-3]=-0.5
    y[0]=-0.5

    # Fit at some higher resolution
    x1 = Numeric.arange(0.,2*pi+pi/40,pi/40)
    y1 = poly(x,y,x1)
    
    # Plot the result
    pylab.plot(x1,y1,linewidth=3)
    pylab.plot(x,y,'ro')
    pylab.axis( (0,7,-1.5,2.0) )
    pylab.xlabel('x',fontsize='large')
    pylab.ylabel('f(x)',fontsize='large')
    pylab.subplots_adjust(left=0.3,right=0.7,top=0.7,bottom=0.3)
    pylab.show()
