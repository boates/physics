#!/usr/bin/env python

import Numeric, RandomArray, matplotlib, pylab

def lregress(x,y):
    """Fits a straight line to data points (x,y).

    The line to be fit is y = A + Bx

    Returns: A, B, dy, dA, dB
    """
    x = Numeric.array(x,typecode=Numeric.Float)
    y = Numeric.array(y,typecode=Numeric.Float)

    N = len(x)

    DELTA = N * Numeric.sum(x**2) - Numeric.sum(x)**2
    
#    print Numeric.shape(Numeric.sum(x**2)),Numeric.shape(Numeric.sum(y))
    
    A = ( Numeric.sum(x**2)*Numeric.sum(y)-Numeric.sum(x)*Numeric.sum(x*y) ) / DELTA

    B = ( N*Numeric.sum(x*y) - Numeric.sum(x)*Numeric.sum(y) )/ DELTA

    sigma_y = Numeric.sqrt( (1./(N-2)) * Numeric.sum((y - A - B*x)**2) )

    sigma_A = sigma_y * Numeric.sqrt( Numeric.sum(x**2)/DELTA )

    sigma_B = sigma_y * Numeric.sqrt( N/DELTA )

    return A,B,sigma_y,sigma_A,sigma_B


# Test program
if __name__== '__main__':

    # Use latex text formatting
    matplotlib.rc('text', usetex=True)

    # Create a test data set
    x = Numeric.arange(20)
    dy = RandomArray.standard_normal(20)
    y = 0.5*x+3 + dy
 
    A,B,dy,dA,dB = lregress(x,y)
    y2 = A + B*x
    print A, B, dy ,dA,dB

#    pylab.plot(x,y,'o')
#    pylab.plot(x,y2,linewidth=3)
#    pylab.xlabel('x',fontsize='large')
#    pylab.ylabel('y',fontsize='large')

#    pylab.text(2, 13, 'A = %.1f $\pm$ %.1f (true: 3.0)' % (A,dA) )
#    pylab.text(2, 12, 'B = %.2f $\pm$ %.2f (true: 0.5)' % (B,dB) )
#    pylab.text(2, 11, 'dy = %.1f (true: 1.0)' % (dy) )

#    pylab.axis( (0,20,2,15) )

#    pylab.subplots_adjust(left=0.3,right=0.7,top=0.7,bottom=0.3)
#    pylab.show()
