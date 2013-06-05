#!/usr/bin/env python
"""
Read Poten_XX.dat file from force-matching code and
fit using n'th order polynomial, extrapolating to
a finite value at r=0 and zero at r--->large
"""
import os, sys, commands, glob, optparse, numpy
from scipy import optimize
from scipy import linalg

def trim(fname):
    """
    Read in the Poten_XX.dat file and trim the beginning

    fname: name of Poten_XX.dat file
    return: trimmed V and r lists
    """
    pp = open(fname,'r')
    lines = pp.readlines()
    pp.close()
    # Remove the "flat" region of pp for small r
    V, r  = [], []
    v0 = lines[0].split()[1]  # get plateau value
    for line in lines:
        ri, vi = line.split()
        if vi != v0:
            V.append(float(vi))
            r.append(float(ri))

    # Remove first few values, where potential begins plateauing
    V = V[4:]
    r = r[4:]

    return V, r

def tail(V, r, rcut):
    """
    Create exponential tail for long range V(r)
    Form: V_tail(r) = b*e^(-a*x)

    V: potential list
    r: distance list
    rcut: distance to begin V_tail(r)
    return: V_tail, r_tail
    """
    # Need the derivative of V(r) at rcut
    dVdr, delta = [], []
    for i in range(len(r)):
        if i == 0:
            dV = V[i+1] - V[i]
            dr = r[i+1] - r[i]
        elif i == len(r)-1:
            dV = V[i] - V[i-1]
            dr = r[i] - r[i-1]
        else:
            dV = V[i+1] - V[i-1]
            dr = r[i+1] - r[i-1]
        dVdr.append( dV/dr )
        delta.append( (rcut - r[i])**2 )

    # Find new rcut
    new_rcut = r[delta.index(min(delta))]

    # c2 is the value of dV/dr at rcut
    c2 = dVdr[delta.index(min(delta))]

    # c1 is the value of V(rcut)
    c1 = V[delta.index(min(delta))]

    # Fitting parameters a and b have analytic forms
    a = - c2 / c1
    b = c1 * numpy.e**( a*new_rcut )

    dr = r[1]-r[0]
    r_tail = numpy.arange(new_rcut,max(r),dr)
    V_tail = b*numpy.e**(-a*r_tail)

    print "r > "+str(new_rcut)+":"
    print "V(r) = "+str(b)+"*e^(-"+str(a)+"*r)"

    return V_tail, r_tail

def head(V, r, t, a=None):
    """
    Create quadratic/quartic head for short range V(r)
    Form: V_head(r) = a - b*r**2
      or: V_head(r) = a - b*r**2 - c*r**4

    V: potential list
    r: distance list
    t: type of fit (quartic, quadratic, or linear)
    return: V_head, r_head
    """
    # c2 is the derivative of V(r) at r_min
    dr = r[1] - r[0]
    c2 = (V[1] - V[0]) / dr
    # c1 is the value of V(r_min)
    c1 = V[0]

    # if a is defined, use a quartic
    if t == 'quartic':
        # Fitting parameters a and b have analytic forms
        b = 2.0*(a-c1)/r[0]**2 + c2/(2.0*r[0])
        c = -c2/(2.0*r[0]**3) - (a-c1)/r[0]**4

        r_head = numpy.arange(0.0,r[0],dr)
        V_head = a - b*r_head**2 - c*r_head**4

        print "r < "+str(r[0])+":"
        print "V(r) = "+str(a)+" - "+str(b)+"*r^2 - "+str(c)+"*r^4"

    # If a isn't defined, use a quadratic
    elif t == 'quadratic':
        # Fitting parameters a and b have analytic forms
        a = c1 - (c2/2.0)*r[0]
        b = -c2 / (2.0*r[0])

        r_head = numpy.arange(0.0,r[0],dr)
        V_head = a - b*r_head**2

        print "r < "+str(r[0])+":"
        print "V(r) = "+str(a)+" - "+str(b)+"*r^2"

    elif t == 'linear':
        # Fitting parameters a and b have analytic forms
        a = c1 - c2*r[0]
        b = -c2

        r_head = numpy.arange(0.0,r[0],dr)
        V_head = a - b*r_head

        print "r < "+str(r[0])+":"
        print "V(r) = "+str(a)+" - "+str(b)+"*r"
    
    return V_head, r_head
    
def fitter(x,y,n):
    """
    Polynomial fitter from polyfit.py
    """
    # Convert to numpy.array
    x = numpy.array(x)
    y = numpy.array(y)
    
    # Build nth order basis
    basis=[]
    for i in range(n):
        basis.append(x**i)

    # Basis array transposed
    M = numpy.transpose(numpy.array(basis))
    a, resid, rank, sigma = linalg.lstsq(M,y)

    # Build nth order fit
    fit = numpy.zeros(len(x))

    for i in range(n):
        fit += a[i]*x**i

    return fit, a

def init_parser():

    usage = 'usage: %prog [options]'
    version="%prog 1.0"
    parser = optparse.OptionParser(usage=usage,version=version)

    parser.add_option('-f', '--fname', type='str', dest='fname', default=None,
                      help='Filename for pair potential file [default=None]')
    parser.add_option('-n', '--order', type='int', dest='order', default=10,
                      help='Order for polynomial fit of V(r) [default=10]')
    parser.add_option('-r', '--rcut', type='float', dest='rcut', default=4.50,
                      help='Cut-off distance for V_tail(r) to begin [default=4.50]')
    parser.add_option('-t', '--typefit', type='string', dest='t', default='linear',
                      help='Type of fit for V_head (linear, quadratic, or quartic) \
                            -a must be defined if using quartic. [default="linear"]')
    parser.add_option('-a', '--afit', type='float', dest='afit', default=None,
                      help='y-intercept for V_head, quartic fit if given, \
                            quadratic fit otherwise [default=None]')
    return parser
                    
def main():

    # Initiate option parser
    parser = init_parser()
    options, args = parser.parse_args()
    
    # Read in variables to avoid always calling dictionary
    fname = options.fname
    order = options.order+1
    rcut  = options.rcut
    t     = options.t
    a_fit = options.afit

    # Check that a filename was given
    if fname == None:
        print ' Must specify a Poten_XX.dat file with the -f flag --- exiting...'
        sys.exit(1)

    # Make sure -t was given correctly, if not - set to linear
    if t not in ['linear','quadratic','quartic']:
        print '\n -t specified incorrectly, use -h for more info - proceeding with default: t="linear"\n'
        t = 'linear'
    elif a_fit != None and t != 'quartic':
        print '\n  >  "a" was specified, but "t" was not set to quartic'
        print '     setting t="quartic" and proceeding with given "a"\n'
        t = 'quartic'
    # Check that if quartic V_head is desired that 'a' is given
    elif t == 'quartic':
        if a_fit == None:
            print '  >  quartic specified for V_head, but no value given for "a"'
            print '     must specify using -a flag, see -h for more info - exiting...'
            sys.exit(1)

    # Remove the plateau at beginning of V(r)
    V, r = trim(fname)

    # Do the polynomial fitting of V(r)
    fit, a = fitter(r,V,order)
    dr = (r[1]-r[0]) / 2.0
    new_r = numpy.arange(min(r)-20*dr,max(r),dr)
    new_V = numpy.zeros(numpy.shape(new_r))
    V_string = "V(r) = "
    for i in range(len(a)):
        new_V += a[i]*new_r**i
        V_string += str(a[i])+'*r^'+str(i)+' + '
    r, V = new_r, new_V

    # Get the short range "head" of V(r)
    V_head, r_head = head(V, r, t, a=a_fit)

    # Get the exponentially decaying end of V(r)
    V_tail, r_tail = tail(V, r, rcut)

    # Print formula for V(r) (middle part)
    print "all other r:\n"+V_string[:-3]

    # Construct total V(r) from fitted head, middle, and tail potentials
    r_total = list(r_head) + list(r)[:list(r).index(r_tail[0])] + list(r_tail)
    V_total = list(V_head) + list(V)[:list(r).index(r_tail[0])] + list(V_tail)

    # Write V_total to file
    out = open(fname.replace('dat','fit'),'w')
    for i in range(len(r_total)):
        out.write(str(r_total[i])+' '+str(V_total[i])+'\n')
    out.close()
    

if __name__ == '__main__':
    main()
