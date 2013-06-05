#!/usr/bin/env python
"""
Calculate the strucutre factor from g(r)
"""
import sys, Numeric
import FFT

def caluclate_SF(hr,r,natom,alat,verbose=False):
    """
    Perform Fourier Transform and calculate the strucure factor
    """
    # Structure Factor Caluclation
    N  = float(len(r))   # get # of data points
    dr = r[1] - r[0]     # get distance spacing
    R  = N*dr            # define the period (total time)
    dQ = 1./R            # define frequency step
    Q  = Numeric.arange(N,typecode=Numeric.Float)*dQ
    H  = FFT.fft(hr)     # calculate the fft
    print N,dr,R,dQ

    SF = 1.0 + (natom/alat**3)*H.real#(Numeric.conjugate(H)*H).real
    SF = Numeric.array( list(SF)[:int(len(SF)/2.0)] )

    if verbose == True:
        hr2 = FFT.inverse_fft(H).real
        out = open('hr2.dat','w')
        for i in range(len(hr2)):
            out.write(str(r[i])+' '+str(hr2[i])+' '+str(Q[i])+' '+str(H.real[i])+'\n')
        out.close()

    return SF, Q                        

def smoothing(hr,r,verbose=False):
    """
    Smoothing algorithm for h(r)
    """
    if hr[-1] > hr[-2] and hr[-1] > 0.0 and hr[-2] > 0.0:
        p = hr[-2]
        prev = []
        k = 0
        while p > 0.0:
            p = hr[-2-k]
            if p > 0.0:
                prev.append(p)
            k += 1
        dr = r[-1] - r[-2]
        prev = [hr[-1],(hr[-1]+hr[-2])/2.0] + list(prev)
        prev[-1] = (prev[-2]+prev[-1])/2.0
        prev.append(prev[-1]/2.0)
        prev.append(prev[-1]/2.0)
        prev.append(prev[-1]/2.0)
        prev.append(prev[-1]/2.0)
        prev.append(0.0)
        r_add = [r[-1]+dr]
        for i in range(len(prev)):
            r_add.append(r_add[i]+dr)

    hr = Numeric.array( list(hr) + list(prev) )
    r = Numeric.array( list(r) + list(r_add) )

    # Write smoothed h(r) to file
    if verbose == True:
        out = open('hr.dat','w')
        for i in range(len(hr)):
            out.write(str(r[i])+' '+str(hr[i])+'\n')
        out.close()

    return hr, r

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1])
	lines = f.readlines()
	f.close()
        natom = int(sys.argv[2])
        alat = float(sys.argv[3])
        smooth = sys.argv[4]
        if smooth == 'y': smooth = True
        elif smooth == 'n': smooth = False
    except:
        print '\nusage: '+sys.argv[0]+' RDF.dat, natom, alat (ang), smooth (y/n)\n'
	sys.exit(0)

    # Read in data and format
    r, gr, integral = [], [], []
    for line in lines:
        row = line.split()
	r.append(float(row[0]))
	gr.append(float(row[1]))
        integral.append(float(row[2]))
    r.pop(-1)
    gr.pop(-1)
    integral.pop(-1)
    r = Numeric.array(r)
    gr = Numeric.array(gr)
    integral = Numeric.array(integral)

    # h(r) = g(r) - 1, calculation and alteration
    if smooth:
        hr, r = smoothing(gr-1,r,verbose=False)
    else: hr = gr - 1

    # Calculate the structure factor
    SF, Q = caluclate_SF(hr,r,natom,alat,verbose=False)

    # Write the structure factor to output file
    out = open('SF.dat','w')
    for i in range(int(len(Q))):
        out.write(str(Q[i])+' '+str(SF[i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
