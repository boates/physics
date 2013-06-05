#!/usr/bin/env python
"""
Introduce thermal noise to a set of configurations for a given basis.
"""
import os, sys, commands
import random, numpy

def basis_sets(basis_type):

    sc  =     numpy.array([[0.0, 0.0, 0.0]])
    
    bcc =     numpy.array([[0.0,0.0,0.0], 
                           [0.5,0.5,0.5]])
    
    fcc =     numpy.array([[0.0,0.0,0.0], 
                           [0.5,0.5,0.0], 
                           [0.5,0.0,0.5], 
                           [0.0,0.5,0.5]])
    
    cI16_JY = numpy.array([[4.0E-02,  4.0E-02,  4.0E-02], 
                           [4.6E-01,  9.6E-01,  5.4E-01], 
                           [9.6E-01,  5.4E-01,  4.6E-01],  
                           [2.9E-01,  2.9E-01,  2.9E-01],  
                           [5.4E-01,  4.6E-01,  9.6E-01], 
                           [2.1E-01,  7.1E-01,  7.9E-01],  
                           [7.1E-01,  7.9E-01,  2.1E-01],  
                           [7.9E-01,  2.1E-01,  7.1E-01],  
                           [5.4E-01,  5.4E-01,  5.4E-01],  
                           [9.6E-01,  4.6E-01,  4.0E-02],  
                           [4.6E-01,  4.0E-02,  9.6E-01],  
                           [4.0E-02,  9.6E-01,  4.6E-01],  
                           [7.9E-01,  7.9E-01,  7.9E-01],  
                           [7.1E-01,  2.1E-01,  2.9E-01],  
                           [2.1E-01,  2.9E-01,  7.1E-01],  
                           [2.9E-01,  7.1E-01,  2.1E-01]])
    
    bc8 =     numpy.array([[0.0E+00,  0.0E+00,  0.0E+00],
                           [5.0E-01,  0.0E+00,  5.0E-01],
                           [0.0E+00,  5.0E-01,  5.0E-01],
                           [5.0E-01,  5.0E-01,  0.0E+00],
                           [5.0E-01,  5.0E-01,  5.0E-01],
                           [0.0E+00,  5.0E-01,  0.0E+00],
                           [5.0E-01,  0.0E+00,  0.0E+00],
                           [0.0E+00,  0.0E+00,  5.0E-01]])
    
    diamond = numpy.array([[0.00, 0.00, 0.00],
                           [0.50, 0.50, 0.00],
                           [0.50, 0.00, 0.50],
                           [0.00, 0.50, 0.50],
                           [0.25, 0.25, 0.25],
                           [0.75, 0.75, 0.25],
                           [0.75, 0.25, 0.75],
                           [0.25, 0.75, 0.75]])
    
    cgN = numpy.array([ [0.067000, 0.067000, 0.067000], [0.567000, 0.933000, 0.567000],
                        [0.933000, 0.433000, 0.433000], [0.567000, 0.433000, 0.933000],
                        [0.500000, 0.500000, 0.500000], [0.933000, 0.433000, 0.000000],
                        [0.500000, 0.067000, 0.933000], [0.067000, 0.933000, 0.433000] ])


    return locals()[basis_type]

def main():

    try:
        x_trans = int(sys.argv[1])
        y_trans = int(sys.argv[2])
        z_trans = int(sys.argv[3])
        alat = float(sys.argv[4])/0.529177
        typat = sys.argv[5]
        nvalence = int(sys.argv[6])
        basis = sys.argv[7]
        lamda = float(sys.argv[8])
        nsteps = int(sys.argv[9])
    except:
        print '\nusage: '+sys.argv[0]+'  x_trans, y_trans, z_trans, alat (ang), typat, nvalence, basis_type, lamda, nsteps\n'
        sys.exit(0)

    a1 = numpy.array([alat,0.0,0.0])
    a2 = numpy.array([0.0,alat,0.0])
    a3 = numpy.array([0.0,0.0,alat])
    
    basis_atoms = basis_sets(basis)
    
    natom = x_trans*y_trans*z_trans*len(basis_atoms)
    rs = (alat)*(3/(4*natom*nvalence*numpy.pi))**(1.0/3.0)

    print '\n rs = '+str(rs)+'\n'
    
    sigma_x, sigma_y, sigma_z = lamda, lamda, lamda

    bravais_lattice = []
    for i in range(x_trans):
         for j in range(y_trans):
              for k in range(z_trans):
                   for l in range(len(basis_atoms)):
                        bravais_lattice.append(i*a1 + j*a2 + k*a3 + basis_atoms[l])
    
    filename = 'TRAJEC-' + '%.4f' % lamda + '.xyz'    
    outputFile = open(filename,'w')
    
    for tstep in range(nsteps):
    
         outputFile.write(repr(len(bravais_lattice)) + '\n')
         outputFile.write(repr(tstep + 1) + '\n')
         
         for triple in bravais_lattice:
             
              xcart = alat/float(x_trans)*triple[0] + random.gauss(0.0, sigma_x)
              ycart = alat/float(y_trans)*triple[1] + random.gauss(0.0, sigma_y)
              zcart = alat/float(z_trans)*triple[2] + random.gauss(0.0, sigma_z)
    
              outputFile.write(typat+' '+repr(xcart)+'   '+repr(ycart)+'   '+repr(zcart)+'\n')
                   
    outputFile.close()


if __name__ == '__main__':
    main()
