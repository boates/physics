#!/usr/bin/env python
"""
Revised version of Isaac's xred-generator.py code
"""
import sys, random
import Numeric

def main():
    """
    Return desired atomic coordinates
    """
    try:
        num_x_trans = int(sys.argv[1])
        num_y_trans = int(sys.argv[2])
        num_z_trans = int(sys.argv[3])
        basis_type = sys.argv[4]
    except IndexError:
        print '\nusage: ' + sys.argv[0] + ' x_trans, y_trans, z_trans, basis_type\n'
        sys.exit(0)

    a = 1.0
    if basis_type == 'hcp':
        a1 = Numeric.array([1.0,0.0,0.0]) *a
        a2 = Numeric.array([-0.5,3**0.5/2.0,0.0]) *a
        a3 = Numeric.array([0.0,0.0,(8.0/3.0)**0.5]) *a
    else:
        a1 = Numeric.array([a,0.0,0.0])
        a2 = Numeric.array([0.0,a,0.0])
        a3 = Numeric.array([0.0,0.0,a])
    
    sc  = Numeric.array([ [0.0,0.0,0.0] ])
    hcp = Numeric.array([ [0.0,0.0,0.0] ])
    bcc = Numeric.array([ [0.0,0.0,0.0], [0.5,0.5,0.5] ])
    fcc = Numeric.array([ [0.0,0.0,0.0], [0.5,0.5,0.0], [0.5,0.0,0.5], [0.0,0.5,0.5] ])
    cI16_JY = Numeric.array([ [4.0E-02,  4.0E-02,  4.0E-02], [4.6E-01,  9.6E-01,  5.4E-01],
                              [9.6E-01,  5.4E-01,  4.6E-01], [2.9E-01,  2.9E-01,  2.9E-01],
                              [5.4E-01,  4.6E-01,  9.6E-01], [2.1E-01,  7.1E-01,  7.9E-01],
                              [7.1E-01,  7.9E-01,  2.1E-01], [7.9E-01,  2.1E-01,  7.1E-01],
                              [5.4E-01,  5.4E-01,  5.4E-01], [9.6E-01,  4.6E-01,  4.0E-02],
                              [4.6E-01,  4.0E-02,  9.6E-01], [4.0E-02,  9.6E-01,  4.6E-01],
                              [7.9E-01,  7.9E-01,  7.9E-01], [7.1E-01,  2.1E-01,  2.9E-01],
                              [2.1E-01,  2.9E-01,  7.1E-01], [2.9E-01,  7.1E-01,  2.1E-01] ])
    tetragonal_nitrogen = Numeric.array([ [0.098,0.098,0.0], [0.902,0.902,0.0],
                                          [0.402,0.598,0.5], [0.598,0.402,0.5] ])
    Nepsilon = Numeric.array([  [0.0477444, 0.0394043, 0.949856], [0.952256, 0.960596, 0.0501445],
                                [0.531101, 0.438327, 0.442201],   [0.435612, 0.359518, 0.54249],
                                [0.172472, 0.502245, 0.925462],   [0.22403, 0.399988, 0.031043],
                                [0.648548, 0.489154, 0.910909],   [0.700106, 0.386897, 0.0164902],
                                [0.0960269, 0.157856, 0.472166],  [0.970638, 0.227801, 0.589518],
                                [0.0374346, 0.633818, 0.458946],  [0.912046, 0.703763, 0.576298],
                                [0.540886, 0.00789884, 0.7523],   [0.446018, 0.90098, 0.706547],
                                [0.606758, 0.94395, 0.272419],    [0.51189, 0.837031, 0.226665]   ])
    cgN = Numeric.array([ [0.067000, 0.067000, 0.067000], [0.567000, 0.933000, 0.567000],
                          [0.933000, 0.433000, 0.433000], [0.567000, 0.433000, 0.933000],
                          [0.500000, 0.500000, 0.500000], [0.933000, 0.433000, 0.000000],
                          [0.500000, 0.067000, 0.933000], [0.067000, 0.933000, 0.433000] ])

    Nchain = Numeric.array([ [0.000000, 0.000000, 0.000000], [0.000000, 0.500000, 0.1300000],
                             [0.500000, 0.000000, 0.500000], [0.500000, 0.500000, 0.6300000] ])

#    Nchain = Numeric.array([ [0.000000, 0.000000, 0.000000], [0.000000, 1.275000, 0.680000],
#                             [1.055000, 0.000000, 2.616500], [1.055000, 1.275000, 3.297000] ])

    
    basis_atoms = locals()[basis_type]
    
    bravis_lattice = []
    
    for i in range(num_x_trans):
         for j in range(num_y_trans):
              for k in range(num_z_trans):
                   for l in range(len(basis_atoms)):
                        bravis_lattice.append(i*a1 + j*a2 + k*a3 + basis_atoms[l])
    
    filename = 'xred.dat'
    outputFile = open(filename,'w')
    
    for triple in bravis_lattice:
             
         xcart = a/float(num_x_trans)*(triple[0])
         ycart = a/float(num_y_trans)*(triple[1])
         zcart = a/float(num_z_trans)*(triple[2])
    
         outputFile.write(' % .8f'  %  xcart + '  % .8f' % ycart + '  % .8f' %  zcart + '\n')
                   
    outputFile.close()

if __name__ == '__main__':
    main()
