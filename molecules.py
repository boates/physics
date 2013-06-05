#!/usr/bin/env python
"""
molecules.py v2.0
Author: Brian Boates

usage: molecules.py molecules.in

Loop over an xyz file and group atoms into connected polymer lists
& output a histogram of molecular species probabilities

Note: atom indexing begins at ZERO. i.e. a 64 atom system is indexed from
atom 0 to 63, NOT 1 to 64

New features!
1. input file rather than command line, format is:
   |  TRAJEC.xyz
   |  alat (Ang)
   |  A B C ... (species names)
   |  A-B AB_cut (species pair and cut-off in Ang)
   |  A-A AA_cut
   |  A-C AC_cut
   |  B-B BB_cut
   |  B-C BC_cut
   |  C-C CC_cut
   |  ... etc.

2. multiple cut-offs allowed (see input file template above for examples)

"""
import os, sys, commands, glob, numpy

def factorial(n):
    """
    Return the factorial of an integer n
    """
    f = 1
    for i in range(1,n+1):
        f *= i
    return f

                        
def choose(m,n):
    """
    Return 'm choose n'
    """
    return factorial(m) / (factorial(n)*factorial(m-n))


def common_value(a,b):
    """
    a, b: two lists of integers
    returns True if a and b share a common value, False otherwise
    """
    c = False
    for v in a:
        if v in b:
            c = True
            break
    return c


def count_atoms(molecules):
    """
    molecules: a list of sublists, each sublist containing integer atomic indices
    returns integer value equal to the combined length of all sublists
    """
    sum = 0
    for m in molecules:
        sum += len(m)
    return sum
        

def remove_duplicates(L):
    """
    L: a list of positive integers
    return: a new list where repeated values have been removed from the list
    """
    new = []
    for i in range(len(L)):
        if L[i] not in new:
            new.append(L[i])
    return new
    

def merge(molecules):
    """
    molecules: a list of sublists, each sublist containing integer atomic indices
    return: a list of merged sublists, merging criteria is if the same index appears
            in multiple sublists, they are merged into a new sublist with redundant
            indicies removed using the remove_duplicates() function.
    """
    natom = max( max(molecules) ) + 1   # +1 comes from zero to natom-1 indexing scheme

    # While there are still duplicates
    while count_atoms(molecules) > natom:

        # Loop over molecules
        for i in range(len(molecules)):
            for j in range(len(molecules)):

                # Avoid self-checking and reassigned containers
                if i != j and molecules[i] != [] and molecules[j] != []:

                    # If molecules i & j share a common atom
                    if common_value(molecules[i], molecules[j]):

                        # Connect molecule j to i and remove duplicates
                        molecules[i] = remove_duplicates(molecules[i]+molecules[j])

                        # Reassign molecule j to an empty list
                        molecules[j] = []

        # Remove empty lists
        for i in range( molecules.count([]) ):
            molecules.pop( molecules.index([]) )

    return molecules


def pbc_round(input_value):
    i = int(input_value)
    if (abs(input_value-i) >= 0.5):
        if (input_value > 0): i+= 1
        if (input_value < 0): i-= 1
    return i


def main():

    # Retrieve user input
    try:
        fin  = open(sys.argv[1],'r')  # open given input file
        fxyz = fin.readline().strip()  # xyz filename is first line in input
        xyz  = open(fxyz,'r')  # open xyz file
        alat = float(fin.readline())  # lattice constant (ang) is 2nd input line
        types = fin.readline().split()  # atom species separated by spaces is 3rd input line
        ntypes = len(types)  # determine ntypes from types
        npairs = choose(ntypes,2)+ntypes  # determine number of pairs

        # The "rest" of the input file should contain the cut-offs in angstroms - check
        rest = fin.readlines()
        if len(rest) != npairs:
            print "\n Wrong number of cut-offs provided --- expected "+napirs+", found "+len(rest)
            print "\n --- Revise input file\n"
            sys.exit(1)
        else:
            # The cutoffs dictionary contains all 'A-B' (and 'B-A') pairs and their cutoff in ang
            cutoffs = {}
            for line in rest:
                row = line.split()
                atom1 = row[0].split('-')[0]  # atom1 species
                atom2 = row[0].split('-')[1]  # atom2 species
                cutoffs[atom1+'-'+atom2] = float(row[1])  # add atom1-atom2 to dict.
                if atom1 != atom2:
                    cutoffs[atom2+'-'+atom1] = float(row[1])  # add atom2-atom1 to dict.
        fin.close()
    except:
        print '\n usage: '+sys.argv[0]+'  molecules.in\n'
        sys.exit(1)

    # Determine nsteps and natom
    natom   =  int(commands.getoutput("head -1 "+fxyz))
    nsteps  =  int( commands.getoutput("wc -l "+fxyz+" | awk '{print $1}'")) / (natom+2)
    nmolecules = 0
    
    # Create a dictionary to sum occurences of all different molecular species
    tracker = {}

    # Perform the analysis
    print "Progress =",
    for i in range(nsteps):

        # Print progress every 100 steps
        if i % 100 == 0:
            print str(round(i/float(nsteps)*100.))+"% -",

        # Read in the current snapshot
        xyz.readline() # natom
        xyz.readline() # tstep
        typat, x, y, z = [], [], [], []
        for j in range(natom):
            row = xyz.readline().split()
            typat.append(row[0])
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(float(row[3]))

        # Create molecules list
        molecules = []
        for j in range(natom):

            # Append a new molecule for the current atom
            molecules.append( [j] )

            # Find atoms within cutoff distance of atom j
            for k in range(natom):

                # Be sure not to check j with itself
                if j != k:

                    # Calculate the distance between atoms j & k
                    dx = x[k] - x[j]
                    dy = y[k] - y[j]
                    dz = z[k] - z[j]
                    dx -= alat * pbc_round( dx/alat )
                    dy -= alat * pbc_round( dy/alat )
                    dz -= alat * pbc_round( dz/alat )
                    d = (dx**2 + dy**2 + dz**2)**0.5

                    # If the distance is less than the cutoff for given pair
                    pair = typat[j]+'-'+typat[k]
                    if d <= cutoffs[pair]:

                        # Append atom k to atom j's molecule
                        molecules[j].append(k)

        #============================================#
        # We now have our initial molecule groupings #
        #--------------------------------------------#
        #     Now we must merge connected groups     #
        #============================================#

        # Merging and duplicate removing takes place in defined functions
        merged = merge(molecules)

        #========================================#
        # We now have our final merged groupings #
        #----------------------------------------#
        #     Now we can see what we have :)     #
        #========================================#

        # Loop over molecules
        for molecule in molecules:

            # Get atomic species for molecule
            s = ''
            for atom in molecule:
                s += typat[atom]
            species = ''
            for atom in types:
                species += atom+str(s.count(atom))

            # If this is the first occurence of this species
            if species not in tracker.keys():
                tracker[species] = 1

            # Otherwise the species is already present
            else:
                tracker[species] += 1

        # Keep a count of total number of molecules found for normalization
        nmolecules += len(molecules)

    xyz.close()

    print "100%"

    #=================================#
    #      Analysis is complete!      #
    #---------------------------------#
    # Now normalize and write to file #
    #=================================#

    # Sort molecules by fraction
    sorted = [ [t[1],t[0]] for t in tracker.items() ]
    sorted.sort()
    sorted = sorted[::-1]

    # Write molecular fractions to output file
    out1 = open('molecules.dat','w')
    out2 = open('molecules.avg','w')
    out1.write('# species, fraction\n')
    out2.write('# species, Nfound/nsteps\n')
    for pair in sorted:
        out1.write(pair[1]+'\t'+str( float(pair[0]) / nmolecules )+'\n')
        out2.write(pair[1]+'\t'+str( float(pair[0]) / nsteps )+'\n')
    out1.close()
    out2.close()


if __name__ == '__main__':
    main()
