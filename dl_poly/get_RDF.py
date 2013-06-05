#!/usr/bin/env python

# Use for DL_POLY RDFDAT files

import os, sys, commands

def writelines(fin, fout, n):
    """ Writes next n lines from fin to fout. """
    for i in range(int(n)):
        fout.write(fin.readline())

def inputstr(comment, defaultvalue):

    # Gets a string from the user with a default value
    print comment, "(default = " + defaultvalue + "): "
    temp = sys.stdin.readline().strip()
    if temp != '':
        return temp
    else:
        return defaultvalue                                                

try:
    atom1 = sys.argv[1]
    atom2 = sys.argv[2]
except IndexError:
    print '\n--Atoms not specified--\n'
    atom1 = inputstr('Give atom1', 'all')
    atom2 = inputstr('Give atom2', 'all')

try:
    plot = sys.argv[3]
except IndexError:
    plot = inputstr('Would you like to plot the results? (y/n),', 'y')
    print '...'
    
try:
    f = open('RDFDAT','r')
except IOError:
    print '\nRDFDAT does not exist--stopping\n'
    sys.exit(0)
    
header = f.readline()
nrdf, npt = f.readline().split()
fnames = []
atom_choices = []

for i in range(int(nrdf)):

    chosen_atoms = atom1+'_'+atom2
    ATOMS = f.readline()
    atoms = ATOMS.split()[0]+'_'+ATOMS.split()[1]
    fname = 'RDF_'+atoms+'.dat'
    atom_choices.append('atom1 = '+ATOMS.split()[0]+'\tatom2 = '+ATOMS.split()[1])

    out = open(fname,'w')

    if chosen_atoms == atoms:
        writelines(f,out,npt)
        fnames.append(fname)
        break
    elif atom1 == 'all' and atom2 == 'all':
        writelines(f,out,npt)
        fnames.append(fname)
    elif atom2 == 'all' and atom1 == atoms.split('_')[0]:
        writelines(f,out,npt)
        fnames.append(fname)
    elif atom1 == 'all' and atom2 == atoms.split('_')[1]:
        writelines(f,out,npt)
        fnames.append(fname)
    else:
        for i in range(int(npt)):
            r.readline()
    out.close()

f.close()

if fnames == []:
    print '\n--No RDF was found in RDFDAT with the given atom specifications--'
    nchoice = str(len(atom_choices))
    display_choices = inputstr('Display the '+nchoice+' possible choices? (y/n)','y')
    if display_choices == 'y':
        for choice in atom_choices: print choice
    plot = 'n'

if plot == 'y':
    p = open('plot_RDF.tmp','w')
    p.write('plot '+fnames.__repr__()[1:-1])
    p.close()
    os.system('plot_RDF.sh')
