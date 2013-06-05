#!/usr/bin/env python
import os,sys,glob,commands

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
    except:
        print '\n usage: '+sys.argv[0]+' TRAJEC.mol(snapshot) \n'
        sys.exit(1)

    # Read header info
    natom = int(f.readline().split()[-1])
    nC    = int(f.readline().split()[-1])
    nO    = int(f.readline().split()[-1])
    f.readline()
    alat = float(f.readline().split()[-1])
    blat = float(f.readline().split()[-1])
    clat = float(f.readline().split()[-1])
    f.readline()
    f.readline()
    f.readline()

    # Read rest of file
    lines = f.readlines()
    f.close()

    # Parse line data
    typat, x, y, z, m1, m2 = [], [], [], [], [], []
    for line in lines:
        row = line.split()
        typat.append(row[0])
        x.append(row[1])
        y.append(row[2])
        z.append(row[3])
        m1.append(int(row[4])-1)
        m2.append(int(row[5])-1)

    # Loop over C atoms and print each molecule in order
    out = open('CO2_reordered.xyz','w')
    out.write(str(natom)+'\n1\n')
    for i in range(nC):
        mol  = typat[i]+' '+x[i]+' '+y[i]+' '+z[i]+'\n'
        mol += typat[m1[i]]+' '+x[m1[i]]+' '+y[m1[i]]+' '+z[m1[i]]+'\n'
        mol += typat[m2[i]]+' '+x[m2[i]]+' '+y[m2[i]]+' '+z[m2[i]]+'\n'
        out.write(mol)
    out.close()


if __name__ == '__main__':
    main()
