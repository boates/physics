#!/usr/bin/env python
"""
Create a two_phase POSCAR file from two POSCAR files
"""
import os, sys, commands

def main():

    # Retrieve user input
    try:
       f_A = sys.argv[1]
       f_B = sys.argv[2]
    except:
        print '\n usage: '+sys.argv[0]+' POSCAR_A POSCAR_B\n'
        sys.exit(0)

    # Open the two POSCAR files
    A = open(f_A,'r')
    B = open(f_B,'r')

    # Read and process header information
    header_A = A.readline()
    header_B = B.readline()
    alat_line = A.readline()
    B.readline()
    rprim1 = A.readline()
    rprim2 = A.readline()
    rprim3 = A.readline()
    B.readline()
    B.readline()
    B.readline()
    natom_A = int( A.readline().strip() )
    natom_B = int( B.readline().strip() )
    natom = natom_A + natom_B
    direct = A.readline()
    B.readline()

    # Read in coordinates
    Ax, Ay, Az = [], [], []
    for i in range(natom_A):
        row = A.readline().split()
        Ax.append(float(row[0]))
        Ay.append(float(row[1]))
        Az.append(float(row[2]))
    Bx, By, Bz = [], [], []
    for i in range(natom_B):
        row = B.readline().split()
        Bx.append(float(row[0]))
        By.append(float(row[1]))
        Bz.append(float(row[2]))

    # Open combined POSCAR file & write header
    out = open('POSCAR_TP','w')
    out.write('two_phase: '+header_A.strip()+' '+header_B.strip()+'\n')
    out.write(alat_line)
    out.write(rprim1.replace('1','2')+rprim2+rprim3)
    out.write('  '+str(natom)+'\n')
    out.write(direct)

    # Modify and write coordinates
    for i in range(natom_A):
        out.write('  '+str(Ax[i]/2.0)+'  '+str(Ay[i])+'  '+str(Az[i])+'\n')
    for i in range(natom_B):
        out.write('  '+str((Bx[i]+1.0)/2.0)+'  '+str(By[i])+'  '+str(Bz[i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
