#!/usr/bin/env python
"""
gibbs_extract.py
Author: Brian Boates

Extracts all the thermodynamic variables fro rs directories
needed to calculate the Gibbs free energy
"""
import os, sys, commands, glob, numpy

def main():

    # Retrieve user input
    try:
        sFile = open(sys.argv[1],'r')
        natom = int(sys.argv[2])
        eShift = float(sys.argv[3])
        pShift = float(sys.argv[4])
    except:
        # Energy and pressure corrections are the amount to ADD to the detected averages
        # This should be the average of (E_more_kpt/atom - E_less_kpt/atom) (similar for P)
        # Hence, if more kpts leads to a "favorable" shift, the number should given should be negative
        print '\n usage: '+sys.argv[0]+' S_per_atom_vs_P.fit natom E_correction(eV/atom) P_correction(GPa)\n'
        sys.exit(0)

    # rs file
    os.system("grep Final */*/energy.blocker | sed s/'\/analysis'/' '/g | awk '{print $1}' > rs")

    # E file
    os.system("grep Final */*/energy.blocker | awk '{print $4/"+str(natom)+"+"+str(eShift)+"}' > E")

    # P file
    os.system("grep Final */*/pressure.blocker | awk '{print $4+"+str(pShift)+"}' > P")

    os.system("paste P E > pe.dat")

    # V file
    f = open('rs','r')
    rs = f.readlines()
#    rs = ['1.85']
    f.close()
    vFile = open('V','w')
    for i in range(len(rs)):
        alat = float(commands.getoutput("head -5 "+rs[i].strip()+"/POSCAR | tail -4 | head -1"))
        ax, ay, az = commands.getoutput("head -5 "+rs[i].strip()+"/POSCAR | tail -3 | head -1").split()
        a = numpy.array([float(ax), float(ay), float(az)])*alat
        bx, by, bz = commands.getoutput("head -5 "+rs[i].strip()+"/POSCAR | tail -2 | head -1").split()
        b = numpy.array([float(bx), float(by), float(bz)])*alat
        cx, cy, cz = commands.getoutput("head -5 "+rs[i].strip()+"/POSCAR | tail -n-1").split()
        c = numpy.array([float(cx), float(cy), float(cz)])*alat

        v = numpy.dot( a, numpy.cross(b,c) ) / natom

        vFile.write(str(v)+'\n')
    vFile.close()

    # T file
    T = float(commands.getoutput("grep TEBEG */INCAR | head -1 | awk '{print $3}'"))
    os.system("grep Final */*/pressure.blocker | sed s/'\/analysis'/' '/g | awk '{print $1-$1+"+str(T)+"}' > T")

    # S file
    lines = sFile.readlines()
    sFile.close()
    pfit, sfit = [], []
    for line in lines:
        row = line.split()
        pfit.append(float(row[0]))
        sfit.append(float(row[1]))

    p = commands.getoutput("head -100 P").split()
    for i in range(len(p)):
        p[i] = float(p[i])

#    if max(p) > max(pfit) or min(p) < min(pfit):
#        print '\n Pressure range for the fit entropy does not encompass all pressures found - exiting...\n'
#        sys.exit(0)

    s = open('S','w')
    for i in range(len(p)):
        dp = []
        for j in range(len(pfit)):
            dp.append(abs(pfit[j]-p[i]))
        k = dp.index(min(dp))
        s.write(str(sfit[k])+'\n')
    s.close()

    # S file from VDOS integration
    os.system('head */analysis/entropy.dat | grep -v analysis | grep "\." > S_int')

    # gibbs.dat, G, and pg.dat files
    os.system("paste rs E P V T S | awk '{print $1,$2,$3,$4,$5,$6,($2+$3*$4*0.0062415097-$5*$6)}' > gibbs.dat")
    os.system("awk '{print $7}' gibbs.dat > G")
    os.system("paste P G > pg.dat")
    os.system("polyfit.py pg.dat 2 10 150")
    os.system("mv polyfit.dat pg.fit")

    # gibbs_int.dat, G_int, and pg_int.dat files
    os.system("paste rs E P V T S_int | awk '{print $1,$2,$3,$4,$5,$6,($2+$3*$4*0.0062415097-$5*$6)}' > gibbs_int.dat")
    os.system("awk '{print $7}' gibbs_int.dat > G_int")
    os.system("paste P G_int > pg_int.dat")
    os.system("polyfit.py pg_int.dat 2 10 150")
#    os.system("polyfit.py pg_int.dat 2 10 700")
    os.system("mv polyfit.dat pg_int.fit")

    # PV, TS, and TS_int files
    os.system("paste rs P T S | awk '{print $2,$3*$4}' > TS")
    os.system("paste rs P T S_int | awk '{print $2,$3*$4}' > TS_int")

    # H and ph.dat files
    os.system("paste E P V | awk '{print $1+$2*$3*0.0062415097}' > H")
    os.system("paste P V | awk '{print $1,$1*$2*0.0062415097}' > PV")
    os.system("paste P H > ph.dat")
    os.system("paste P S_int > ps_int.dat")

if __name__ == '__main__':
    main()
