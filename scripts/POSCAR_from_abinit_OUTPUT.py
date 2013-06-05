#!/usr/bin/env python
"""
Create POSCAR file from abinit OUTPUT file
"""
import os, sys, commands, glob

def main():

    natom = commands.getoutput("grep natom OUTPUT | tail -n-1 | awk '{print $2}'")

    alat = float(commands.getoutput("grep acell OUTPUT | awk '{print $2}'"))
    blat = float(commands.getoutput("grep acell OUTPUT | awk '{print $3}'"))
    clat = float(commands.getoutput("grep acell OUTPUT | awk '{print $4}'"))

    try:
        ax, ay, az = commands.getoutput("grep rprim OUTPUT").split()[1:]
        bx, by, bz = commands.getoutput("grep -A1 rprim OUTPUT | tail -n-1").split()
        cx, cy, cz = commands.getoutput("grep -A2 rprim OUTPUT | tail -n-1").split()
        ax, ay, az = float(ax)*alat, float(ay)*blat, float(az)*clat
        bx, by, bz = float(bx)*alat, float(by)*blat, float(bz)*clat
        cx, cy, cz = float(cx)*alat, float(cy)*blat, float(cz)*clat
    except:
        ax, by, cz = alat, blat, clat
        ay, az, bx, bz, cx, cy = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        
    a0 = 0.5291772
    header = commands.getoutput("grep 'space group' OUTPUT")

    out = open('POSCAR','w')
    out.write(header+'\n')
    out.write('    1.00000000\n')
    out.write('     '+str(ax*a0)+'  '+str(ay*a0)+'  '+str(az*a0)+'\n')
    out.write('     '+str(bx*a0)+'  '+str(by*a0)+'  '+str(bz*a0)+'\n')
    out.write('     '+str(cx*a0)+'  '+str(cy*a0)+'  '+str(cz*a0)+'\n')
    out.write('  '+natom+'\n')
    out.write('Direct\n')
    out.close()

    os.system("grep xred OUTPUT | awk '{print $2/1.0,$3/1.0,$4/1.0}' >> POSCAR")
    os.system("grep xred -A"+str(int(natom)-1)+" OUTPUT | grep -v xred | awk '{print $1/1.0,$2/1.0,$3/1.0}' >> POSCAR")

    
if __name__ == '__main__':
    main()
