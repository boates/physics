#!/usr/bin/env python
"""
Report space groups of of structures that lay within a given range
of a given comparison enthalpy.
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        H_ref = float(sys.argv[1])
        dH    = float(sys.argv[2])
        tol   = sys.argv[3]
    except:
        print '\n usage: '+sys.argv[0]+'  H_ref(eV/atom)  dH(meV/atom)  tolerance(e.g. 0.01)\n'
        sys.exit(0)

    # Detect present directories
    dirs_raw = commands.getoutput("tail */OUTCAR | grep -B10 Voluntary | grep OUTCAR | awk '{print $2}'").split()
    dirs = [d.split('/')[0] for d in dirs_raw]
    if 'cannot' in dirs:
        print '\n No completed data detected - exiting... \n'
        sys.exit(0)

    # Perform analysis
    cwd = os.getcwd()
    out = open('competing_structures_tol'+tol+'.dat','w')
    out.write('Using a tolerance of '+tol+' and H_ref='+str(H_ref)+' (eV/atom)\n')
    out.write('# dir\tspgrp\tH(eV/atom)\tdH(meV/atom)\n')
    rej = open('rejected_structures_tol'+tol+'.dat','w')
    rej.write('Using a tolerance of '+tol+', H_ref = '+str(H_ref)+' eV/atom, & dH = '+str(dH)+' meV/atom\n')
    rej.write('# dir\tspgrp\tH(eV/atom)\tdH(meV/atom)\n')
    all = open('all_structures_tol'+tol+'.dat','w')
    all.write('Using a tolerance of '+tol+' and H_ref='+str(H_ref)+' (eV/atom)\n')
    all.write('# dir\tspgrp\tH(eV/atom)\tdH(meV/atom)\n')
    for d in dirs:
        os.chdir(d)
        if 'CONTCAR' in glob.glob('*'):
            tail = commands.getoutput('tail -n-1 OUTCAR')
            if 'Voluntary' in tail.split():
                os.system('run_findsym.py CONTCAR '+tol)
                spgrp = commands.getoutput("grep 'Space Group' findsym_tol"+tol+".out | awk '{print $5}'")
                H = commands.getoutput('enthalpy_from_vasp_relax.py').split()[-2]
                diff = (float(H) - H_ref) * 1000.0
                if (abs(diff) < dH) or (diff < 0):
                    out.write(d+'\t'+spgrp+'\t'+H+'\t'+str(diff)+'\n')
                    out.flush()
                else:
                    rej.write(d+'\t'+spgrp+'\t'+H+'\t'+str(diff)+'\n')
                    rej.flush()
                all.write(d+'\t'+spgrp+'\t'+H+'\t'+str(diff)+'\n')
                all.flush()
        os.chdir(cwd)
    out.close()
    rej.close()
    all.close()


if __name__ == '__main__':
    main()
