#!/usr/bin/env python
"""
Loop generate_random_cell.py and submit to queue
(intended for use on clusters where optimization
is to be run)
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        natom  = sys.argv[1]
        nruns  = int(sys.argv[2])
        launch = int(sys.argv[3])
    except:
        print '\n usage: '+sys.argv[0]+' natom nruns starting_index_for_dirname\n'
        sys.exit(0)

    # Generate list of new directories for calculations
    dir_nums = range(launch,launch+nruns)
    dirs = ['0'+str(d) for d in dir_nums]

    # Create existing list of directories for checks
    old_raw = glob.glob(natom+'_atoms/0*')
    old = [o.split('/')[-1] for o in old_raw]

    # Proceed with generations and submissions
    os.chdir(natom+'_atoms')
    cwd = os.getcwd()    
    for d in dirs:
        if d not in old:
            os.mkdir(d)
            os.system('cp -rf template/* '+d)
            os.system('generate_random_cell.py -n '+natom)
            os.system('mv -f random.log '+d)
            os.system('mv -f random.POSCAR '+d+'/POSCAR')
            os.system('sleep 2s')
            os.chdir(d)
            os.system('qsub vasp.mahone2.s')
            os.chdir(cwd)
        else:
            print '\n Directory', d, 'already exists - skipping...\n'
        

if __name__ == '__main__':
    main()
