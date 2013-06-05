#!/bin/bash\
### Script to do the analyzing of Lithium melt curve simulations
### acell required as command line argument (angstroms)
acell=$1
temps=" \
0100K \
0200K \
0300K \
0400K \
0450K \
0500K \
0510K \
0520K \
0530K \
0540K \
0550K \
0560K \
0570K \
0580K \
0590K \
0600K \
0650K \
0700K \
0750K \
0800K \
0850K \
0900K \
0950K \
1000K \
"

for n in $temps; do

    # Enter appropriate directory
    cd $n

    # Calculate everything you want
    get_RDF.py all all n
    get_pvt.py
    convert_HISTORY_xyz.x
    # Unwrapping requires input: make a file
    echo HISTORY.xyz > unwrap_PBC.in
    echo $acell >> unwrap_PBC.in
    unwrap_PBC.x < unwrap_PBC.in
    echo unwrapped.xyz > msd.in
    msd.x < msd.in
    
    # Remove temporary input files
    rm unwrap_PBC.in, msd.in

    # Rename data appropriately
    mv RDF_LI_LI.dat $n.rdf
    mv pvt.dat $n.pvt
    mv HISTORY.xyz $n.xyz
    mv unwrapped.xyz $n.unwrapped.xyz
    mv msd.dat $n.msd

    # Exit the directory
    cd ../

done

exit 0
