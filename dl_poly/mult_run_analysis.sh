#! /bin/bash

# Script to enter a bunch of run[#] directories
# and run run_STATIS.sh, diffusion_from_statis.py
# and get_RDF.py

runs="
1 \
2 \
3 \
4 \
5 \
6 \
7 \
8 \
9 \
10 \
"

for run in $runs; do
    
    cd run$run
    run_STATIS.sh
    diffusion_from_statis.py
    get_RDF.py all all n
    cd ../

done

averaging_diffusion.py 10

exit 0
