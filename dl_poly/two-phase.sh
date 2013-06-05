#!/bin/bash

# Script to perform automatic two-phase simulations over a range of specified temperatues

temps="\
400 \
410 \
420 \
430 \
440 \
450 \
460 \
470 \
480 \
490 \
500 \
510 \
520 \
530 \
540 \
550 \
560 \
570 \
580 \
590 \
"

for t in $temps; do

    mkdir $t
    cp 0600K_preparation/REV* $t
    cp third_stage_inputs/* $t
    cd $t

    /home/boates/software/dl_poly/copy
    sed s/TEMPERATURE/$t.0/g CONTROL.template > CONTROL
    qsub dl.par

    cd ../

done

exit 0
    
