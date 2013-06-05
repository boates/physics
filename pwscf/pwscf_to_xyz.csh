#!/bin/csh

# give $1 as pwscf input file name, i.e. input.gamma.pw

set acell=`grep celldm $1 | head -1 | awk '{print $3}' | sed s/,//g`
set ratio=`grep celldm $1 | tail -1 | awk '{print $3}' | sed s/,//g`
set natoms=`grep nat $1 | tail -1 | awk '{print $3}' | sed s/,//g`
set atomname=`grep -A1 ATOMIC_SPECIES $1 | tail -1 | awk '{print $1}'`

grep -A$natoms ATOMIC_POSITIONS output/*.out | awk '/N/{print $2," ",$3," ",$4}' | sed s/"\-"/" \-"/g | awk -f /home/boates/software/convert_coor_TRAJEC.awk acell=$acell ratio=$ratio natoms=$natoms atomname=$atomname > TRAJEC.xyz
