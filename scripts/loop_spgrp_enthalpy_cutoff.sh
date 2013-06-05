#!/bin/bash

H_ref=$1
dH=$2

if [ -z $H_ref ]
then
    echo " usage: ./loop_spgrp_enthalpy_cutoff.sh H_ref(eV) dH(meV) "
    exit 1
fi
if [ -z $dH ]
then
    echo " usage: ./loop_spgrp_enthalpy_cutoff.sh H_ref(eV) dH(meV) "
    exit 1
fi

for i in `seq 0 9`;
do
    spgrp_enthalpy_cutoff.py $H_ref $dH 0.0$i
done

exit 0
