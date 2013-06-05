#!/bin/bash

# in .../Cmcm/3x3x3/ or in .../I2_13/1x1x1/ etc.

N=$1

if [ -z "$N" ]
then
    echo -e "\n usage: entropy.sh natom_unit_cell\n"
    exit
fi

# For T < 1000K
#grep -A51 Avogadro */anaddb.out | grep " "$T"\.1 " | sed s/"GPa"/" "/g | awk '{print $1,$6*6.2415097E18/6.0221415E23/'$N'.0}' > P_vs_S_0$T\K.dat
#####################################################################

grep -A100000 Avogadro anaddb.out_THERMO | grep -v Avogadro | awk '{print $1,$4*6.2415097E18/6.0221415E23/'$N'.0}' > S.dat
awk '{print $1,$2*$1}' S.dat > TS.dat

exit 0
