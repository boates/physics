#!/bin/bash

# Make sure you are in the kpt directory on IRM
# i.e. in .../Cmcm/3x3x3/ or in .../I2_13/1x1x1/

T=$1
N=$2

if [ -z "$T" ]
then
    echo -e "\n usage: entropy.sh temp(K) natom_unit_cell\n"
    exit
fi
if [ -z "$N" ]
then
    echo -e "\n usage: entropy.sh temp(K) natom_unit_cell\n"
    exit
fi

# For T < 1000K
#grep -A51 Avogadro */anaddb.out | grep " "$T"\.1 " | sed s/"GPa"/" "/g | awk '{print $1,$6*6.2415097E18/6.0221415E23/'$N'.0}' > P_vs_S_0$T\K.dat

#polyfit_irm.py P_vs_S_0$T\K.dat 2 50 550
#mv -f polyfit.dat P_vs_S_0$T\K.fit
#####################################################################

grep -A100000 Avogadro */anaddb.out_THERMO | grep ""$T"\.1 " | sed s/"GPa"/" "/g | awk '{print $1,$5*6.2415097E18/6.0221415E23/'$N'.0}' > P_S_$T\K.dat
awk '{print $1,$2*'$T'.0}' P_S_$T\K.dat > P_TS_$T\K.dat

#polyfit_irm.py P_vs_S_$T\K.dat 2 60 120
#mv -f polyfit.dat P_vs_S_$T\K.fit


exit 0
