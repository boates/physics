#!/bin/bash

mkdir analysis/

grep volume OUTCAR  | grep -v ion | head -n-1 | awk '{print $5}' > volume.dat
grep TOTEN OUTCAR | grep 'free  e' | awk '{print $5}' > energy.dat
grep pressure OUTCAR | awk '{print $4/10.0}' > pressure.dat
grep pressure OUTCAR | awk '{print $4*0.00062415097}' > pressure_eVang3.dat

paste energy.dat pressure_eVang3.dat volume.dat > EPV.dat

awk '{print $1+$2*$3}' EPV.dat > enthalpy.dat

xyz_from_vasp_CO2.x

#echo 'TRAJEC.xyz' > unwrap.in
#echo $2,$2,$2 >> unwrap.in
#unwrap_PBC.x < unwrap.in
#mv TRAJEC.xyz wrapped.xyz
#mv unwrapped.xyz TRAJEC.xyz

select_snapshot.py TRAJEC.xyz 20 > log

mv -f snapshot.xyz final.xyz

nn_dist.py 200 $1 $1 $1 final.xyz final >> log
nn_dist.py 200 $1 $1 $1 TRAJEC.xyz >> log

echo 'TRAJEC_final.cnn' > angle.in
echo '200' >> angle.in
CO2_angles.x < angle.in >> log
mv CO2_angles.dat CO2_angles_final.dat
echo 'TRAJEC.cnn' > angle.in
echo '200' >> angle.in
CO2_angles.x < angle.in >> log
rm -f angle.in

mv -f *.dat *.xyz *.cnn *.hist analysis/
