#!/bin/bash

echo "TRAJEC.cnn" > in.in
echo 32 >> in.in
echo 64 >> in.in

CO2_molecule_detector.x < in.in > out.out

echo "TRAJEC.mol" > in.in

CO2_track_angles.x < in.in > out.out

echo "32" > in.in

CO2_velocity_angle_analysis.x < in.in > out.out

rm -f in.in out.out

exit 0
