#!/bin/bash
echo 'TRAJEC.cnn' > all.in
echo 'TRAJEC_C.cnn' >> all.in
echo '0' >> all.in
echo '200' >> all.in

CO2_orientation.x < all.in > all.out

echo 'TRAJEC.cnn' > one.in
echo 'TRAJEC_C.cnn' >> one.in
echo '1' >> one.in
echo '200' >> one.in

CO2_orientation.x < one.in > one.out

rm -f all.in all.out one.in one.out
