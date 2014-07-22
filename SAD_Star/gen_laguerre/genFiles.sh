#!/bin/bash

# This script will generate the files containing abscissa and weight of a quadrature

# Parameters:
#############################
quadName="gauss-laguerre"
orderMin=1400
orderMax=2000
#############################

for i in `seq $orderMin $orderMax` ; do
  echo $i
  name=$quadName"_n"$i
  ./a.out $i 0. 0. 1 $name > tmp
  rm $name"_r.txt"
done

rm tmp
echo ""
echo "Files generated !"
echo ""

