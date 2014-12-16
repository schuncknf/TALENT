#!/bin/bash

compiler=gfortran

for file in ./*.f90 ; do
$compiler -c $file
done
for file in ./*.f ; do
$compiler -c $file
done
for file in ./*.F90 ; do
$compiler -c $file
done
# renorm-library.f renorm-main.f90 enorm-renorm-ndelta.f90 renorm-nocore.f90 renorm-phaseshift.f90 renorm-potentials.f renorm-vkrg.f90 renorm-vlowk.f90 -o exe.x
