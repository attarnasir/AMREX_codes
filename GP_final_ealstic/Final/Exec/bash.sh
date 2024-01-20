#!/bin/bash

#make clean

#make phasecount=3

make

g++ -o Replace Replace.cpp

./Replace input_multiphase.in Filling_multiphase.in
#./Replace Input.in Filling.in
#./Replace Input_TP_13k.in Filling_TP.in
#./Replace Input_3d_3p.in Filling_3d_3p.in
#./Replace Input_3d_2p.in Filling_3d_2p.in
#./Replace Input_tdb_new_NiAlMo.in Filling.in

mpirun -np 4 ./main2d.gnu.MPI.ex input2.in

#mpirun -np 4 ./main2d.gnu.MPI.ex input13k.in

#mpirun -np 4 ./main3d.gnu.MPI.ex input2.in
