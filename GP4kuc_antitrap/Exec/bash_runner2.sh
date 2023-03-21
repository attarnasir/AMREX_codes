#!/bin/bash
#script By- Swaroop Sampad Pradhan
# Code to run :   bash bash_runner2.sh <inputfilename.in>

echo "How many cores do you want to run this script over?"
read num_cores


echo "Do you want to run with CUDA? (y/n)"
read answer2
echo $answer2
if [[ "$answer2" == "y" ]] || [[ "$answer2" == "Y" ]]; then
  export CUDA_STATUS=1
else
  export CUDA_STATUS=0
fi



if [ $CUDA_STATUS == "1" ]; then
  drivin_file=./main2d.gnu.MPI.CUDA.ex
else
  drivin_file=./main2d.gnu.MPI.ex 
fi   



make

g++ -o Replace Replace.cpp

./Replace

mpirun -np $num_cores $drivin_file $1


