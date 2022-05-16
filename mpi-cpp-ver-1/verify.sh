echo building
mpic++ -o VRFY head.h head.cpp verify.cpp func.cpp para.cpp

echo running
mpirun -n 8 ./VRFY
