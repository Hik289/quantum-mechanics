echo compiling
make
echo running
mpirun -n 8 ./TQ
