numThd=12
dir=$(cd $(dirname $0); pwd)

echo running mesr...
mpirun -n ${numThd} ${dir}/exe/CYC -d ${dir}/data
