numThd=12
dir=$(cd $(dirname $0); pwd)

if [ ! -d ${dir}/data/log ]
then
	mkdir ${dir}/data/log
	echo data/log/ created
else
	rm -r ${dir}/data/log
	mkdir ${dir}/data/log
	echo data/log/ refreshed
fi

if [ ! -d ${dir}/data/stts ]
then
	mkdir ${dir}/data/stts
	echo data/stts/ created
else
	rm -r ${dir}/data/stts
	mkdir ${dir}/data/stts
	echo data/stts/ refreshed
fi

if [ ! -d ${dir}/data/equ_check ]
then
	mkdir ${dir}/data/equ_check
	echo data/equ_check/ created
else
	rm -r ${dir}/data/equ_check
	mkdir ${dir}/data/equ_check
	echo data/equ_check/ refreshed
fi

if [ ! -d ${dir}/data/res ]
then
	mkdir ${dir}/data/res
	echo data/res/ created
else
	rm -r ${dir}/data/res
	mkdir ${dir}/data/res
	echo data/res/ refreshed
fi

if [ ! -d ${dir}/data/neq ]
then
	mkdir ${dir}/data/neq
	echo data/neq/ created
else
	rm -r ${dir}/data/neq
	mkdir ${dir}/data/neq
	echo data/neq/ refreshed
fi

echo running init...
mpirun -n ${numThd} ${dir}/exe/INIT -d ${dir}/data
