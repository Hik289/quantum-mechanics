dir=$(cd $(dirname $0); pwd)

if [ ! -d ${dir}/exe ]
then
	mkdir ${dir}/exe
fi

echo building init...
mpic++ -o ${dir}/exe/INIT -I ${dir}/src/lib/ ${dir}/src/init.cpp

echo building cyc...
mpic++ -o ${dir}/exe/CYC -I ${dir}/src/lib/ ${dir}/src/cyc.cpp

#echo building avrg_neq_grn...
#icpc -o ${dir}/exe/AVRG_NEQ_GRN ${dir}/src/avrg_neq_grn.cpp

#echo building prep_ac...
#icpc -o ${dir}/exe/PREP_AC ${dir}/src/prep_ac_xgmx.cpp
