dir=$(cd $(dirname $0); pwd)

function init_dir(){
	if [ ! -d $1 ]
	then
		mkdir $1
		echo $1 created
	else
		rm -r $1
		mkdir $1
		echo $1 refreshed
	fi
}

init_dir ${dir}/data/prep_ac_res

echo running prep_ac...
${dir}/exe/PREP_AC -d ${dir}/data
