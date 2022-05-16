#!/bin/bash

#########################check folders########################
dir=$(cd $(dirname $0); pwd)

if [ ! -d ${dir}/exe ]
then
	mkdir ${dir}/exe
fi

if [ ! -d ${dir}/data ]
then
	mkdir ${dir}/data
fi

if [ ! -d ${dir}/data/log ]
then
	mkdir ${dir}/data/log
fi


###########################compiling###########################
g++ -o ${dir}/exe/MAIN -I ${dir}/src/lib/ ${dir}/src/main.cpp



############################running############################
valgrind --tool=memcheck --leak-check=full  ${dir}/exe/MAIN -d ${dir}/data
