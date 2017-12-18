#!/bin/bash

ARRAY=(1380)
NUC=(ADP)
PROT=(prot deprot)
CD=`pwd`

function MakePDBFiles {
#Take 1 ns snapshots from 35 - 85 ns (50)
DIR=/home/birac/Desktop/LRA/propka/XTCs
XTC=$DIR/${NUC[j]}.100.${ARRAY[i]}_MD_${PROT[k]}.xtc
TPR=$DIR/${NUC[j]}.100.${ARRAY[i]}_MD_${PROT[k]}.tpr
echo Protein | $GMX51 trjconv -f ${XTC} -s ${TPR} -b 35000 -e 85000 -skip 10 -sep -o ${NUC[j]}.${ARRAY[i]}.${PROT[k]}.pdb >& trj
# ADP.1380.prot45.pdb
sed '/REMARK/d' *pdb -i
sed '/TITLE/d' *pdb -i
sed '/CRYST/d' *pdb -i
sed '/MODEL/d' *pdb -i
}

function RunProPKA {
pdbs=(*pdb)
for (( m=0; m<${#ARRAY[@]}; m++ )); do
	file=${pdbs[m]::-4}
	propka31 -f ${file}.pdb > ${file}.out
done
}

function PrepLRA {
#pKa = 1/2 * ({pKa(c)}p + {pKa(c)}d)
pdbs=(*pdb)
rm -f ${NUC[j]}.${ARRAY[i]}.${PROT[k]}.pkas
for (( m=0; m<${#ARRAY[@]}; m++ )); do
        res=`echo "${ARRAY[i]} - 1353" | bc`
        pka=`grep "^LYS  $res F" ${file}.out | awk '{print $4}'`
	echo $pka >> ${NUC[j]}.${ARRAY[i]}.${PROT[k]}.pkas
done
}

function DoLRA {
echo hello
}

for (( j=0; j<${#NUC[@]}; j++ )); do
	mkdir -p ${NUC[j]}
	cd ${NUC[j]}
	CD2=`pwd`
	for (( i=0; i<${#ARRAY[@]}; i++ )); do
		mkdir -p ${ARRAY[i]}
		cd $CD2
		cd ${ARRAY[i]}
		for (( k=0; k<${#PROT[@]}; k++ )); do
			MakePDBFiles
			RunProPKA
			PrepLRA
		done
i#		DoLRA
	done
done
