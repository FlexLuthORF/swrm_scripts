#!/bin/bash
#Iter over list of sample numbers and submit changeo script to qsub
#requires generation fasta AND fmt7 file from igblast as well as database in folder

workingDir=/home/datier01/bCellDevelopment
NPROC=12

mkdir ${workingDir}/logs
SAMPLES=samples.txt

cat $SAMPLES | while read SAMPLE
do
	S=${SAMPLE}
	DATE=$(date +%Y-%m-%d)
	JOBNAME=${S}_${DATE}
	echo sample ${S} named ${JOBNAME} logged in ${LOGDIR}
	sbatch --time=72:00:00 -c 1 -p compute -o ${workingDir}/logs/jobs_${S}.txt -e ${workingDir}/logs/jobs_${S}.err -n 10 --wrap="sh ${workingDir}/PipelineFixUMI_igblast3_PreTigger_PreNovel_6Sept2023.sh ${S}"
	#qsub -q long -l nodes=1:ppn=${NPROC},walltime=164:00:00 -v S=${S} -N ${JOBNAME} -e ${LOGDIR} -o ${LOGDIR} ${workingDir}/PipelineFixUMI_igblast3_PreTigger_PreNovel_6Sept2023.sh
done
