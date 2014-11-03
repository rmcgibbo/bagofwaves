#!/bin/bash

#PBS -S /bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -e /dev/null
#PBS -o /dev/null
#PBS -V
#PBS -q MP

cd $PBS_O_WORKDIR
PDIR=$HOME/projects/bagofwaves
NO_OF_CORES=`cat $PBS_NODEFILE | egrep -v '^#'\|'^$' | wc -l | awk '{print $1}'`

for i in `seq $NO_OF_CORES`; do
    mongo-task-submit -e $PDIR/.env -t $PDIR/task.yaml --loop \
         > mongo-task-submit.$PBS_JOB_ID.$i.stdout \
	2> mongo-task-submit.$PBS_JOB_ID.$i.stderr &
done

wait