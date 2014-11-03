#!/bin/bash

#SBATCH --partition=normal
#SBATCH --job-name=bagofwaves
#SBATCH --output=calculation.%j.log
#SBATCH --time=2-00:00:00
#SBATCH --time-min=12:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmcgibbo@gmail.com
#SBATCH --ntasks-per-node=8

cd $SLURM_SUBMIT_DIR
PDIR=$HOME/projects/bagofwaves

for i in `seq $SLURM_CPUS_ON_NODE`; do
    mongo-task-submit -e $PDIR/.env -t $PDIR/task.yaml --loop \
         > mongo-task-submit.$SLURM_JOB_ID.$i.stdout \
	2> mongo-task-submit.$SLURM_JOB_ID.$i.stderr &
done

wait
