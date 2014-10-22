#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -t 1-16
#$ -l h_rt=48:00:00,h_vmem=1.00G

PDIR=$HOME/projects/bagofwaves

mongo-task-submit -e $PDIR/.env -t $PDIR/task.yaml --loop \
    > mongo-task-submit.$SGE_TASK_ID.stdout \
    2> mongo-task-submit.$SGE_TASK_ID.stderr
