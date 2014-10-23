#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m abe
#$ -M rmcgibbo+jobs@gmail.com
#$ -l h_rt=1:00:00,h_vmem=4.1G
#$ -t 1-4
#$ -notify

PDIR=$HOME/projects/bagofwaves

mongo-task-submit -e $PDIR/.env -t $PDIR/task.yaml --loop \
    > mongo-task-submit.$SGE_TASK_ID.stdout \
    2> mongo-task-submit.$SGE_TASK_ID.stderr