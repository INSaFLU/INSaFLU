#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -q all.q
#$ -N job_name_job_name

echo $HOSTNAME > /tmp/sge_${JOB_ID}.out
date >> /tmp/sge_${JOB_ID}.out
echo "start waiting" >> /tmp/sge_${JOB_ID}.out
sleep 2m
date >> /tmp/sge_${JOB_ID}.out
echo "end" >> /tmp/sge_${JOB_ID}.out
