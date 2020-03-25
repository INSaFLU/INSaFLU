#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -q all.q
#$ -N job_name_wait
#$ -hold_jid job_name_job_name

echo $HOSTNAME > /tmp/sge_wait_${JOB_ID}_${JOB_NAME}.out
date >> /tmp/sge_wait_${JOB_ID}_${JOB_NAME}.out
echo "start waiting" >> /tmp/sge_wait_${JOB_ID}_${JOB_NAME}.out
sleep 2m
date >> /tmp/sge_wait_${JOB_ID}_${JOB_NAME}.out
echo "end" >> /tmp/sge_wait_${JOB_ID}_${JOB_NAME}.out
