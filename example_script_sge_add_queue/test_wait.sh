for i in `seq 0 13`; do echo $i; qsub example_job.sh; sleep 5s; done;
qsub example_job_wait.sh

