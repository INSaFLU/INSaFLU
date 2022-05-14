

Change second line in the files "grid_add_host_fast.txt", "grid_add_host_queue_1.txt" and "grid_add_host_queue_2.txt" and replace 'brazil' to your computer name. Mine is 'brazil'

Get your computer name:
$ uname -a
Linux brazil 4.15.0-50-generic #54-Ubuntu SMP Mon May 6 18:46:08 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux


Add your name to manage list
$ sudo qconf -am <your name>

If get error about "qconf" no found or "SGE_ROOT" not set, do this:

### this need to be improved
$ sudo locate settings.sh
$ sudo chmod a+x /opt/sge/default/common/settings.sh
$ /opt/sge/default/common/settings.sh
$ env | grep sge
$ sudo su
# PATH=   ##### /opt/sge/bin:/opt/sge/bin/lx-amd64:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin ### the path from last env
# SGE_ROOT=  ### /opt/sge   ### the path from last env
# export SGE_ROOT;
# qconf -am <your name>


Then run:

$ qconf -Ahgrp grid_add_all_hosts.txt

### show all groups
$ qconf -shgrpl

### if you want to delete a group namem
$ qconf -dhgrp <a group name>

$ qconf -Aq grid_add_queue_all.txt
$ qconf -Aq grid_add_queue_fast.txt
$ qconf -Aq grid_add_queue_queue_1.txt
$ qconf -Aq grid_add_queue_queue_2.txt

### if you want remove a queue
$ qconf -dq <a queue name>

### show all queues
$ qstat -f

### show all info
$ qconf -sq <queue name>

### edit queues
$ qconf -mq <queue name>  


### change the default schedule_interval from 0:0:15 to 0:0:5. This setting specifies how often the scheduler checks for new jobs.
$ qconf -msconf

### Cahnge priority to the job 38127
$ qalter -p -1020 38127


#### show all scripts 
$ ll /opt/sge/default/spool/qmaster/job_scripts/

