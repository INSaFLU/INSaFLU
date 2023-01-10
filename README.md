<p align="center"><img src="static/insa/logo_insaflu_new.png" alt="INSaFLU" width="300"></p>

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)


**INSaFLU-TELEVIR platform (https://insaflu.insa.pt) is a free bioinformatics web-based (but also locally installable) suite that deals with primary sequencing data (Illumina, Ion Torrent and Oxford Nanopore Technologies reads) towards:**

   -	**metagenomics virus detection** (from reads to virus detection) 

   -	**routine genomic surveillance** (from reads to mutations detection, consensus generation, virus classification, alignments, “genotype-phenotype” screening, phylogenetics, integrative Nextstrain phylogeographical and temporal analysis etc). 

**INSaFLU-TELEVIR versatility and functionality is expected to supply public health laboratories and researchers with a user-oriented “start-to-end” bioinformatics framework that can potentiate a strengthened and timely detection and monitoring of viral (emerging) threats.**



<p align="center"><img src="static/insa/main_modules.png" alt="INSaFLU" width="800"></p>

- **Online tool:** https://insaflu.insa.pt
- **Documentation / Tutorial:** https://insaflu.readthedocs.io/en/latest/
- **Code:** https://github.com/INSaFLU/INSaFLU
- **Easy local installation:** https://github.com/INSaFLU/docker



## How to cite

If you use INSaFLU in your work, please cite Borges V, Pinheiro M et al. Genome Medicine (2018) 10:46, [https://doi.org/10.1186/s13073-018-0555-0](https://doi.org/10.1186/s13073-018-0555-0)

# Bioinformatics pipeline

## Main contributors
Miguel Pinheiro, João Dourado Santos, Daniel Sobral, Vitor Borges 

# Installation

For an easy and rapid installation using docker see here [https://github.com/INSaFLU/docker](https://github.com/INSaFLU/docker).

This installation is oriented for Ubuntu Server 16.04 and Centos 7.X.
There are several steps and packages to install, so, please, be patient. First, it is necessary to install and configure all bioinformatics software, then the database, batch-queuing system and, finally, the web site.

The user "flu_user" is used in all operations and it is going to be the user to run the apache web server.

## General packages

###Some general packages to install in Ubuntu 18.04: 

```
$ sudo apt install binutils libproj-dev gdal-bin dos2unix parallel
$ sudo apt install postgresql-10
$ sudo apt install postgresql-10-postgis-2.4
$ sudo apt install postgresql-10-postgis-scripts
$ sudo apt install python3
$ sudo apt install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl
```

###Some general packages to install in Centos 7.X: 

```
$ sudo yum install gdal gdal-devel dos2unix parallel
$ sudo yum install postgis-10
$ sudo yum install postgresql-devel
$ sudo yum install python3
$ sudo yum install perl-Time-Piece perl-XML-Simple perl-Digest-MD5 git java perl-CPAN perl-Module-Build
$ sudo cpan -i Bio::Perl
# sudo yum install https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-1.x86_64.rpm

```

:warning: Important, must be blast 2.7.1 version because of Abricate. Newer versions fails on database creation.

## Bioinformatics software

The software can be installed in this directory "/usr/local/software/insaflu". If you choose other directory it is necessary to edit the file "constants/software_names.py" and set the variable "DIR_SOFTWARE".
 
```
$ sudo mkdir -p /usr/local/software/insaflu
$ sudo chown flu_user:flu_user /usr/local/software/insaflu
```

Software to install:

* [IGVTools](https://software.broadinstitute.org/software/igv/igvtools) 2.3.98
* [SPAdes](http://cab.spbu.ru/software/spades/) 3.11.1
* [Abricate](https://github.com/tseemann/abricate) 0.8-dev
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.11.9
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) 0.27
* [Bamtools](https://github.com/pezmaster31/bamtools) 2.5
* [Prokka](https://github.com/tseemann/prokka) 1.12
* [Mauve](http://darlinglab.org/mauve/mauve.html) 2.4.0, Feb 13 2015
* [Mafft](https://mafft.cbrc.jp/alignment/software/) 7.313
* [seqret](http://emboss.sourceforge.net/download/) (EMBOSS) 6.6.0.0
* [FastTreeDbl](http://www.microbesonline.org/fasttree/) 2.1.10 Double precision
* [freebayes](https://github.com/ekg/freebayes) v1.1.0-54-g49413aa - Also need some scripts available in freebays
* [Snippy](https://github.com/tseemann/snippy) 3.2-dev
	* samtools 1.3
	* bgzip 1.3
	* tabix 1.3
	* snpEff 4.1l - Important, it's necessary to use this version. Recent versions have a problem when variants involve more than one base.
	* freebayes v1.1.0-54-g49413aa
		
Some scripts to install:

* [convertAlignment.pl](https://github.com/lskatz/lyve-SET/blob/master/scripts/convertAlignment.pl) 
	* this script need to be installed in <SoftwareNames.DIR_SOFTWARE>/scripts/convertAlignment.pl
* [Fastq-tools](https://github.com/dcjones/fastq-tools) 0.8	

:warning: Important, edit fastqc file `$ vi <install software path>/FastQC/0.11.9/FastQC/fastqc` and change the line `my $memory = 250 * $threads;` to `my $memory = 1000 * $threads;`

:warning: Important, copy the file `bin/snippy-vcf_to_tab` to `bin/snippy-vcf_to_tab_add_freq` and do this change:

```
$ cd /usr/local/software/insaflu/snippy/bin
$ cp snippy-vcf_to_tab snippy-vcf_to_tab_add_freq
$ vi snippy-vcf_to_tab_add_freq

and change the line 57 from:
print join("\t", qw(CHROM POS TYPE REF ALT EVIDENCE), @ANNO), "\n";
to
print join("\t", qw(CHROM POS TYPE REF ALT FREQ), @ANNO), "\n";
```

:warning: Important, change shebang in spades.py file from `#!/usr/bin/env python` to `#!/usr/bin/env python3`. 

:warning: Important, change snippy script to allow snpEff 4.1 version

	#xpto@brazil:/usr/local/software/insaflu/snippy/bin$ diff snippy snippy~
	90c90
	< parse_version( 'snpEff -version',     4.1, qr/(\d+\.\d+)/           );
	---
	> parse_version( 'snpEff -version',     4.3, qr/(\d+\.\d+)/           );
	
## Database PostgreSQL

	* postgresql 10.X
		* create a database and a user. Then reflect these names in ".env" file in root path of web site.

## Sun Grid Engine/Open Grid Engine

	Software:
	* gzip
	* [Sun Grid Engine/Open Grid Engine](https://arc.liv.ac.uk/downloads/SGE/releases)
		* [download 8.1.9 version](https://arc.liv.ac.uk/downloads/SGE/releases/8.1.9/sge_8.1.9.tar.xz)
		* queues that will be created:
			* all.q - generic queue
			* fast.q - to run quick process
			* queue_1.q and queue_2.q - to run slow process


Install SGE/OGE tips
 
```
$ sudo groupadd -g 58 gridware
$ sudo useradd -u 63 -g 58 -d /opt/sge sgeadmin
$ cd ~
$ mkdir sge; cd sge
$ wget https://arc.liv.ac.uk/downloads/SGE/releases/8.1.9/sge_8.1.9.tar.xz
$ tar -xJvf sge_8.1.9.tar.xz
$ cd sge-8.1.9/source
$ scripts/bootstrap.sh

### centos version
$ sudo yum install hwloc-devel openssl-devel pam-devel libXt-devel motif motif-devel readline-devel
### ubuntu
$ sudo apt-get install libhwloc-dev libssl-dev

$ ./aimk -no-java -no-jni

### caveat in last command if something like this `ed.screen.c:(.text+0x247c): undefined reference to `tputs'` appears in the screen
$ cd 3rdparty/qtcsh/LINUXAMD64
    Add  "-lreadline -lncurses" to the end of command that fails
$ cd ../../..
$ ./aimk -no-java -no-jni
### end caveat

$ sudo su
# export SGE_ROOT=/opt/sge
# scripts/distinst -local -allall -noexit
# chown -R sgeadmin:gridware /opt/sge
# cd $SGE_ROOT
# ./install_qmaster
# . /opt/sge/default/common/settings.sh
# ./install_execd

### create a file to set the environment variables to SGE
$ sudo vi /etc/profile.d/sun-grid-engine.sh
## add the follow line to the file
 . /opt/sge/default/common/settings.sh
```


Configure queues

Go to the folder `example_script_sge_add_queue` and change second line in the files `"grid_add_all_hosts.txt"` and replace 'brazil' to your computer name. Mine is 'brazil'

Get your computer name:

```
$ uname -a
Linux brazil 4.15.0-50-generic #54-Ubuntu SMP Mon May 6 18:46:08 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux
```


Add your name to manage list, to obtain permissions to change SGE configurations:

```
$ sudo qconf -am <your name>
```

:warning: If you get an error about `qconf not found` or `SGE_ROOT not set`, do something like this:

```
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
```


Then run:

```
$ cd example_script_sge_add_queue
$ qconf -Ahgrp grid_add_all_hosts.txt
```

Show all groups:

```
$ qconf -shgrpl
$ qconf -shgrp_resolved @allhosts
brazil
```

If you want to delete a group name:

```
$ qconf -dhgrp <a group name>
```

To add the queues:

```
$ qconf -Aq grid_add_queue_all.txt
$ qconf -Aq grid_add_queue_fast.txt
$ qconf -Aq grid_add_queue_queue_1.txt
$ qconf -Aq grid_add_queue_queue_2.txt
```

If you want to delete a queue:

```
$ qconf -dq <a queue name>
```

Show all info

```
$ qconf -sq <queue name>
```

Edit queues. If you want to change slots change the number in 'slots'.

```
$ qconf -mq <queue name> 
```

Change the default `schedule_interval` from `0:0:15` to `0:0:5`. This setting specifies how often the scheduler checks for new jobs.

```
$ qconf -msconf
```

After the OGE/SGE configuration you need to have these queue names in your system.

```
$ qstat -f
queuename                      qtype resv/used/tot. load_avg arch          states
---------------------------------------------------------------------------------
all.q@brazil               BIP   0/0/2          1.19     lx26-amd64    
---------------------------------------------------------------------------------
fast.q@brazil              BIP   0/0/1          1.19     lx26-amd64        
---------------------------------------------------------------------------------
queue_1.q@brazil           BIP   0/0/1          1.19     lx26-amd64    
---------------------------------------------------------------------------------
queue_2.q@brazil           BIP   0/0/1          1.19     lx26-amd64 
```

:warning: `brazil` is the name of the computer where the installation is. You have other certainly. The computer name need to be in `/etc/hosts` with the IP address and not with `localhost` to SGE work properly.
Example:

```
$ cat /etc/hosts
127.0.0.1	localhost
::1     ip6-localhost ip6-loopback

192.168.1.14	brazil
```

Of course you have a different IP address from '192.168.1.14'

More [help to configure queues](https://peteris.rocks/blog/sun-grid-engine-installation-on-ubuntu-server/). 

## INSaFLU website


```

$ sudo mkdir -p /usr/local/web_site
$ sudo mkdir -p /var/log/insaFlu
$ sudo chown flu_user:flu_user /usr/local/web_site
$ sudo chown flu_user:flu_user /var/log/insaFlu
$ cd /usr/local/web_site
$ git clone https://github.com/INSaFLU/INSaFLU.git
$ cd INSaFLU
$ sudo pip3 install -r requirements.txt
$ cp .env_model .env
```

Edit the file ".env" and config all variables. Define also a backend to the email. I have defined a posix server. 

To create the database

```
$ python3 manage.py migrate
```

To create a super user, it is going to be the administrator user account

```
$ python3 manage.py createsuperuser
```

To join all files, in "static_all" path, that is necessary to run the web site and then read default databases. All data that belong to databases are in "<insaflu path>/static_all/db/..."

```
$ python3 manage.py collectstatic
$ python3 manage.py load_default_files

```



Test if all bioinformatic tolls are installed  

```
$ cd /usr/local/web_site
$ python3 manage.py test constants.tests_software_names

```

Test everything

```
$ cd /usr/local/web_site
$ python3 manage.py test

```

:warning: All tests must pass otherwise something is not working properly.

If all tests passed you can test immediately it is working:

```
$ cd /usr/local/web_site
$ python3 manage.py runserver
```

Go to your internet explorer and write the ip of the computer where the web site is installed "<ip server>:8000". If it is in same computer can be "localhost:8000".
If it is working let's go to install in a Apache web server. If you prefer, can be in a Nginx web server too.
 

## Apache web server

###Config apache2 in Centos 7.X:


Add `flu_user` to the `apache` group and add `insaflu.conf` to apache2.

```
$ sudo usermod -a -G flu_user apache
$ sudo vi /etc/httpd/conf.d/insaflu.conf

<VirtualHost *:80>

	# General setup for the virtual host, inherited from global configuration

	ServerName insaflu.pt

        Alias /media /usr/local/web_site/INSaFLU/media
        Alias /static /usr/local/web_site/INSaFLU/static_all
        <Directory "/usr/local/web_site/INSaFLU/static_all">
                Require all granted
        </Directory>
        <Directory "/usr/local/web_site/media">
                Options FollowSymLinks
                AllowOverride None
                Require all granted
        </Directory>

        #### for log files
        <Directory "/var/log/insaFlu">
                Require all granted
        </Directory>

        <Directory "/usr/local/web_site/INSaFLU/fluwebvirus">
            <Files "wsgi.py">
                Require all granted
            </Files>
        </Directory>
	
	WSGIDaemonProcess flu_user.insa.pt user=flu_user group=flu_user python-path=/usr/local/web_site/INSaFLU/fluwebvirus;/usr/lib/python3.<minor version of your python>/site-packages
        WSGIProcessGroup flu_user.insa.pt
        WSGIScriptAlias / /usr/local/web_site/INSaFLU/fluwebvirus/wsgi.py

# Use separate log files for the SSL virtual host; note that LogLevel
# is not inherited from httpd.conf.
ErrorLog /var/log/httpd/insaflu_error.log
TransferLog /var/log/httpd/insaflu_transfer.log
LogLevel warn

</VirtualHost> 

$ sudo yum install httpd-devel mod_wsgi
$ sudo updatedb

## START small caveat...
$ mv /etc/httpd/modules/mod_wsgi.so /etc/httpd/

OR
$ locate mod_wsgi.so 
$ sudo ln -s <last hit for the locate> /etc/httpd/modules/mod_wsgi.so
## END small caveat...

$ sudo systemctl restart httpd
$ sudo systemctl status httpd
```

###Config apache2 in Ubuntu 16.X:

Add `flu_user` to the `apache` group and add `insaflu.conf` to apache2.

```
$ sudo usermod -a -G flu_user apache
$ sudo apt install libapache2-mod-wsgi-py3
$ sudo vi /etc/apache2/sites-available/insaflu.conf

<VirtualHost *:80>

	# General setup for the virtual host, inherited from global configuration

	ServerName insaflu.pt

        Alias /media /usr/local/web_site/media
        Alias /static /usr/local/web_site/static_all
        <Directory "/usr/local/web_site/static_all">
                Require all granted
        </Directory>
        <Directory "/usr/local/web_site/media">
                Options FollowSymLinks
                AllowOverride None
                Require all granted
        </Directory>

        #### for log files
        <Directory "/var/log/insaFlu">
                Require all granted
        </Directory>

        <Directory "/usr/local/web_site/insaflu">
            <Files "wsgi.py">
                Require all granted
            </Files>
        </Directory>
	
	WSGIDaemonProcess flu_user.insa.pt user=flu_user group=flu_user python-path=/usr/local/web_site/insaflu;/usr/lib/python3.<minor version of your python>/site-packages
	WSGIProcessGroup flu_user.insa.pt
	WSGIScriptAlias / /usr/local/web_site/insaflu/wsgi.py

# Use separate log files for the SSL virtual host; note that LogLevel
# is not inherited from httpd.conf.
ErrorLog /var/log/apache2/insaflu_error.log
TransferLog /var/log/apache2/insaflu_transfer.log
LogLevel warn


</VirtualHost> 

$ sudo a2ensite insaflu.conf
$ sudo systemctl restart apache2
$ sudo systemctl status apache2
```


:warning:  Add "AddType application/octet-stream .bam" to httpd.conf in "IfModule mime_module" element, for the IGV viewer.

## Create users without access to INSaFLU web page

Go to your internet explorer and put this address `http://127.0.0.1:80/admin/`
Make the authentication with your superuser credentials and in `AUTHENTICATION AND AUTHORIZATION` you can create new accounts. 

## Remove files from file system removed by the user on web site

You can remove the original fastq.gz files from system because they are not used anymore. The Trimmomatic result fastq files are the ones that are going to be used. 
You can can also remove files that belong to the samples, references, uploaded in batch and project samples that were deleted in web site by the users.
This operation will save several GB in your hard drives. 

:warning:By default, only files with 10 days after been removed in web site will be removed in file system. 
:warning:The original fastq.gz files will be removed after 10 days of being processed by Trimmomatic. 

To identify the files that can be removed:

```
$ cd <where your INSaFLU is installed>
$ python3 manage.py run_remove_files --only_identify_files true
```
A log file will be created with this information in `/var/log/insaflu/remove_files.log`

To remove the files permanently from file system:
:warning: The files can't be recovered.

```
$ cd <where your INSaFLU is installed>
$ python3 manage.py run_remove_files --only_identify_files false
```

Tip:

You can create a cron job to run this task every week.
