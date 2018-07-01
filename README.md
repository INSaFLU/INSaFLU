<p align="center"><img src="static/insa/logo_insaflu_new.png" alt="INSaFLU" width="300"></p>


[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)


# INSaFLU
INSaFLU (“INSide the FLU”) is an influenza-oriented bioinformatics free web-based platform for an effective and timely whole-genome-sequencing-based influenza laboratory surveillance.

INSaFLU is freely available at [https://insaflu.insa.pt](https://insaflu.insa.pt)
Documentation (latest) for each INSaFLU module is provided at [http://insaflu.readthedocs.io/](http://insaflu.readthedocs.io/)

## Synopsis

INSaFLU (“INSide the FLU”) is a bioinformatics free web-based suite that deals with primary NGS data (reads) towards the automatic generation of the output data that are actually the core first-line “genetic requests” for effective and timely influenza laboratory surveillance (e.g., type and sub-type, gene and whole-genome consensus sequences, variants annotation, alignments and phylogenetic trees).

## Main features

Highlights / Main advantages
* open to all, free of charge, user-restricted accounts
* applicable to NGS data collected from any amplicon-based schema
* allows advanced, multi-step software intensive analyses in a user-friendly manner without previous advanced training in bioinformatics
* allows integrating data in a cumulative manner, thus fitting the analytical dynamics underlying the continuous epidemiological surveillance during flu epidemics
* outputs are provided in nomenclature-stable and standardized formats and can be explored in situ or through multiple compatible downstream applications for data analysis and visualization

Main outputs
INSaFLU yields:
* influenza type and subtype/lineage 
* gene and whole-genome consensus sequences
* annotation of variants and intra-host minor variants
* gene, protein and genome alignments
* gene- and genome-scale phylogenetic trees

Other features:
INSaFLU also automatically provides:
* raw NGS data quality analysis and improvement
* a rapid snapshot of whole-genome backbone of each virus (draft assembled contigs are assigned to each viral segment and to close related reference influenza viruses). 
* coverage statistics
* detection of putative mixed infections

## How to cite

If you use INSaFLU in your work, please cite Borges V, Pinheiro M et al. Genome Medicine (2018) 10:46, [https://doi.org/10.1186/s13073-018-0555-0](https://doi.org/10.1186/s13073-018-0555-0)

# Bioinformatics pipeline

## Authors
Miguel Pinheiro
Vitor Borges 

# Installation

This installation is oriented for Ubuntu Server 16.04 and Centos 7.X.
There are several steps and packages to install, so, please, be patient. First, it is necessary to install and configure all bioinformatics software, then the database, batch-queuing system and, finally, the web site.

The user "flu_user" is used in all operations and it is going to be the user to run the apache web server.

## General packages

###Some general packages to install in Ubuntu 16.X: 

```
$ sudo apt install binutils libproj-dev gdal-bin
$ sudo apt install postgis*
$ sudo apt install bioperl
$ sudo apt install python3
$ sudo apt install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl
```

###Some general packages to install in Centos 7.X: 

```
$ sudo yum install gdal gdal-devel 
$ sudo yum install postgis
$ sudo yum install python3
$ sudo yum install perl-Time-Piece perl-XML-Simple perl-Digest-MD5 git java perl-CPAN perl-Module-Build
$ sudo cpan -i Bio::Perl
```

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
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.11.5
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) 0.27
* [Bamtools](https://github.com/pezmaster31/bamtools) 2.5
* [Prokka](https://github.com/tseemann/prokka) 1.2
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

:warning: Important, change snippy script to allow snpEff 4.1 version

	#xpto@brazil:/usr/local/software/insaflu/snippy/bin$ diff snippy snippy~
	90c90
	< parse_version( 'snpEff -version',     4.1, qr/(\d+\.\d+)/           );
	---
	> parse_version( 'snpEff -version',     4.3, qr/(\d+\.\d+)/           );
	
## Database PostgreSQL

	* postgresql 9.X
		* create a database and a user. Then reflect these names in ".env" file in root path of web site.

## Operation System

	Software necessary:
	* gzip
	* [Sun grid engine/Open Grid Engine](http://gridscheduler.sourceforge.net/)
		* [download 6.2 version](https://sourceforge.net/projects/gridscheduler/files/SGE6.2u5p2/)
		* create the queues:
			* all.q - generic queue
			* fast.q - run quick process
			* queue_1.q and queue_2.q - run slow process


Install SGE tips
 
```
$ SGE_ROOT=/home/software/sge
$ ./aimk -no-java -no-jni
$ scripts/distinst -all -local -noexit
% cd $SGE_ROOT
% ./install_qmaster
% ./install_execd
```

Configure queues with this [help.](https://peteris.rocks/blog/sun-grid-engine-installation-on-ubuntu-server/)

After the SGE configuration you need to have these queue names in your system.

```
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

## INSaFLU website


```

$ sudo mkdir -p /usr/local/web_site
$ sudo mkdir -p /var/log/insaflu
$ sudo chown flu_user:flu_user /usr/local/web_site
$ sudo chown flu_user:flu_user /var/log/insaflu
$ cd /usr/local/web_site
$ git clone https://github.com/INSaFLU/INSaFLU.git
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

To join all files, in "static_all" path, that is necessary to run the web site

```
$ python3 manage.py collectstatic

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
## From IUS repo
$ sudo yum install python3<minor version of your python>u-mod_wsgi
$ sudo vi /etc/httpd/conf.d/insaflu.conf

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

$ sudo a2ensite insaflu
$ sudo systemctl restart apache2
$ sudo systemctl status apache2
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

## Create users without access to INSaFLU web page

Go to your internet explorer and put this address `http://127.0.0.1:8000/admin/`
Make the authentication with your superuser credentials and in `AUTHENTICATION AND AUTHORIZATION` you can create new accounts. 


