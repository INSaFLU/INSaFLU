#!/bin/bash
## Upadte pangolin
##
## $ crontab -e
## 0 1 * * * /usr/local/web_site/INSaFLU/update_pangolin.sh

set -e
## set the directory where Insaflu WebServer is
cd /insaflu_web/INSaFLU
echo `pwd`
python3 manage.py update_pangolin
