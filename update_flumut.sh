#!/bin/bash
## Update flumut
## 
## $ crontab -e
## 0 1 * * * /usr/local/web_site/INSaFLU/update_flumut.sh

set -e 
## set the directory where Insafluy Webserver is
cd /insaflu_web/INSaFLU
python3 manage.py update_flumut
