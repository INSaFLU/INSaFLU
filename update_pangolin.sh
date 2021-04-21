## Upadte pangolin
set -e
## set the directory where Insaflu WebServer is
cd /insaflu_web/INSaFLU
echo `pwd`
python3 manage.py update_pangolin
