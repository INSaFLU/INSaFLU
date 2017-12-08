
cd static_all
rm -r admin*

cd ..
python3 manage.py collectstatic

