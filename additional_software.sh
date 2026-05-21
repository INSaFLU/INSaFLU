sudo apt install libpq-dev
sudo apt install build-essential python3 python3.10-dev python3.10-venv python3.10-distutils
sudo apt install binutils libproj-dev gdal-bin dos2unix parallel

# Install PostgreSQL 16
sudo apt install dirmngr ca-certificates software-properties-common apt-transport-https lsb-release curl -y
curl -fSsL https://www.postgresql.org/media/keys/ACCC4CF8.asc | gpg --dearmor | sudo tee /usr/share/keyrings/postgresql.gpg > /dev/null
echo deb [arch=amd64,arm64,ppc64el signed-by=/usr/share/keyrings/postgresql.gpg] http://apt.postgresql.org/pub/repos/apt/ $(lsb_release -cs)-pgdg main | sudo tee /etc/apt/sources.list.d/postgresql.list

echo deb [arch=amd64,arm64,ppc64el signed-by=/usr/share/keyrings/postgresql.gpg] http://apt.postgresql.org/pub/repos/apt/ $(lsb_release -cs)-pgdg-testing main | sudo tee /etc/apt/sources.list.d/postgresql.list
sudo apt update
sudo apt install postgresql-client-16 postgresql-16
###
apt-get update \
&& apt-cache showpkg postgresql-$PG_MAJOR-postgis-$POSTGIS_MAJOR \
&& apt-get install -y --no-install-recommends \
# ca-certificates: for accessing remote raster files;
#   fix: https://github.com/postgis/docker-postgis/issues/307
ca-certificates \
\
postgresql-$PG_MAJOR-postgis-$POSTGIS_MAJOR=$POSTGIS_VERSION \
postgresql-$PG_MAJOR-postgis-$POSTGIS_MAJOR-scripts \
&& rm -rf /var/lib/apt/lists/*


### Install additional packages
sudo apt install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl
