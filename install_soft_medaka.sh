#!/bin/bash
set -e

cd /usr/local/software/insaflu
pip3 install virtualenv
virtualenv medaka --python=python3 --prompt "(medaka 1.2.1) "
. medaka/bin/activate
pip3 install medaka==1.2.1
cd medaka

#Install minimap
# cd 
mkdir -p extra_software
cd extra_software
git clone https://github.com/lh3/minimap2.git
cd minimap2/
make

#Install HTSLIB
cd ..
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar -vxjf htslib-1.9.tar.bz2
cd htslib-1.9
make
cd ..
rm htslib-1.9.tar.bz2

#Install SAMTOOLS
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -vxjf samtools-1.9.tar.bz2
cd samtools-1.9
make
cd ..
rm samtools-1.9.tar.bz2

#Install BCFTools
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar -vxjf bcftools-1.9.tar.bz2
cd bcftools-1.9
make
cd ..
rm bcftools-1.9.tar.bz2

### links
cd ../bin
ln -s ../extra_software/minimap2/minimap2 .
ln -s ../extra_software/htslib-1.9/bgzip .
ln -s ../extra_software/htslib-1.9/tabix .
ln -s ../extra_software/samtools-1.9/samtools .
ln -s ../extra_software/bcftools-1.9/bcftools .

