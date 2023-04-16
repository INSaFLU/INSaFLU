#! /bin/bash

cd -- "$(dirname -- "$BASH_SOURCE")"

cat *.fastq.gz > concat.fastq.gz
