#!/bin/bash
FS="/"
THRESHOLD=90
OUTPUT=($(LC_ALL=C df -P ${FS}))
CURRENT=$(echo ${OUTPUT[11]} | sed 's/%//')
[ $CURRENT -gt $THRESHOLD ] && echo "Warning: INSAflu file system usage $CURRENT%" | mail -s "file system: $CURRENT%" insaflu@insa.min-saude.pt
[ $CURRENT -gt 50 ] && echo "Warning: INSAflu file system usage $CURRENT%" | mail -s "file system: $CURRENT%" monsantopinheiro@gmail.com
