#!/bin/bash

PROJECT_ID=$1

/usr/bin/python3 manage.py run_televir_background \
--project_id $PROJECT_ID \
--out_dir /tmp/insaFlu/test_bg \
--log_dir ../tmp_log \
--metagenomics \
--max_threads 6

