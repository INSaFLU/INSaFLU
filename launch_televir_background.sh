#!/bin/bash

SAMPLE_ID=$1

/usr/bin/python3 manage.py !/bin/bash

PROJECT_ID=$1

/usr/bin/python3 manage.py run_televir_background \
--project_id $PROJECT_ID \
--out_dir /tmp/insaFlu/test_bg \
--log_dir ../tmp_log \
--max_threads 4

