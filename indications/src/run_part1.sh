#!/bin/sh

BASEDIR='/scratch/ias41/ae_code'

PROJECT_DIR=$BASEDIR/indications
MYSQL_DEFAULTS=$BASEDIR/psm_aeolus/mysql_info.txt

jupyter nbconvert --to script get_aeolus_aes_for_hier_analysis.ipynb
python ./get_aeolus_aes_for_hier_analysis.py --mysql_defaults $MYSQL_DEFAULTS  --project_dir $PROJECT_DIR
