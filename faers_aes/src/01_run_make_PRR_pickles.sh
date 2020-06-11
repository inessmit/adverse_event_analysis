#!/bin/sh

BASEDIR='/scratch/ias41/ae_code'

PROJECT_DIR=$BASEDIR/faers_aes
MAPPED_COMPOUNDS_DB=$BASEDIR/compound_mapping/results/201903_mapped_compounds_calculon.db

PSM_RESULTS_DIR=$BASEDIR/psm_aeolus/results/data

jupyter nbconvert --to script make_PRR_pickles.ipynb 
python ./make_PRR_pickles.py --project_dir $PROJECT_DIR --mapped_compounds_db $MAPPED_COMPOUNDS_DB --psm_results_dir $PSM_RESULTS_DIR --background_condition 'all_random_controls'

