#!/bin/sh

BASEDIR='/scratch/ias41/ae_code'
PROJECT_DIR=$BASEDIR/indications
MAPPED_COMPOUNDS_DB=$BASEDIR/compound_mapping/results/201903_mapped_compounds_calculon.db
UMLS_CONSO_FILE=$PROJECT_DIR/data/MRCONSO_2019AA.RRF

jupyter nbconvert --to script rxnorm2indications.ipynb
python ./rxnorm2indications.py --project_dir $PROJECT_DIR --mapped_compounds_db $MAPPED_COMPOUNDS_DB --mrconso_file $UMLS_CONSO_FILE
