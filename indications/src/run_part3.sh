#!/bin/sh

BASEDIR='/scratch/ias41/ae_code'

PROJECT_DIR=$BASEDIR/indications
MAPPED_COMPOUNDS_DB=$BASEDIR/compound_mapping/results/201903_mapped_compounds_calculon.db
MESH2MEDDRA_DF=$PROJECT_DIR/data/mesh2meddra_df.txt
MEDDRA_HIERARCHY_OUTPUT=$PROJECT_DIR/data/meddra_terms_indications_hierarchy_output_22_1.xlsx
ALL_AES_MEDDRA_HIER_FILE=$PROJECT_DIR/data/all_aeolus_aes_hierarchy_output_22_1.xlsx
PREVENT_DICT=$PROJECT_DIR/data/rxnorm2may_prevent.pkl
TREAT_DICT=$PROJECT_DIR/data/rxnorm2may_treat.pkl

jupyter nbconvert --to script molregno2indications_via_hlt.ipynb 
python ./molregno2indications_via_hlt.py --project_dir $PROJECT_DIR --mapped_compounds_db $MAPPED_COMPOUNDS_DB --meddra_hier_output $MEDDRA_HIERARCHY_OUTPUT --mesh2meddra_df $MESH2MEDDRA_DF --all_aes_meddra_hier_file $ALL_AES_MEDDRA_HIER_FILE --may_prevent_dict $PREVENT_DICT --may_treat_dict $TREAT_DICT

