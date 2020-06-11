#!/bin/sh

#### SETUP ####

# Specify path to base directory
BASEDIR='/scratch/ias41/ae_code/psm_aeolus'

# Specify file with MySQL and FAERS_AEOLUS connection details
MYSQL_DEFAULTS=$BASEDIR/mysql_info.txt

# Number of CPU cores to use
NUMBER_OF_CORES=4

# Specify file with \n separated AEOLUS drug concept ids (used to restrict analysis, e.g. to drugs of interest/with other data). 
# Otherwise will run for all drugs in AEOLUS database.
#AEOLUS_IDS_LIST=$BASEDIR/data/aeolus_ids_with_bioact_for_ae_detection.txt
#AEOLUS_IDS_LIST=$BASEDIR/data/aeolus_ids_wo_bioact.txt
#AEOLUS_IDS_LIST=$BASEDIR/data/example_aeolus_ids.txt

#### END OF SETUP ####

mkdir -p $BASEDIR/results

echo "Starting SQL queries for correlations..."

mysql --defaults-file=$MYSQL_DEFAULTS < $BASEDIR/src/isr_drugs.sql > $BASEDIR/data/isr_drugs.txt
mysql --defaults-file=$MYSQL_DEFAULTS < $BASEDIR/src/primaryid_drugs.sql > $BASEDIR/data/primaryid_drugs.txt
mysql --defaults-file=$MYSQL_DEFAULTS < $BASEDIR/src/isr_indications.sql > $BASEDIR/data/isr_indications.txt
mysql --defaults-file=$MYSQL_DEFAULTS < $BASEDIR/src/primaryid_indications.sql > $BASEDIR/data/primaryid_indications.txt
mysql --defaults-file=$MYSQL_DEFAULTS < $BASEDIR/src/drug_names.sql > $BASEDIR/data/drug_names.txt
mysql --defaults-file=$MYSQL_DEFAULTS < $BASEDIR/src/indication_names.sql > $BASEDIR/data/indication_names.txt

echo "Finished SQL queries.."
echo "Starting to insert results in sqlite database.. "

python ./create_sqlite_db.py --basedir $BASEDIR
echo "Finished inserting into the database"

echo "Starting to analyse drug-adverse event associations, using $NUMBER_OF_CORES cores"

if [ -z ${AEOLUS_IDS_LIST+x} ]; then 
   python ./faers_integrated.py --number_of_cores $NUMBER_OF_CORES --basedir $BASEDIR --mysql_defaults $MYSQL_DEFAULTS --min_nr_cases 5 
else 
   python ./faers_integrated.py --number_of_cores $NUMBER_OF_CORES --basedir $BASEDIR --mysql_defaults $MYSQL_DEFAULTS --min_nr_cases 5 --aeolus_ids_list $AEOLUS_IDS_LIST
fi

echo "Finished"
