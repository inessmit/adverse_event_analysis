import psm_functions
import logging
import os
import argparse
import multiprocessing
import multiprocessing_logging
import pymysql

# Command line argparse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--aeolus_ids_list", default=None, type=str, help="filepath to file with \n-separated aeolus drug concept IDs for drugs to run AE detection on.")
parser.add_argument("--number_of_cores", type=int, help="number of cores to use for multiprocessing")
parser.add_argument("--basedir", type=str, help="path to basedir")
parser.add_argument("--mysql_defaults", type=str, help="file to mysql details")
parser.add_argument("--min_nr_cases", type=int, default=5, help='minimum number of reports with adverse event')
parser.add_argument("--do_plots", action="store_true", help="Provide this flag to create figures of PS distributions and difference plots")

args = parser.parse_args()

# File locations and directories 
basedir = args.basedir

# Set up logging
if not os.path.exists(basedir + '/logs'):
    os.mkdir(basedir + '/logs')
logging.basicConfig(filename=basedir + '/logs/faers_integrated.log', filemode='a', level=logging.INFO, format='%(asctime)s %(message)s')

# Use multiprocessing_logging module to collect logs from multiprocessing and write to same log file
multiprocessing_logging.install_mp_handler()

logging.info('script started')

indication_names_file = basedir + '/data/indication_names.txt'
drug_names_file = basedir + '/data/drug_names.txt'
primaryid_drugs = basedir + '/data/primaryid_drugs.txt'
primaryid_indications = basedir + '/data/primaryid_indications.txt'
isr_drugs = basedir + '/data/isr_drugs.txt'
isr_indications = basedir + '/data/isr_indications.txt'

correlations_drugs = basedir + '/data/correlations.db'
correlations_indications = basedir + '/data/correlations.db'

# Make directory for figures
if not os.path.exists(basedir + '/results/figures'):
    os.mkdir(basedir + '/results/figures')
if not os.path.exists(basedir + '/results/data'):
    os.makedirs(basedir + '/results/data')
data_path = basedir + '/results/data/'
savefig_path = basedir + '/results/figures/'

# loading dicts, this loads the info for all drugs
drug_names, indication_names, drug2reports, indication2reports, report2attributes_primaryid, report2attributes_isr = psm_functions.load_all_dictionaries(drug_names_file = drug_names_file, indication_names_file = indication_names_file
                                            , primaryid_drugs = primaryid_drugs, primaryid_indications = primaryid_indications
                                            , isr_drugs = isr_drugs, isr_indications = isr_indications)

logging.info('loaded various dicts')

# Open file with RxNorm concept IDs if specified, otherwise use all drugs
if args.aeolus_ids_list != None:
    with open(args.aeolus_ids_list, 'r') as f:
        aeolus_ids_result = [int(i) for i in f.read().splitlines()]
elif args.aeolus_ids_list == None:
    aeolus_conn = pymysql.connect(read_default_file=args.mysql_defaults, unix_socket='/var/run/mysqld/mysqld.sock')
    cursor = aeolus_conn.cursor()
    aeolus_ids_query = """select distinct con.concept_id 
from FAERS_AEOLUS.concept con
join FAERS_AEOLUS.standard_case_drug drug on drug.standard_concept_id = con.concept_id
where con.vocabulary_id = 'RxNorm'
and con.concept_class_id = 'Ingredient'
and con.standard_concept = 'S';"""
    cursor.execute(aeolus_ids_query)
    aeolus_ids_result = [i[0] for i in cursor.fetchall()]
    with open(basedir + '/data/aeolus_ids_for_ae_detection.txt', 'w') as f:
        f.write('\n'.join(str(i) for i in aeolus_ids_result))
    aeolus_conn.close()

# 
all_reports_isr = set(report2attributes_isr.keys())
all_reports_primaryid = set(report2attributes_primaryid.keys())

# Prepare dictionary with reports per AE looks like: all_ae_dict = {isrs: [], primaryids: []}
report2ae_isrs, report2ae_primaryids, ae_concept_id2pt = psm_functions.make_report2ae_dicts(db_defaults=args.mysql_defaults)
logging.info('Loaded AE dictionary')

# Define function for doing all adverse event (AE) detection tasks
def do_all_ae_tasks(drug):
    logging.info('{} starting drug'.format(str(drug)))
    
    # get correlated drugs and indications
    correlated_drugs, correlated_indications = psm_functions.retrieve_features(drug_id=drug,
                                                                                 inds_corrs_db_path=correlations_indications,
                                                                                 drugs_corrs_db_path=correlations_drugs)
    
    # Exit function in case there are no significant features
    if len(correlated_drugs) == 0 and len(correlated_indications) == 0:
        message = '{} quitting function because no features were found'.format(str(drug))
        logging.info(message)
        return
    
    random_controls = dict()
    random_controls['primaryids'] = all_reports_primaryid - drug2reports[drug]['primaryids']
    random_controls['isrs'] = all_reports_isr - drug2reports[drug]['isrs']
    #limited_random_controls = psm_functions.limit_control_nr(random_controls) - used as trial as faster
    
    # The below controls are only for PSM and selected to have at least one of the correlated features
    # So not random
    cases, controls = psm_functions.define_case_control_reports(drug, indication2reports, drug2reports,
                                                                  correlated_indications, correlated_drugs)
    
    # Exit function if number of case reports is lower than threshold
    if len(cases['isrs']) + len(cases['primaryids']) < args.min_nr_cases:
        message = '{} quitting function because number of case reports below threshold of {}'.format(str(drug), str(args.min_nr_cases))
        logging.info(message)
        return
    
    # limit controls, define vectors
    limited_controls = psm_functions.limit_control_nr(controls)
    # 'valid' I believe meant that the report has at least one of the correlated features (as in Tatonetti ea)
    valid_case_reports, valid_control_reports, attributes_vector, label_vector = psm_functions.prepare_regression_input(drug, cases, limited_controls,
                                                                               report2attributes_primaryid,
                                                                               report2attributes_isr, correlated_drugs,
                                                                               correlated_indications)

    # get predicted scores
    predicted = psm_functions.return_predicted_by_model(attributes_vector, label_vector)
    
    # finished prepare objects elements separately
    
    # implement PS model, calculate propensity scores and make plots of PS distributions 
    if args.do_plots:
        psm_functions.do_ps_plots(drug, predicted, attributes_vector, label_vector, valid_case_reports, valid_control_reports, drug_names, savefig_path = savefig_path, data_path=data_path)
        logging.info('{} did ps model and plots'.format(str(drug)))
    
    # divide reports into bins and do sampling procedure
    sampled_controls, sampled_cases = psm_functions.sample_controls_and_cases(valid_case_reports, valid_control_reports, predicted, label_vector)
    
    if (len(sampled_cases['primaryids']) + len(sampled_cases['isrs'])) == 0 or (len(sampled_controls['primaryids']) + len(sampled_controls['isrs'])) == 0:
        logging.info('{} quitting function because not enough sampled reports'.format(str(drug)))
        return
    logging.info('{} did sampling'.format(str(drug)))
    
    # calculate differences in drug use and indications 
    differences_drugs = psm_functions.calculate_differences(correlated_list = correlated_drugs, item = 'drugs', names_dict = drug_names, 
                                                           top_correlated = 15, attributes_vector = attributes_vector, 
                                                           cases = cases, controls = random_controls, sampled_cases = sampled_cases, sampled_controls = sampled_controls, 
                                                           report2attributes_isr = report2attributes_isr, report2attributes_primaryid = report2attributes_primaryid)
    
    differences_indications = psm_functions.calculate_differences(correlated_list = correlated_indications, item = 'indications', names_dict = indication_names, 
                                                           top_correlated = 15, attributes_vector = attributes_vector, 
                                                           cases = cases, controls = random_controls, sampled_cases = sampled_cases, sampled_controls = sampled_controls, 
                                                           report2attributes_isr = report2attributes_isr, report2attributes_primaryid = report2attributes_primaryid)
    logging.info('{} calculated difference drugs and indications before and after'.format(str(drug)))
    
    # save plots of differences in drug use and indications
    if args.do_plots:
        psm_functions.plot_differences(drug = drug, item = 'drugs', differences_list = differences_drugs, savefig_path = savefig_path, data_path = data_path, drug_names = drug_names)
    
        psm_functions.plot_differences(drug = drug, item = 'indications', differences_list = differences_indications, savefig_path = savefig_path, data_path = data_path, drug_names = drug_names)
        logging.info('{} plotted differences'.format(str(drug)))
    
    
    # AE DETECTION WITH PROPENSITY SCORE MATCHING
    logging.info('{} starting PSM PRR calculations'.format(str(drug)))
    all_aes = psm_functions.combine_results(report_dict_for_cases=sampled_cases, report_dict_for_controls=sampled_controls, report2ae_isrs=report2ae_isrs, report2ae_primaryids=report2ae_primaryids, ae_concept_id2pt=ae_concept_id2pt, drug=drug)
    
    # Exit if there are no overlapping AEs between cases (exposed) and controls (non-exposed)
    if len(all_aes) == 0:
        logging.info('{} no AE results'.format(str(drug)))
    
    elif len(all_aes) > 0:
        logging.info('{} retrieved adverse events'.format(str(drug)))
    
        # save plot of AEs
        psm_functions.save_all_ae_results_csv(drug = drug, all_aes = all_aes, data_path = data_path, drug_names = drug_names, file_ending='_ae_associations.txt')
    
        logging.info('{} saved csv of AEs and finished drug'.format(str(drug)))

    # AE DETECTION WITHOUT PROPENSITY SCORE MATCHING (100,000/all randomly selected control reports)
    # retrieve AEs for reports
    logging.info('{} starting no PSM PRR calculations'.format(str(drug)))
    all_aes = psm_functions.combine_results(report_dict_for_cases=cases, report_dict_for_controls=random_controls, report2ae_isrs=report2ae_isrs, report2ae_primaryids=report2ae_primaryids, ae_concept_id2pt=ae_concept_id2pt, drug=drug)
    
    # Exit if there are no overlapping AEs between cases (exposed) and controls (non-exposed)
    if len(all_aes) == 0:
        logging.info('{} No AE results'.format(str(drug)))
    
    elif len(all_aes) > 0:
        logging.info('{} retrieved adverse events - no PSM'.format(str(drug)))
    
        # save plot of AEs
        psm_functions.save_all_ae_results_csv(drug = drug, all_aes = all_aes, data_path = data_path, drug_names = drug_names, file_ending='_no_PSM_ae_associations.txt')
    
        logging.info('{} saved csv of AEs and finished drug - no PSM'.format(str(drug)))
    
# Use multiprocessing
#701322(memantine), 1586226(cerivastatin)
if __name__ == '__main__':
    with multiprocessing.Pool(args.number_of_cores) as mypool:
        mypool.map(do_all_ae_tasks, aeolus_ids_result)

