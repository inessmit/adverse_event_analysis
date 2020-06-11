#!/usr/bin/env python3

"""Find target-ae associations"""

import pandas as pd
import pickle
import itertools
import sklearn.metrics
import logging
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import datetime
import os
import argparse
import LR_functions_inactive as LR_functions
from statsmodels.stats import multitest
import multiprocessing
import multiprocessing_logging

# Argparse
parser = argparse.ArgumentParser()
parser.add_argument("-n", "--number_of_cores", default=1, type=int, help="number of cores to use for multiprocessing")
parser.add_argument("--min_nr_compounds_ae", default = 5, type=int, help="Minimum number of compounds with AE (to filter AE-target associations for consideration)")
parser.add_argument("--min_nr_compounds_active", default = 5, type=int, help="Minimum number of compounds active at a target (to filter AE-target associations for consideration)")
parser.add_argument("-t", "--targets", type=str, help="comma-separated uniprots to process, if none supplied will do all targets available")
parser.add_argument("-d", "--dirname", type=str, help='name for the output/experiment directory, will be newly made with date in name')
parser.add_argument("-ae_pickle", "--ae_pickle", type=str, help='path to molregno2ae pickle with aes')
parser.add_argument("-bioact_file", "--bioact_file", type=str, help='path to file of bioactivity data to use (integrated with plasma data already)')
parser.add_argument("-mtc", "--mtc_method", default='fdr_by', type=str, help="method for multiple testing correction for statsmodels multitest multipletests function, e.g. fdr_bh, fdr_by")
parser.add_argument('--assume_inactive', action="store_true", help="Provide this flag to assume negative bioactivity for unmeasured combinations")
parser.add_argument("-sign", "--sign_threshold", default=0.05, type=float)
parser.add_argument("-lr", "--lr_threshold", default=2, type=float)
parser.add_argument('--includes_predictions', action="store_true", help="Provide this flag if bioact_df includes predictions")

args = parser.parse_args()

# File locations
basedir = '/scratch/ias41/ae_code'

# Make directory for the outputs
experiment_name = args.dirname

current_date = datetime.date.today().strftime("%Y%m%d")
my_output_dir = basedir + '/ae_target_links/output/{}_{}'.format(current_date, experiment_name)

os.mkdir(my_output_dir)
os.mkdir(my_output_dir + '/plots')
os.mkdir(my_output_dir + '/data')
os.mkdir(my_output_dir + '/logs')

# Set up logging
logging.basicConfig(filename=my_output_dir + '/logs/{}_{}.log'.format(current_date, args.dirname), filemode='a', level=logging.INFO, format='%(asctime)s %(message)s')

# Use multiprocessing_logging module to collect logs from multiprocessing and write to same log file
multiprocessing_logging.install_mp_handler()

# Open drug2ae dictionaries
with open(args.ae_pickle, 'rb') as f:
    molregno2aes_all = pickle.load(f)
    
# Find compounds with no significantly associated AEs
no_info_molregnos = [molregno for molregno in molregno2aes_all if len(molregno2aes_all[molregno])==0]
nr_compounds_without_aes = len(no_info_molregnos)
# Restrict to drugs with at least one AE
molregno2aes = molregno2aes_all.copy()
for molregno in no_info_molregnos:
    del(molregno2aes[molregno])
assert nr_compounds_without_aes == len(molregno2aes_all) - len(molregno2aes)

# Open ChEMBL target names
target_names = pd.read_csv(basedir + '/ae_target_links/data/target_names.txt', sep='\t')

# Open bioactivity data
bioact_data = pd.read_csv(args.bioact_file, sep='\t')

# Assume inactivities?
if args.assume_inactive:
    bioact_data = bioact_data.loc[bioact_data['integrated_plasma_activity']==1,:]

# Restrict bioactivity to dataset to those compounds in the current dataset (sider/faers) and at least one AE
bioact_df = bioact_data.loc[bioact_data['parent_molregno'].isin(molregno2aes.keys()),:]

def do_LRs_mtc_plots(my_target):
    """For a target, loop over AEs and compute Likelihood Ratios, assemble results, do multiple testing correction on p-values (positive LRs only) and save as txt file. For highly significant results, save bar plot.
    kwargs: my_target -- str uniprot accession"""
    
    # If we are not assuming inactivity, then quit if all compounds are active
    # If all compounds are inactive it's already filtered out by having min_n active compounds filter
    if args.assume_inactive == False:
        inactive = bioact_df.loc[(bioact_df['accession']==my_target)&(bioact_df['integrated_plasma_activity']==0),:]
        if len(inactive) < 1:
            logging.info('{}: All compounds are active (no inactives), not analysed'.format(my_target))

    # Determine compounds with data for this target and AEs with data for these compounds
    molregnos_for_current_target = list(bioact_df.loc[bioact_df['accession']==my_target,'parent_molregno'])
    aes_for_these_drugs = set([i for i in itertools.chain(*[molregno2aes[molregno] for molregno in molregnos_for_current_target])])
    
    logging.info('{} : Now starting'.format(my_target))
    
    try:
        target_name = list(target_names.loc[target_names['accession']==my_target,'pref_name'])[0]
        organism_string = list(target_names.loc[target_names['accession']==my_target,'target_organism'])[0].replace(' ', '_')
    except IndexError:
        target_name = my_target # If not target name, set name equal to uniprot
        organism_string='Homo_sapiens'
    
    # Loop over AEs and append to results for this target
    data_list = []
    for my_ae in aes_for_these_drugs:
        data = LR_functions.compute_lr(my_target=my_target, bioact_df=bioact_df, ae_dict=molregno2aes, my_ae=my_ae, assume_inactive=args.assume_inactive, includes_predictions=args.includes_predictions)
        if data:
            data_list.append(data)
    
    # Check if there are results for any adverse events
    if len(data_list) == 0:
        logging.info('{}: No adverse events to report (not enough data)'.format(my_target))
        return 
    
    # If there are results, save text file, do multiple testing correction
    if len(data_list) >= 1:
        column_list = ['accession', 'nr compounds', 'nr compounds with AE','ae_hit_rate', 'nr compounds without AE', 'nae_hit_rate' , 'nr compounds active', 'nr compounds inactive', 'Adverse Event', 'Likelihood Ratio', 'p-value', 'activity_vector', 'ae_vector', 'molregnos', 'active_molregnos']
        if args.includes_predictions:
            column_list.append('predicted_vector')
        result_df = pd.DataFrame.from_records(data = data_list, columns = column_list)

        # Save version of the file with everything
        result_df.to_csv(my_output_dir + '/data/' + my_target + '_' + organism_string + '_likelihood_ratios_all.txt', sep = '\t', index=False)
                
        # Restrict to only AEs with positive LR before multiple testing correction
        result_df = result_df.loc[result_df['Likelihood Ratio']>1,:]
        if len(result_df) > 1:
            pvalues = result_df['p-value']
            multitest_results = multitest.multipletests(pvalues, alpha=0.25, method=args.mtc_method)
            result_df['corrected p-value'] = multitest_results[1]
            result_df.sort_values(by = 'corrected p-value', inplace=True)
            result_df.to_csv(my_output_dir + '/data/' + my_target + '_' + organism_string + '_likelihood_ratios.txt', sep = '\t', index=False)
            logging.info('{}: Finished LR and multiple testing correction'.format(my_target))    
        else:
            logging.info('{}: No adverse events to report (not positive LRs to analyse/not enough data)'.format(my_target))
            return
    
    # Do plots only for highly significant associations
    def do_plot(x):
        if x['Likelihood Ratio'] > 2.5 and x['corrected p-value'] < args.sign_threshold and x['nr compounds with AE']>=10 and x['nr compounds active']>=10:
            
            ae_vector = x['ae_vector']
            activity_vector = x['activity_vector']
            my_ae = x['Adverse Event']
            
            LR_functions.bar_plot(my_target=my_target, my_ae=my_ae, activity_vector = activity_vector, ae_vector=ae_vector, my_output_dir = my_output_dir, target_name = target_name, organism_string=organism_string)
        
    #result_df.apply(do_plot, axis=1)
    #logging.info('{}: Finished plots if any'.format(my_target))    

# Determine targets to loop over
if args.targets == None:
    targets = bioact_df['accession'].drop_duplicates()
elif args.targets:
    targets = [i for i in args.targets.split(',')]

# Do the LR calculation, plotting, multiple testing correction for current targets    
if __name__ == '__main__':
    with multiprocessing.Pool(args.number_of_cores) as mypool:
        mypool.map(do_LRs_mtc_plots, targets)

# With data file saved, get counts about the dataset of target-ae combinations tested with LR/included
# Gather all LR data
all_associations = LR_functions.find_all_associations(output_dir=my_output_dir, min_n=args.min_nr_compounds_active)
positive_associations, significant_associations = LR_functions.find_associations(output_dir=my_output_dir, min_n=args.min_nr_compounds_active, lr=args.lr_threshold, pv=args.sign_threshold, target_info=target_names)

# Dataset counts of all associations
with open(my_output_dir + '/dataset_counts.txt', 'a') as f:
    f.write('All associations for which LR was calculated (min_n={}):\n'.format(args.min_nr_compounds_active) + str(LR_functions.do_dataset_counts(all_associations, min_n=args.min_nr_compounds_active)))

# Dataset counts of positive associations
with open(my_output_dir + '/dataset_counts.txt', 'a') as f:
    f.write('\nAll positive associations for which LR was calculated (min_n={}):\n'.format(args.min_nr_compounds_active)+ str(LR_functions.do_dataset_counts(positive_associations, min_n=args.min_nr_compounds_active)))
    
# Save copy of all significant associations
file_info = 'sign_associations_LR{}_sign{}_minAE{}_minACT{}'.format(str(args.lr_threshold), str(args.sign_threshold), str(args.min_nr_compounds_ae), str(args.min_nr_compounds_active))
significant_associations.to_csv(my_output_dir + '/' + file_info + '.txt', sep='\t', index=False)
    
# Group significant associations by accession to get count of AE per target
target_counts = pd.DataFrame(significant_associations.groupby('accession')['Adverse Event'].count())
target_counts.reset_index(inplace=True, drop=False)
target_counts.columns = ['accession', 'Nr. of associated AEs']    
target_counts = target_counts.merge(target_names, on='accession', how='left')
target_counts = target_counts[['accession', 'target_organism', 'pref_name', 'Nr. of associated AEs']]
    
# Save copy of significant targets and counts
file_info2 = 'sign_targets_counts_LR{}_sign{}_minAE{}_minACT{}'.format(str(args.lr_threshold), str(args.sign_threshold), str(args.min_nr_compounds_ae), str(args.min_nr_compounds_active))
target_counts.sort_values(by='Nr. of associated AEs', ascending=False).to_csv(my_output_dir + '/' + file_info2 + '.txt', sep='\t', index=False)

# Share of predictions if applicable, for whole dataset and per target
if args.includes_predictions:
    with open(my_output_dir + '/share_of_predictions.txt', 'w') as f:
        f.write(str(LR_functions.calculate_share_of_predictions(bioact_df=bioact_df)))
    share_per_target_df = LR_functions.calculate_share_of_predictions_per_target(bioact_df=bioact_df)
    share_per_target_df.to_csv(my_output_dir + '/share_of_predictions_per_target.txt', sep='\t')
                
# Prepare general information about bioactivity and AE datasets (does not take account of overlap) used in this experiment and save to file        

# All unique AEs
nr_unique_aes = len(set([i for i in itertools.chain(*[molregno2aes[molregno] for molregno in molregno2aes.keys()])]))
# Nr of compounds with at least one AE
nr_compounds_with_aes = 0
for key in molregno2aes:
    if len(molregno2aes[key]) > 0:
        nr_compounds_with_aes += 1
    
info = [
    'Date: {}'.format(current_date)
    , 'Experiment name: {}'.format(experiment_name)
    , 'Adverse Event dataset: {}'.format(args.ae_pickle)
    , 'Bioactivity dataset: {}'.format(args.bioact_file)
    , 'Assume inactivities: {}'.format(args.assume_inactive)
    , 'Significance threshold: {}'.format(args.sign_threshold)
    , 'Likelihood Ratio threshold: {}'.format(args.lr_threshold)
    , 'Multiple testing correction method: {}'.format(args.mtc_method)
    , 'Minimum number of compounds active at a target: {}'.format(args.min_nr_compounds_active)
    , 'Minimum number of compounds with adverse event: {}'.format(args.min_nr_compounds_ae)
    , 'General bioactivity dataset - Number of unique compounds with at least one bioactivity: {} compounds'.format(len(bioact_df['parent_molregno'].drop_duplicates()))
    , 'General bioactivity dataset - Number of unique targets: {} targets'.format(len(bioact_df['accession'].drop_duplicates()))
    , 'General adverse event dataset - Number of unique adverse events: {} adverse events'.format(str(nr_unique_aes))
    , 'General adverse event dataset - Number of compounds with at least one adverse event: {} compounds'.format(str(nr_compounds_with_aes))
    , 'Number of compounds in AE dataset but without significantly associated AEs (were excluded): {}'.format(str(nr_compounds_without_aes)) 
]
with open(my_output_dir + '/conditions.txt', 'w') as f:
    f.write('\n'.join(info))
