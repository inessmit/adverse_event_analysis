#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""Extract adverse events, PRR, and chi squared value from PSM results, both with PSM and without, save those above PRR>3 and chi2 >4 to pickle, full data to df/txt file."""


# In[23]:


import argparse
import pandas as pd
import pickle
import os
import collections
import sqlite3 as sqlite
import re
import datetime
import itertools


# In[ ]:


# Command line argparse arguments

parser = argparse.ArgumentParser()
parser.add_argument("--project_dir", help="directory")
parser.add_argument("--mapped_compounds_db", type=str, help="database with mapped compounds (from compound_mapping dir)")
parser.add_argument("--psm_results_dir", type=str, help="directory with ae_associations.txt files from psm_aeolus")
parser.add_argument("--background_condition", type=str, help="describes which background reports used for no PSM situation")

args = parser.parse_args()
project_dir = args.project_dir


# In[ ]:


report_information = collections.OrderedDict()


# In[ ]:


# Construct aeolus2molregno dir

# Open sqlite database with aeolus to molregno mappings
conn = sqlite.connect(args.mapped_compounds_db)
cur = conn.cursor()

# Find aeolus to molregno mappings 
parent_query = """
select distinct aeolus_concept, mapped_parent_molregno
from compound_structures
where mapped_parent_molregno is not null
"""
parent_aeolus_ids = [i for i in cur.execute(parent_query).fetchall()]

conn.close()

# Set up dictionaries
aeolus2molregno = dict()
molregno2aeolus = dict()

# aeolus2molregno should be one molregno per aeolus only
for item in parent_aeolus_ids:
    aeolus2molregno[item[0]] = item[1]
    
# sometimes there are molregnos that are mapped to multiple aeolus_concepts
# Create a set of aeolus_ids per molregno
for item in parent_aeolus_ids:
    try: 
        molregno2aeolus[item[1]].add(item[0])
    except KeyError:
        molregno2aeolus[item[1]] = {item[0]}


# In[ ]:


report_information['Number of aeolus mapped to molregnos in mapping db'] = '{} aeolus ids mapped to {} molregnos'.format(len(aeolus2molregno), len(molregno2aeolus))


# In[ ]:


psm_ae_dict = dict()
no_psm_ae_dict = dict()

for molregno in molregno2aeolus.keys():
    psm_ae_dict[molregno] = {'aes': set(), 'dfs': []}

for molregno in molregno2aeolus.keys():
    no_psm_ae_dict[molregno] = {'aes': set(), 'dfs': []}


# In[ ]:


def extract_ae_data(filenames, report_situation, ae_dict):
    """Add AEs to dictionary per molregno, and dfs to list per molregno"""
    
    aeolus_ids_with_ae_file = []
    aeolus_ids_with_ae_file_but_no_PRRs_calculated = []

    for molregno in ae_dict.keys():
    
        aeolus_ids = molregno2aeolus[molregno]

        # iterate over the multiple aeolus_ids in case there are more than one and combine dataframes
        for aeolus_id in aeolus_ids:
            try:
                filename = [f for f in filenames if f.startswith(str(aeolus_id))][0]
            except IndexError:
                continue
            ae_df = pd.read_csv(args.psm_results_dir + '/{}'.format(filename), sep = '\t') 

            if len(ae_df.loc[~ae_df['PRR'].isnull()]) == 0:
                aeolus_ids_with_ae_file_but_no_PRRs_calculated.append(aeolus_id)

            if len(ae_df.loc[~ae_df['PRR'].isnull()]) > 0:
                aeolus_ids_with_ae_file.append(aeolus_id)
                ae_df['aeolus_id'] = aeolus_id
                ae_df['molregno'] = molregno
                ae_dict[molregno]['dfs'].append(ae_df.loc[~ae_df['PRR'].isnull()])

            for ae in ae_df.loc[(ae_df['PRR']>2) & (ae_df['chi-squared statistic']>4), 'adverse event']:
                ae_dict[molregno]['aes'].add(ae.upper())
    
    assert len(set(aeolus_ids_with_ae_file) & set(aeolus_ids_with_ae_file_but_no_PRRs_calculated)) == 0
    molregnos_with_ae_file = set([aeolus2molregno[int(i)] for i in aeolus_ids_with_ae_file])
    molregnos_without_PRR_calculated = set([aeolus2molregno[int(i)] for i in aeolus_ids_with_ae_file_but_no_PRRs_calculated])
        
    report_information['{}: Number of mapped aeolus/molregno with AE output file and PRRs calculated'.format(report_situation)] = '{} aeolus ids / {} molregnos'.format(len(aeolus_ids_with_ae_file), len(molregnos_with_ae_file))
    report_information['{}: Number of mapped aeolus/molregno with AE output file but no PRR calculated'.format(report_situation)] = '{} aeolus ids / {} molregnos'.format(len(aeolus_ids_with_ae_file_but_no_PRRs_calculated), len(molregnos_without_PRR_calculated))
        


# In[ ]:


# The way I select the filenames looks confusing but is correct
# I had no mentioning of 'PSM' as default for PSM situation, later added 'No PSM' to filename for that situation
filenames = [i for i in os.listdir(args.psm_results_dir) if ('ae_associations' in i and 'PSM' not in i)]
extract_ae_data(filenames=filenames, report_situation='PSM', ae_dict=psm_ae_dict)

filenames_no_psm = [i for i in os.listdir(args.psm_results_dir) if ('ae_associations' in i and 'PSM' in i)]
extract_ae_data(filenames=filenames_no_psm, report_situation='No PSM',ae_dict=no_psm_ae_dict)


# In[ ]:


# Remove items without data
psm_ae_dict_small = psm_ae_dict.copy()
for key, value in psm_ae_dict.items():
    if len(value['dfs']) == 0 and len(value['aes']) == 0:
        del(psm_ae_dict_small[key])
        
no_psm_ae_dict_small = no_psm_ae_dict.copy()
for key, value in no_psm_ae_dict.items():
    if len(value['dfs']) == 0 and len(value['aes']) == 0:
        del(no_psm_ae_dict_small[key])


# In[ ]:


# Save AE associations PSM to text file
df_lists = [psm_ae_dict_small[key]['dfs'] for key in psm_ae_dict_small.keys()] 
all_dfs_combined = pd.concat(list(itertools.chain(*df_lists)))
all_dfs_combined['adverse event'] = all_dfs_combined['adverse event'].map(lambda x: x.upper())
all_dfs_combined.to_csv(args.project_dir + f'/results/all_ae_info_PSM_{args.background_condition}.txt', sep='\t', index=False)

# Report counts PSM
psm_nr_molregnos = len(all_dfs_combined['molregno'].drop_duplicates())
psm_nr_aeolus_ids = len(all_dfs_combined['aeolus_id'].drop_duplicates())
psm_nr_aes = len(all_dfs_combined['adverse event'].drop_duplicates())
report_information['PSM all PRRs counts'] = f'number of molregnos: {psm_nr_molregnos}, aeolus_ids: {psm_nr_aeolus_ids}, adverse events: {psm_nr_aes}'

# Save AE associations no PSM to text file
df_lists = [no_psm_ae_dict_small[key]['dfs'] for key in no_psm_ae_dict_small.keys()] 
all_dfs_combined = pd.concat(list(itertools.chain(*df_lists)))
all_dfs_combined['adverse event'] = all_dfs_combined['adverse event'].map(lambda x: x.upper())
all_dfs_combined.to_csv(args.project_dir + f'/results/all_ae_info_no_PSM_{args.background_condition}.txt', sep='\t', index=False)

# Report counts no PSM
no_psm_nr_molregnos = len(all_dfs_combined['molregno'].drop_duplicates())
no_psm_nr_aeolus_ids = len(all_dfs_combined['aeolus_id'].drop_duplicates())
no_psm_nr_aes = len(all_dfs_combined['adverse event'].drop_duplicates())
report_information['No PSM all PRRs counts'] = f'number of molregnos: {no_psm_nr_molregnos}, aeolus_ids: {no_psm_nr_aeolus_ids}, adverse events: {no_psm_nr_aes}'


# In[ ]:


# Save AE dictionaries as pickles
psm_molregno2aes = dict()
for key,value in psm_ae_dict_small.items():
    psm_molregno2aes[key] = psm_ae_dict_small[key]['aes']
    
no_psm_molregno2aes = dict()
for key,value in no_psm_ae_dict_small.items():
    no_psm_molregno2aes[key] = no_psm_ae_dict_small[key]['aes']

current_date = datetime.date.today().strftime("%Y%m%d")
with open(args.project_dir + f'/results/{current_date}_PSM_molregno2aes_PRR3_chi4_faers_{args.background_condition}.pkl', 'wb') as f:
    pickle.dump(psm_molregno2aes, f)
    
with open(args.project_dir + f'/results/{current_date}_no_PSM_molregno2aes_PRR3_chi4_faers_{args.background_condition}.pkl', 'wb') as f:
    pickle.dump(no_psm_molregno2aes, f)

psm_all_aes =  set([i for i in itertools.chain(*[psm_molregno2aes[molregno] for molregno in psm_molregno2aes.keys()])])
report_information['PSM number of molregnos with significantly associated AEs before restricting to min5'] = '{}, linked to {} unique AEs'.format(len(psm_molregno2aes), len(psm_all_aes))
no_psm_all_aes =  set([i for i in itertools.chain(*[no_psm_molregno2aes[molregno] for molregno in no_psm_molregno2aes.keys()])])
report_information['No PSM number of molregnos with significantly associated AEs before restricting to min5'] = '{}, linked to {} unique AEs'.format(len(no_psm_molregno2aes), len(no_psm_all_aes))


# In[22]:


def restrict_min5drugs_per_ae(molregno2ae_dict):
    """Remove AEs with less than 5 significantly associated drugs and return new dictionary"""

    all_aes = set([i for i in itertools.chain(*[molregno2ae_dict[molregno] for molregno in molregno2ae_dict.keys()])])

    # Make new dictionary reversed and make empty set for each AE
    ae2molregnos = dict()
    for ae in all_aes:
        ae2molregnos[ae] = set()
    # Loop over molregnos in molregno2aes dictionary, and add the molregno to the ae2molregno dictionary
    for molregno in molregno2ae_dict:
        for ae in molregno2ae_dict[molregno]:
            ae2molregnos[ae].add(molregno)
    
    # Loop over ae2molregnos dictionary and find those AEs with less than 5 compounds
    aes_without_5_drugs = set()
    for ae in ae2molregnos:
        if len(ae2molregnos[ae]) < 5:
            aes_without_5_drugs.add(ae)
    
    # Make a new dictionary and remove those AEs with less than 5 compounds from the sets in molregno2ae
    restricted_molregno2ae_dict = molregno2ae_dict.copy()
    for ae in aes_without_5_drugs:
        for molregno in restricted_molregno2ae_dict:
            try:
                restricted_molregno2ae_dict[molregno].remove(ae)
            except KeyError:
                continue
    
    return restricted_molregno2ae_dict


# In[ ]:


# Restrict to AEs with min5 drugs
restricted_psm_molregno2aes = restrict_min5drugs_per_ae(molregno2ae_dict = psm_molregno2aes)
restricted_no_psm_molregno2aes = restrict_min5drugs_per_ae(molregno2ae_dict = no_psm_molregno2aes)

# Find and remove compounds with no significantly associated AEs
psm_no_info = [key for key, value in restricted_psm_molregno2aes.items() if len(value)==0]
no_psm_no_info = [key for key, value in restricted_no_psm_molregno2aes.items() if len(value)==0]

for molregno in psm_no_info:
    del(restricted_psm_molregno2aes[molregno])
for molregno in no_psm_no_info:
    del(restricted_no_psm_molregno2aes[molregno])

# Report info
all_aes_restricted_psm_molregno2aes =  set([i for i in itertools.chain(*[restricted_psm_molregno2aes[molregno] for molregno in restricted_psm_molregno2aes.keys()])])
report_information['PSM number of molregnos with significantly associated AEs after restricting to min5'] = '{}, linked to {} unique AEs'.format(len(restricted_psm_molregno2aes), len(all_aes_restricted_psm_molregno2aes))
all_aes_restricted_no_psm_molregno2aes =  set([i for i in itertools.chain(*[restricted_no_psm_molregno2aes[molregno] for molregno in restricted_no_psm_molregno2aes.keys()])])
report_information['No PSM number of molregnos with significantly associated AEs after restricting to min5'] = '{}, linked to {} unique AEs'.format(len(restricted_no_psm_molregno2aes), len(all_aes_restricted_no_psm_molregno2aes))

# Save as pickles
current_date = datetime.date.today().strftime("%Y%m%d")
with open(args.project_dir + f'/results/{current_date}_PSM_molregno2aes_PRR2_chi4_faers_min5drugs_{args.background_condition}.pkl', 'wb') as f:
    pickle.dump(restricted_psm_molregno2aes, f)
    
with open(args.project_dir + f'/results/{current_date}_no_PSM_molregno2aes_PRR2_chi4_faers_min5drugs_{args.background_condition}.pkl', 'wb') as f:
    pickle.dump(restricted_no_psm_molregno2aes, f)


# In[ ]:


# Write report with counts
with open(args.project_dir + f'/results/report_{args.background_condition}.txt', 'w') as f:
    for k, v in report_information.items():
        f.write(k + ': ' + v + '\n')


# In[24]:


# import pickle
# with open('/scratch/ias41/ae_code/faers_aes/results/20191204_no_PSM_molregno2aes_PRR3_chi4_faers.pkl', 'rb') as f:
#     no_psm_molregno2aes = pickle.load(f)
    
# with open('/scratch/ias41/ae_code/faers_aes/results/20191204_PSM_molregno2aes_PRR3_chi4_faers.pkl', 'rb') as f:
#     psm_molregno2aes = pickle.load(f)

