#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pickle
import argparse
import analysis_functions
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import itertools


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument("--faers_dir", type=str, help="output directory for FAERS")
parser.add_argument("--sider_dir", type=str, help="output directory for SIDER")
parser.add_argument("--dest_dir", type=str, help="destination directory for the figures, in basedir/analysis")

args = parser.parse_args()


# In[ ]:


basedir = '/scratch/ias41/ae_code'


# In[ ]:


destination_dir = basedir + f"/analysis/results/{args.dest_dir}"


# In[ ]:


with open(basedir + '/analysis/data/dirs_info.pkl', 'rb') as f:
    dirs = pickle.load(f)


# In[ ]:


# Files and info needed for the analysis

# Target information
target_info = pd.read_csv(basedir + '/ae_target_links/data/target_names.txt', sep='\t')
target_info = target_info.loc[target_info['accession_organism']=='Homo sapiens',:]

# Target classification
chembl_target_classification = pd.read_csv(basedir + '/analysis/data/target_classification_all_levels_r.txt', sep = '\t')

# MedDRA hierchy
meddra_hier = pd.read_excel(basedir + '/analysis/data/all_faers_and_sider_aes_hier_output.xlsx', skiprows=4)
meddra_hier_selection = meddra_hier.loc[meddra_hier['Primary SOC']=='Y',[' Term','HLT','SOC']].drop_duplicates()
meddra_hier_selection['HLT'] = meddra_hier_selection['HLT'].apply(lambda x: x.upper())

# ATC codes
all_atc_codes_loc = basedir + '/faers_sider_comparison/data/atc_all.txt'
small_molecule_atc_codes_loc = basedir + '/faers_sider_comparison/data/atc_small_molecules.txt'

# Previously reported associations
# Known associations, merge with known hierarchy HLT
known_associations = pd.read_excel(basedir + '/prev_reported_safety_associations/data/safety_meddra_annotated_effects.xls')
known_associations['Annotated MedDRA PT'] = known_associations['Annotated MedDRA PT'].apply(lambda x: x.upper())
known_meddra_hier = pd.read_excel(basedir + '/prev_reported_safety_associations/data/safety_meddra_annotated_effects_for_hierarchy_output.xlsx', skiprows=4)
known_meddra_hier['PT'] = known_meddra_hier['PT'].apply(lambda x: x.upper())
known_meddra_hier[' Term'] = known_meddra_hier[' Term'].apply(lambda x: x.upper())
known_meddra_hier['HLT'] = known_meddra_hier['HLT'].apply(lambda x: x.upper())
known_meddra_hier_selection = known_meddra_hier.loc[known_meddra_hier['Primary SOC']=='Y',['PT','HLT',' Term']].drop_duplicates()
known_merged = known_associations.merge(known_meddra_hier_selection, left_on='Annotated MedDRA PT', right_on=' Term')

hlt_manual = pd.read_excel(basedir + '/prev_reported_safety_associations/data/safety_meddra_manually_annotated_hlt_effects.xls', index=False)
hlt_manual.rename(columns={'Annotated MedDRA HLT': 'HLT'}, inplace=True)
hlt_manual['HLT'] = hlt_manual['HLT'].apply(lambda x: x.upper())
hlt_manual.drop(columns=['Annotated MedDRA HLT Code'])

known_merged = pd.concat([known_merged, hlt_manual], sort=False).reset_index(drop=True)
known_tuples = set([(x[1]['Accession'], x[1]['HLT']) for x in known_merged.iterrows()])


# In[ ]:


faers_data = dirs[args.faers_dir]
sider_data = dirs[args.sider_dir]

sider_all = analysis_functions.find_all_associations(basedir + '/ae_target_links/output/' + sider_data['dir'], min_n=sider_data['min_n'])
faers_all = analysis_functions.find_all_associations(basedir + '/ae_target_links/output/' + faers_data['dir'], min_n=faers_data['min_n'])

# Load main info from directories
faers_pos, faers_sign =  analysis_functions.find_associations(basedir + '/ae_target_links/output/' + faers_data['dir'], min_n=faers_data['min_n'], lr=faers_data['lr'], pv=faers_data['pv'], target_info=target_info)
sider_pos, sider_sign = analysis_functions.find_associations(basedir + '/ae_target_links/output/' + sider_data['dir'], min_n=sider_data['min_n'], lr=sider_data['lr'], pv=sider_data['pv'], target_info=target_info)

# Calculate Positive predictive value
faers_sign['PPV'] = faers_sign.apply(lambda x: analysis_functions.calculate_ppv(x), axis=1)
sider_sign['PPV'] = sider_sign.apply(lambda x: analysis_functions.calculate_ppv(x), axis=1)

# Find top associations
faers_top = faers_sign.loc[(faers_sign['Likelihood Ratio']>=faers_data['top_lr'])&(faers_sign['corrected p-value']<=faers_data['top_pv'])&(faers_sign['PPV']>=faers_data['ppv'])&(faers_sign['ae_hit_rate']>=faers_data['hr']),:]
sider_top = sider_sign.loc[(sider_sign['Likelihood Ratio']>=sider_data['top_lr'])&(sider_sign['corrected p-value']<=sider_data['top_pv'])&(sider_sign['PPV']>=sider_data['ppv'])&(sider_sign['ae_hit_rate']>=sider_data['hr']),:]


# In[ ]:


def calculate_prevalence(x):
    ae_prevalence = (x['nr compounds with AE'] / x['nr compounds'])
    return ae_prevalence


# In[ ]:


sider_sign['ae_prevalence'] = sider_sign.apply(calculate_prevalence, axis=1)
faers_sign['ae_prevalence'] = faers_sign.apply(calculate_prevalence, axis=1)


# In[ ]:


def calculate_improvement(x): 
    improvement =x['PPV'] - x['ae_prevalence']
    return improvement
for df in [faers_sign, sider_sign]:
    df['Value-added PPV'] = df.apply(calculate_improvement, axis=1)


# In[8]:


annotations_df = pd.read_csv(destination_dir + '/annotations_sign_va-PPV.txt', sep='\t')
annotations = []
for row in annotations_df.iterrows():
    annotations.append((row[1]['Annotation'], (float(row[1]['x']),float(row[1]['y'])), (float(row[1]['x_text']),float(row[1]['y_text']))))


# In[ ]:


annotations_custom_df = pd.read_csv(destination_dir + '/annotations_sign_custom_va-PPV.txt', sep='\t')
annotations_custom = []
for row in annotations_custom_df.iterrows():
    annotations_custom.append((row[1]['Annotation'], (float(row[1]['x']),float(row[1]['y'])), (float(row[1]['x_text']),float(row[1]['y_text']))))


# In[ ]:


analysis_functions.make_sign_target_overview_table(faers_top, 'FAERS', sider_top, 'SIDER', chembl_target_classification, known_merged, output_loc = destination_dir)
print('done comparison table of significant targets')

analysis_functions.do_hit_rate_ppv_plot(df2=faers_sign, df2_name='FAERS', df2_color='sandybrown', df1=sider_sign, df1_name='SIDER', df1_color='steelblue', output_dir=destination_dir, annotations=annotations, additional_filename='_significant')
analysis_functions.do_hit_rate_ppv_plot(df2=faers_sign, df2_name='FAERS', df2_color='sandybrown', df1=sider_sign, df1_name='SIDER', df1_color='steelblue', output_dir=destination_dir, additional_filename='_significant_no_annot')
analysis_functions.do_hit_rate_ppv_plot(df2=faers_sign, df2_name='FAERS', df2_color='sandybrown', df1=sider_sign, df1_name='SIDER', df1_color='steelblue', output_dir=destination_dir, annotations=annotations_custom, additional_filename='_significant_custom_annot')

print('done hit rate PPV plot')

analysis_functions.do_cumulative_recall(df1=faers_pos, df1_name='FAERS', df1_color='sandybrown', df2=sider_pos, df2_name='SIDER', df2_color='darkslateblue', meddra_hier=meddra_hier_selection, sort_by_columns=['corrected p-value','Likelihood Ratio'], ascending=[True,False], known_merged=known_merged, output_loc=destination_dir)
print('done cumulative recall')

analysis_functions.do_target_class_bar_plot(df1=sider_all, df1_name='SIDER', df1_color='steelblue', df2=faers_all, df2_name='FAERS', df2_color='sandybrown', chembl_target_classification=chembl_target_classification, safety_targets=set(known_merged['Accession']), df3_name='Previously reported safety targets', df3_color='darkgrey', output_loc=destination_dir)
analysis_functions.do_target_class_bar_plot_without_safety(df1=sider_all, df1_name='SIDER', df1_color='steelblue', df2=faers_all, df2_name='FAERS', df2_color='sandybrown', chembl_target_classification=chembl_target_classification, output_loc=destination_dir)
print('done target class bar plots')

analysis_functions.plot_actives_per_target(assoc_df1=faers_all, dataset1_name='FAERS', df1_color='sandybrown', assoc_df2=sider_all, dataset2_name='SIDER', df2_color='steelblue', output_loc=destination_dir)

with open(destination_dir + '/actives_count_describe.txt', 'w') as f:
    f.write('Number of active compounds per protein (FAERS):\n' + str(faers_pos[['accession','nr compounds active']].drop_duplicates().describe()) + f"\nmedian: {np.median(faers_pos[['accession','nr compounds active']]['nr compounds active'])}")
    f.write('\nNumber of active compounds per protein (SIDER):\n' + str(sider_pos[['accession','nr compounds active']].drop_duplicates().describe()) + f"\nmedian: {np.median(sider_pos[['accession','nr compounds active']]['nr compounds active'])}")

print('done histograms of actives per target')
            
analysis_functions.do_atc_bar_plot_three_datasets(assoc_df2=faers_all, df2_name='FAERS', df2_color='sandybrown', assoc_df1=sider_all, df1_name='SIDER', df1_color='steelblue', all_atc_codes_loc=all_atc_codes_loc, small_molecule_atc_codes_loc=small_molecule_atc_codes_loc, output_loc=destination_dir)
print('done combined ATC plot')            


# In[ ]:


# Write number of targets and significant AE-associations
sign_combined = pd.concat([faers_sign, sider_sign])
sign_combined['Adverse Event'] = sign_combined['Adverse Event'].apply(lambda x: x.upper())
nr_sign_aes = len(set(sign_combined['Adverse Event']))
nr_sign_targets = len(set(sign_combined['accession']))

top_combined = pd.concat([faers_top, sider_top])
top_combined['Adverse Event'] = top_combined['Adverse Event'].apply(lambda x: x.upper())
nr_top_aes = len(set(top_combined['Adverse Event']))
nr_top_targets = len(set(top_combined['accession']))

with open(destination_dir + '/nr_of_aes_counts.txt', 'w') as f:
    f.write(f"Number of unique significant targets: {nr_sign_targets}")
    f.write(f"\nNumber of unique significant AEs: {nr_sign_aes}")
    f.write(f"\nTop is defined as PPV => {faers_data['ppv']} (SIDER) or {sider_data['ppv']} (FAERS) and AE hit rate => {faers_data['hr']} (FAERS) or {sider_data['hr']} (SIDER).")
    f.write(f"\nNumber of unique top targets: {nr_top_targets}")
    f.write(f"\nNumber of unique top AEs: {nr_top_aes}")


# In[ ]:





# In[21]:


combined = pd.concat([faers_sign, sider_sign])


# In[22]:


combined['Adverse Event'] = combined['Adverse Event'].apply(lambda x: x.upper())


# In[23]:


len(set(combined['Adverse Event']))


# ### Do hit rate -PPV plot of known associations

# In[ ]:


# Calculate Positive predictive value
faers_merged = faers_sign.merge(meddra_hier_selection, left_on='Adverse Event', right_on=' Term')
sider_merged = sider_sign.merge(meddra_hier_selection, left_on='Adverse Event', right_on=' Term')

def find_known(row):
    if ((row['accession'],row['HLT'])) in known_tuples:
        return 1
    else:
        return 0
    
faers_merged['known'] = faers_merged.apply(find_known, axis=1)
sider_merged['known'] = sider_merged.apply(find_known, axis=1)

faers_known = faers_merged.loc[faers_merged['known']==1]
sider_known = sider_merged.loc[sider_merged['known']==1]


# In[ ]:


annotations_known_df = pd.read_csv(destination_dir + '/annotations_known_va-PPV.txt', sep='\t')
annotations_known = []
for row in annotations_known_df.iterrows():
    annotations_known.append((row[1]['Annotation'], (float(row[1]['x']),float(row[1]['y'])), (float(row[1]['x_text']),float(row[1]['y_text']))))


# In[ ]:


analysis_functions.do_hit_rate_ppv_plot(faers_known, 'FAERS', 'orange', sider_known, 'SIDER', 'blue', output_dir=destination_dir, annotations=annotations_known, additional_filename='_known', alpha=0.6)
analysis_functions.do_hit_rate_ppv_plot(faers_known, 'FAERS', 'orange', sider_known, 'SIDER', 'blue', output_dir=destination_dir, additional_filename='_known_no_annot', alpha=0.6)

