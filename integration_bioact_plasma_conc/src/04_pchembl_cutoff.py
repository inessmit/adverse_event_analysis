#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""To compare the effect of integrating with plasma concentrations, here we use cutoff of pchembl 6 instead."""


# In[7]:


import pandas as pd
import numpy as np


# In[8]:


basedir = '/scratch/ias41/ae_code'


# In[9]:


# Bioactivity data
median_bioact = pd.read_csv(basedir + '/bioactivities/results/bioact_medians_ae_drugs.txt', sep='\t')
bioact_slim = median_bioact[['parent_molregno', 'accession', 'summary']]


# In[10]:


def determine_activity(x):
    if x['summary'] == 'inactive':
        return 0
    elif float(x['summary']) >= 6:
        return 1
    else:
        return 0


# In[11]:


# This is not an integrated plasma activity but it's hardcoded in the script.. 
bioact_slim['integrated_plasma_activity'] = bioact_slim.apply(determine_activity, axis=1)


# In[12]:


def restrict_min_n(integrated_df):
    """Return copies of dataframe with bioactivities with less than 5 compounds and less than 5 active compounds removed.
    kwargs: integrated_df -- dataframe with integrated bioact&plasma concentration"""
    
    # Find which targets have less than 5 compounds associated
    targets_without_5_compounds = list()
    for group in integrated_df.groupby('accession'):
        if len(group[1]['parent_molregno'].drop_duplicates()) < 5:
            targets_without_5_compounds.append(group[0])
            
    # Find which targets have less than 5 active compounds associated
    targets_without_5_active_compounds = list()    
    for group in integrated_df.groupby('accession'):
        if len(group[1].loc[group[1]['integrated_plasma_activity']==1,:]) < 5:
            targets_without_5_active_compounds.append(group[0])
               
    chembl_plasma_margin_minimum5 = integrated_df.loc[~integrated_df['accession'].isin(targets_without_5_compounds),:]
    chembl_plasma_margin_minimum5active = integrated_df.loc[~integrated_df['accession'].isin(targets_without_5_active_compounds),:]
    
    return chembl_plasma_margin_minimum5, chembl_plasma_margin_minimum5active


# In[13]:


chembl_data_select_0 = bioact_slim[['accession','parent_molregno','integrated_plasma_activity']].drop_duplicates()
chembl_data_select_0['predicted'] = 0                                                                                                             


# In[14]:


pchembl_df_min5, pchembl_df_min5_active = restrict_min_n(chembl_data_select_0)


# In[15]:


len(chembl_data_select_0), len(pchembl_df_min5), len(pchembl_df_min5_active)


# In[15]:


# pchembl_df_min5_active.rename(columns={'molregno': 'parent_molregno'}, inplace=True)
# pchembl_df_min5_active.to_csv(basedir + '/ae_target_links/data/pchembl5_cutoff_median.txt', sep = '\t', index=False)


# In[ ]:


chembl_data_select_0.to_csv(basedir + '/integration_bioact_plasma_conc/results/pchembl6_cutoff_median.txt', sep = '\t', index=False)


# In[16]:


pchembl_df_min5_active.to_csv(basedir + '/integration_bioact_plasma_conc/results/pchembl6_cutoff_median_min5active.txt', sep = '\t', index=False)


# ### Target predictions

# In[91]:


measured_targets = set(bioact_slim['accession'])


# In[98]:


measured_pairs = set([(i,j) for i,j in zip(bioact_slim['parent_molregno'],bioact_slim['accession'])])
def find_measured(x):
    if (x['Compound'],x['Target']) in measured_pairs:
        return 1
    else:
        return 0 


# In[94]:


tp = pd.read_csv(basedir + '/bioactivities/data/pidgin_input.smi_out_predictions_20200108-164123_ad70pr0.7.txt', sep='\t')
tp.set_index('Compound', inplace=True)


# In[95]:


# Reformat dataframe

all_compound_target_combinations = [(compound,target) for compound in tp.index for target in set([i.split('_')[0] for i in tp.columns])]
conc_data = dict()
for item in all_compound_target_combinations:
    compound = item[0]
    target = item[1]
    conc_data[(compound,target)] = {'Target': target, 'Compound': compound, 7: np.nan, 6: np.nan, 5: np.nan, 4: np.nan}

pconc = {'0.1': 7, '1': 6, '10': 5, '100': 4}
for row in tp.iterrows():
    compound = row[0]
    for column, value in zip(row[1].index, row[1].values):
        target = column.split('_')[0]
        concentration = column.split('_')[1]
        conc_data[(compound,target)][pconc[concentration]] = value
tp_pivoted = pd.DataFrame(list(conc_data.values()))


# In[96]:


tp_pivoted.head()


# In[ ]:


# Restrict to targets with a measurement
tp_pivoted = tp_pivoted.loc[tp_pivoted['Target'].isin(measured_targets),:]


# In[ ]:


# Restrict to where no measurement is available
tp_pivoted['already measured'] = tp_pivoted.apply(find_measured, axis=1)


# In[ ]:


tp_pivoted = tp_pivoted.loc[tp_pivoted['already measured']==0]


# In[32]:


def make_no_data_summary(x):
    if all(np.isnan(i) for i in [x[7], x[6], x[5], x[4]]):
        return 'no information'
    else:
        return np.nan
# Restrict to those rows with at least one active/inactive prediction
tp_pivoted['no information'] = tp_pivoted.apply(make_no_data_summary, axis=1)
tp_pivoted = tp_pivoted.loc[tp_pivoted['no information'].isnull()]
tp_pivoted.drop(labels='no information', axis=1, inplace=True)


# In[30]:


def is_negative_prediction(x):
    if np.isnan(x):
        return False
    if x < 0.4:
        return True
    else: 
        return False
def is_positive_prediction(x):
    if np.isnan(x):
        return False
    if x > 0.6:
        return True
    else: 
        return False


# In[35]:


def find_active_unreliable_predictions(x):
    values = [i for i in [x[7], x[6], x[5], x[4]] if not np.isnan(i)]
    if len(values) == 1:
        return np.nan
    if any(is_positive_prediction(i) for i in values) and any(is_negative_prediction(i) for i in values):
        if [round(i,1) for i in values] == sorted([round(i,1) for i in values]):
            return 'OK'
        else:
            return 'not OK'
    else:
        return np.nan

tp_pivoted['trend'] = tp_pivoted.apply(find_active_unreliable_predictions, axis=1)


# In[38]:


#Exclude cases that are 'not ok'
tp_pivoted = tp_pivoted.loc[tp_pivoted['trend']!='not OK',:]


# In[41]:


tp_pivoted.drop(labels='trend', axis=1, inplace=True)


# In[42]:


tp_pivoted.head()


# In[44]:


tp_pivoted.loc[(tp_pivoted[7]>0.6)&(tp_pivoted[6]<0.4)]


# In[45]:


def determine_prediction_activity_call(x):
    if x[7] >= 0.6 or x[6] >= 0.6:
        return 1
    elif x[6] <= 0.4:
        return 0


# In[51]:


tp_pivoted['integrated_plasma_activity'] = tp_pivoted.apply(determine_prediction_activity_call, axis=1)


# In[ ]:


tp_select = tp_pivoted.loc[~tp_pivoted['integrated_plasma_activity'].isnull(),['Target','Compound','integrated_plasma_activity']].drop_duplicates()
tp_select.columns = ['accession', 'parent_molregno', 'integrated_plasma_activity']
tp_select['predicted'] = 1

chembl_data_select = bioact_slim[['accession','parent_molregno','integrated_plasma_activity']].drop_duplicates()
chembl_data_select['predicted'] = 0                                                                                                             

concat_df = pd.concat([chembl_data_select, tp_select])


# In[61]:


def restrict_min_n(integrated_df):
    """Return copies of dataframe with bioactivities with less than 5 compounds and less than 5 active compounds removed.
    kwargs: integrated_df -- dataframe with integrated bioact&plasma concentration"""
    
    # Find which targets have less than 5 compounds associated
    targets_without_5_compounds = list()
    for group in integrated_df.groupby('accession'):
        if len(group[1]['parent_molregno'].drop_duplicates()) < 5:
            targets_without_5_compounds.append(group[0])
            
    # Find which targets have less than 5 active compounds associated
    targets_without_5_active_compounds = list()    
    for group in integrated_df.groupby('accession'):
        if len(group[1].loc[group[1]['integrated_plasma_activity']==1,:]) < 5:
            targets_without_5_active_compounds.append(group[0])
               
    chembl_plasma_margin_minimum5 = integrated_df.loc[~integrated_df['accession'].isin(targets_without_5_compounds),:]
    chembl_plasma_margin_minimum5active = integrated_df.loc[~integrated_df['accession'].isin(targets_without_5_active_compounds),:]
    
    return chembl_plasma_margin_minimum5, chembl_plasma_margin_minimum5active


# In[ ]:


concat_df_min5, concat_df_min5active = restrict_min_n(concat_df)

concat_df_min5active.rename(columns={'molregno': 'parent_molregno'}, inplace=True)

concat_df.rename(columns={'molregno': 'parent_molregno'}, inplace=True)
concat_df.to_csv(basedir + '/integration_bioact_plasma_conc/results/pchembl6_tp_cutoff_median.txt', sep = '\t', index=False)

concat_df_min5active.to_csv(basedir + '/integration_bioact_plasma_conc/results/pchembl6_tp_cutoff_median_min5active.txt', sep = '\t', index=False)

