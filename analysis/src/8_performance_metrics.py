#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""In addition to PPV, calculate other performance metrics"""


# In[ ]:


import analysis_functions
import pickle
import pandas as pd
import sklearn.metrics


# In[ ]:


basedir = '/scratch/ias41/ae_code'


# In[ ]:


with open(basedir + '/analysis/data/dirs_info.pkl', 'rb') as f:
    dirs = pickle.load(f)


# In[ ]:


#faers_data = dirs['20200110_faers_unbound_margin_pred_005_PRR2']
#sider_data = dirs['20200110_sider_unbound_margin_pred']
faers_data = dirs['20200110_faers_cutoff6_pred_005_PRR2']
sider_data = dirs['20200110_sider_cutoff6_pred']


# In[ ]:


# Target information
target_info_file = basedir + '/ae_target_links/data/target_names.txt'
target_info = pd.read_csv(target_info_file, sep='\t')
target_info = target_info.loc[target_info['accession_organism']=='Homo sapiens',:]


# In[ ]:


# Load main info from directories
faers_pos, faers_sign =  analysis_functions.find_associations(basedir + '/ae_target_links/output/' + faers_data['dir'], min_n=faers_data['min_n'], lr=faers_data['lr'], pv=faers_data['pv'], target_info=target_info)
sider_pos, sider_sign = analysis_functions.find_associations(basedir + '/ae_target_links/output/' + sider_data['dir'], min_n=sider_data['min_n'], lr=sider_data['lr'], pv=sider_data['pv'], target_info=target_info)

for df in faers_pos, sider_pos, faers_sign, sider_sign:
    df['Adverse Event'] = df['Adverse Event'].apply(lambda x: x.upper())
    
# Calculate Positive predictive value
faers_pos['PPV'] = faers_pos.apply(lambda x: analysis_functions.calculate_ppv(x), axis=1)
sider_pos['PPV'] = sider_pos.apply(lambda x: analysis_functions.calculate_ppv(x), axis=1)
faers_sign['PPV'] = faers_sign.apply(lambda x: analysis_functions.calculate_ppv(x), axis=1)
sider_sign['PPV'] = sider_sign.apply(lambda x: analysis_functions.calculate_ppv(x), axis=1)


# In[ ]:


def calculate_prevalence(x):
    ae_prevalence = (x['nr compounds with AE'] / x['nr compounds'])
    return ae_prevalence


# In[ ]:


sider_pos['ae_prevalence'] = sider_pos.apply(calculate_prevalence, axis=1)
sider_sign['ae_prevalence'] = sider_sign.apply(calculate_prevalence, axis=1)
faers_pos['ae_prevalence'] = faers_pos.apply(calculate_prevalence, axis=1)
faers_sign['ae_prevalence'] = faers_sign.apply(calculate_prevalence, axis=1)


# In[ ]:


# FAERS and SIDER files molregnos are strings, need to turn back into lists to analyse drug overlap
def make_into_lists(row):
    return [int(float(i)) for i in row.strip('[]').split(', ')]


# In[ ]:


for df in [faers_pos, faers_sign, sider_pos, sider_sign]:
    for column in ['ae_vector', 'molregnos', 'active_molregnos', 'activity_vector']:
        df[column] = df[column].apply(lambda x: make_into_lists(x))


# In[ ]:


# AE_hit_rate = recall = sensitivity
# PPV = precision


# In[ ]:


def calculate_specificity(x):
    confusion_result = sklearn.metrics.confusion_matrix(x['ae_vector'], x['activity_vector']).ravel()
    tn, fp, fn, tp = confusion_result
   
    specificity = tn / (tn + fp)

    return specificity


# In[ ]:


def calculate_pru(x):
    PRU = (x['PPV'] - x['ae_prevalence']) / (1 - x['ae_prevalence'])
    return PRU


# In[ ]:


def calculate_improvement(x): 
    improvement =x['PPV'] - x['ae_prevalence']
    return improvement


# In[ ]:


for df in [faers_pos, faers_sign, sider_pos, sider_sign]:

    df['specificity'] = df.apply(calculate_specificity, axis=1)
    df['PRU'] = df.apply(calculate_pru, axis=1)
    df['improvement_over_prevalence'] = df.apply(calculate_improvement, axis=1)


# In[ ]:


sider_destination_dir = sider_data['dir']
faers_destination_dir = faers_data['dir']


# In[ ]:


faers_pos.to_csv(basedir + f'/ae_target_links/output/{faers_destination_dir}/pos_assoc_performance.txt', sep='\t', index=False)

faers_sign.to_csv(basedir + f'/ae_target_links/output/{faers_destination_dir}/sign_assoc_performance.txt', sep='\t', index=False)

sider_pos.to_csv(basedir + f'/ae_target_links/output/{sider_destination_dir}/pos_assoc_performance.txt', sep='\t', index=False)

sider_sign.to_csv(basedir + f'/ae_target_links/output/{sider_destination_dir}/sign_assoc_performance.txt', sep='\t', index=False)


# In[ ]:





# In[ ]:




