#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


# In[2]:


pd.set_option('display.max_rows',1000)


# In[3]:


basedir = '/scratch/ias41/ae_code'


# In[4]:


# applicability domain 70, performance filter TSSCV PRAUC 0.7


# In[4]:


tp = pd.read_csv(basedir + '/bioactivities/data/pidgin_input.smi_out_predictions_20200108-164123_ad70pr0.7.txt', sep='\t')
tp.set_index('Compound', inplace=True)


# In[8]:


median_plasma = pd.read_csv(basedir + '/plasma_concentrations/results/molregno2median_plasma_total_unbound.txt', sep='\t')


# In[11]:


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


# In[12]:


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


# In[13]:


def make_no_data_summary(x):
    if all(np.isnan(i) for i in [x[7], x[6], x[5], x[4]]):
        return 'no information'
    else:
        return np.nan
def find_negative_predictions(x):
    if any(is_negative_prediction(i) for i in [x[7], x[6], x[5], x[4]]):
        return 'negative prediction'
    else:
        return np.nan
def find_positive_predictions(x):
    if any(is_positive_prediction(i) for i in [x[7], x[6], x[5], x[4]]):
        return 'positive prediction'
    else:
        return np.nan


# In[68]:


# Restrict to those rows with at least one active/inactive prediction
tp_pivoted['no information'] = tp_pivoted.apply(make_no_data_summary, axis=1)
tp_pivoted = tp_pivoted.loc[tp_pivoted['no information'].isnull()]

# Identify rows with positive and rows with negative predictions
tp_pivoted['negative prediction'] = tp_pivoted.apply(find_negative_predictions, axis=1)
tp_pivoted['positive prediction'] = tp_pivoted.apply(find_positive_predictions, axis=1)


# In[74]:


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

def determine_lowest_active(x):
    
    for conc in [7,6,5,4]:
        if np.isnan(x[conc]):
            continue
        if x[conc] >= 0.6:
            lowest_predicted_active = conc
            return lowest_predicted_active    

def determine_highest_inactive(x):
    highest_predicted_inactive = np.nan
    
    for conc in [7,6,5,4]:
        if np.isnan(x[conc]):
            continue
        if x[conc] <= 0.4:
            highest_predicted_inactive = conc
    
    return highest_predicted_inactive


# 1. Exclude rows with 'not OK' trend (unreliable predicted active)
# 2. if active and inactive predictions: check lowest active prediction, should be lower than plasma conc. or within 1log unit range. Else inactive.
# 3. if only negative predictions: take highest predicted inactive, should be higher than plasma conc or in range. else no data
# 4. if only positive predictions: take lowest predicted active, should be lower than plasma conc or within 1log unit range. else: no data

# In[80]:


def integrate_plasma_conc_and_prediction(x, plasma_column, margin=True):
    
    # There are both positive and negative predictions depending on concentration
    if (x['positive prediction'] == 'positive prediction' and x['negative prediction'] == 'negative prediction'):
        
        # Lowest active prediction should be lower (higher pconc value) than plasma concentration or within 1 log unit range.
        # Logic is that if predicted active at a low concentration, should be active at higher concentration too.
        # Else inactive because there is a negative prediction, so any concentration lower should also be inactive.
        if x['lowest predicted active'] > x[plasma_column]:
            return 1
        elif x[plasma_column] - x['lowest predicted active'] < 1:
            if margin == True:
                return 1
            elif margin == False:
                return 0
        else:
            return 0
        
    # If there are only negative predictions
    if (x['positive prediction'] != 'positive prediction' and x['negative prediction'] == 'negative prediction'):
        # Highest predicted inactive should be higher (lower pconc value) than plasma conc or within 1 log unit range. Else no data.
        # Logic is that if predicted inactive at a concentration higher than plasma concentration, lower concentrations should be inactive too.
        if x['highest predicted inactive'] < x[plasma_column]:
            return 0
        elif x['highest predicted inactive'] - x[plasma_column] < 1:
            if margin == True: 
                return 0
            elif margin == False:
                return np.nan
        else:
            return np.nan
        
    # If there are only positive predictions
    if (x['positive prediction'] == 'positive prediction' and x['negative prediction'] != 'negative prediction'):
        # Lowest active prediction should be lower (higher pconc value) than plasma concentration or within 1 log unit range. Else 
        # Logic is that if predicted active at a low concentration, should be active at higher concentration too.


        if x['lowest predicted active'] > x[plasma_column]:
            return 1
        elif x[plasma_column] - x['lowest predicted active'] < 1:
            if margin == True:
                return 1
            elif margin == False:
                return np.nan
        else:
            return np.nan


# In[55]:


def do_target_prediction_plasma_conc_merge(tp_pivoted, plasma_concentrations_df, plasma_column):
    
    all_negative_predictions = tp_pivoted.loc[(~tp_pivoted['Compound'].isin(set(plasma_concentrations_df['molregno'])))&(tp_pivoted[7]<0.4)&(tp_pivoted[6]<0.4)&(tp_pivoted[5]<0.4)&(tp_pivoted[4]<0.4),:]
    all_negative_predictions.drop(labels='no information', axis=1, inplace=True)
    all_negative_predictions['integrated_plasma_activity'] = 0
    all_negative_predictions.rename(columns={'Compound': 'molregno'}, inplace=True)
    
    tp_merged = tp_pivoted.merge(plasma_concentrations_df, left_on='Compound', right_on='molregno')
    tp_merged.drop(labels=['no information'], axis=1, inplace=True)
    tp_merged['trend'] = tp_merged.apply(find_active_unreliable_predictions, axis=1)
    
    # What is % of rows 'not OK', counter-intuitive dose-response?'
    affected = len(tp_merged.loc[tp_merged['trend']=='not OK',:])
    total_len = len(tp_merged)
    percentage = (affected / total_len) * 100
    with open(basedir + '/integration_bioact_plasma_conc/results/not_ok_rows.txt', 'w') as f:
        f.write(f"Rows with 'not OK', counter-intuitive dose-response make up {affected} out of {total_len} rows, which is {percentage}%")

    #Exclude cases that are 'not ok' 
    tp_merged = tp_merged.loc[tp_merged['trend']!='not OK',:]
    
    tp_merged['lowest predicted active'] = tp_merged.apply(determine_lowest_active, axis=1)
    tp_merged['highest predicted inactive'] = tp_merged.apply(determine_highest_inactive, axis=1)
    
    tp_merged['prediction plasma integrated margin'] = tp_merged.apply(integrate_plasma_conc_and_prediction, plasma_column=plasma_column, margin=True, axis=1)
    tp_merged['prediction plasma integrated no margin'] = tp_merged.apply(integrate_plasma_conc_and_prediction, plasma_column=plasma_column, margin=False, axis=1)
    
    tp_merged_margin = tp_merged.drop(labels='prediction plasma integrated no margin', axis=1)
    tp_merged_margin.rename(columns={'prediction plasma integrated margin': 'integrated_plasma_activity'}, inplace=True)
    tp_merged_margin.dropna(subset=['integrated_plasma_activity'], axis=0, inplace=True)
    tp_merged_margin_plus_negative = pd.concat([tp_merged_margin, all_negative_predictions], sort=False)
    tp_merged_margin_plus_negative.rename(columns={'Target': 'accession'}, inplace=True)

    tp_merged_no_margin = tp_merged.drop(labels='prediction plasma integrated margin', axis=1)
    tp_merged_no_margin.rename(columns={'prediction plasma integrated no margin': 'integrated_plasma_activity'}, inplace=True)
    tp_merged_no_margin.dropna(subset=['integrated_plasma_activity'], axis=0, inplace=True)
    tp_merged_no_margin_plus_negative = pd.concat([tp_merged_no_margin, all_negative_predictions], sort=False)
    tp_merged_no_margin_plus_negative.rename(columns={'Target': 'accession'}, inplace=True)

    return tp_merged_margin_plus_negative, tp_merged_no_margin_plus_negative


# In[ ]:


plasma_total = median_plasma.loc[~median_plasma['median pMolar total plasma concentration'].isnull()]
integrated_total_margin, integrated_total_no_margin = do_target_prediction_plasma_conc_merge(tp_pivoted=tp_pivoted, plasma_concentrations_df=plasma_total, plasma_column='median pMolar total plasma concentration')

integrated_total_margin.to_csv(basedir + '/integration_bioact_plasma_conc/results/total_target_prediction_integrated_plasma_margin.txt', sep='\t', index=False)
integrated_total_no_margin.to_csv(basedir + '/integration_bioact_plasma_conc/results/total_target_prediction_integrated_plasma_no_margin.txt', sep='\t', index=False)

# Not needed any more
del(integrated_total_margin)
del(integrated_total_no_margin)


# In[ ]:


plasma_unbound = median_plasma.loc[~median_plasma['median pMolar unbound plasma concentration'].isnull()]
integrated_unbound_margin, integrated_unbound_no_margin = do_target_prediction_plasma_conc_merge(tp_pivoted=tp_pivoted, plasma_concentrations_df=plasma_unbound, plasma_column='median pMolar unbound plasma concentration')

integrated_unbound_margin.to_csv(basedir + '/integration_bioact_plasma_conc/results/unbound_target_prediction_integrated_plasma_margin.txt', sep='\t', index=False)
integrated_unbound_no_margin.to_csv(basedir + '/integration_bioact_plasma_conc/results/unbound_target_prediction_integrated_plasma_no_margin.txt', sep='\t', index=False)

