#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""Have measured bioactivity data integrated with plasma concentrations, and target predictions integrated with plasma concentrations.
Combine datasets, accepting predictions only where there are no measured bioactivities."""


# In[2]:


import pandas as pd


# In[3]:


basedir = '/scratch/ias41/ae_code'


# In[4]:


measured_unbound_margin = pd.read_csv(basedir + '/integration_bioact_plasma_conc/results/unbound_median_margin.txt', sep='\t')
measured_unbound_no_margin = pd.read_csv(basedir + '/integration_bioact_plasma_conc/results/unbound_median_no_margin.txt', sep='\t')
measured_total_margin = pd.read_csv(basedir + '/integration_bioact_plasma_conc/results/total_median_margin.txt', sep='\t')
measured_total_no_margin = pd.read_csv(basedir + '/integration_bioact_plasma_conc/results/total_median_no_margin.txt', sep='\t')


# In[7]:


predicted_unbound_margin = pd.read_csv(basedir + '/integration_bioact_plasma_conc/results/unbound_target_prediction_integrated_plasma_margin.txt', sep='\t')
predicted_unbound_no_margin = pd.read_csv(basedir + '/integration_bioact_plasma_conc/results/unbound_target_prediction_integrated_plasma_no_margin.txt', sep='\t')
predicted_total_margin = pd.read_csv(basedir + '/integration_bioact_plasma_conc/results/total_target_prediction_integrated_plasma_margin.txt', sep='\t')
predicted_total_no_margin = pd.read_csv(basedir + '/integration_bioact_plasma_conc/results/total_target_prediction_integrated_plasma_no_margin.txt', sep='\t')


# In[7]:


def combine_measured_and_predicted(measured_df, predicted_df):
    measured_targets = set(measured_df['accession'])
    
    #only include targets for which there is a measurement (all measured proteins only)
    predicted_selected = predicted_df.loc[(~predicted_df['integrated_plasma_activity'].isnull())&(predicted_df['accession'].isin(measured_targets)),:]
        
    measured_pairs = set([(i,j) for i,j in zip(measured_df['parent_molregno'],measured_df['accession'])])
    def find_measured(x):
        if (x['molregno'],x['accession']) in measured_pairs:
            return 1
        else:
            return 0
    
    predicted_selected['already measured'] = predicted_selected.apply(find_measured, axis=1)
    
    not_measured = predicted_selected.loc[predicted_selected['already measured']==0]
    not_measured.rename(columns={'molregno': 'parent_molregno'}, inplace=True)
    not_measured['predicted'] = 1
    
    measured_df_copy = measured_df.copy()
    measured_df['predicted'] = 0
    
    combined = pd.concat([measured_df, not_measured], sort=False)
    combined.drop(labels='already measured', axis=1, inplace=True)
    
    combined_selected = combined[['parent_molregno', 'accession', 'integrated_plasma_activity', 'pref_name','predicted']]
    return combined_selected


# In[53]:


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
        temp_df = group[1]
        target = group[0]
        if len(temp_df.loc[temp_df['integrated_plasma_activity']==1,:]) < 5:
            targets_without_5_active_compounds.append(target)
               
    chembl_plasma_margin_minimum5 = integrated_df.loc[~integrated_df['accession'].isin(targets_without_5_compounds),:]
    chembl_plasma_margin_minimum5active = integrated_df.loc[~integrated_df['accession'].isin(targets_without_5_active_compounds),:]
    
    return chembl_plasma_margin_minimum5, chembl_plasma_margin_minimum5active


# In[ ]:


# Do integration

combined_margin_unbound = combine_measured_and_predicted(measured_unbound_margin, predicted_unbound_margin)
combined_no_margin_unbound = combine_measured_and_predicted(measured_unbound_no_margin, predicted_unbound_no_margin)
combined_margin_total = combine_measured_and_predicted(measured_total_margin, predicted_total_margin)
combined_no_margin_total = combine_measured_and_predicted(measured_total_no_margin, predicted_total_no_margin)


unbound_median_margin_min5, unbound_median_margin_min5active = restrict_min_n(combined_margin_unbound)
unbound_median_no_margin_min5, unbound_median_no_margin_min5active = restrict_min_n(combined_no_margin_unbound)
total_median_margin_min5, total_median_margin_min5active = restrict_min_n(combined_margin_total)
total_median_no_margin_min5, total_median_no_margin_min5active = restrict_min_n(combined_no_margin_total)


# In[ ]:


# Save files without min5active
combined_margin_unbound.to_csv(basedir + '/integration_bioact_plasma_conc/results/unbound_median_margin_added_preds_measured_proteins.txt', sep = '\t', index=False)
combined_no_margin_unbound.to_csv(basedir + '/integration_bioact_plasma_conc/results/unbound_median_no_margin_added_preds_measured_proteins.txt', sep = '\t', index=False)
combined_margin_total.to_csv(basedir + '/integration_bioact_plasma_conc/results/total_median_margin_added_preds_measured_proteins.txt', sep = '\t', index=False)
combined_no_margin_total.to_csv(basedir + '/integration_bioact_plasma_conc/results/total_median_no_margin_added_preds_measured_proteins.txt', sep = '\t', index=False)


# In[ ]:


# Save files with min5active
unbound_median_margin_min5active.to_csv(basedir + '/integration_bioact_plasma_conc/results/unbound_median_margin_min5active_added_preds_measured_proteins.txt', sep = '\t', index=False)
unbound_median_no_margin_min5active.to_csv(basedir + '/integration_bioact_plasma_conc/results/unbound_median_no_margin_min5active_added_preds_measured_proteins.txt', sep = '\t', index=False)
total_median_margin_min5active.to_csv(basedir + '/integration_bioact_plasma_conc/results/total_median_margin_min5active_added_preds_measured_proteins.txt', sep = '\t', index=False)
total_median_no_margin_min5active.to_csv(basedir + '/integration_bioact_plasma_conc/results/total_median_no_margin_min5active_added_preds_measured_proteins.txt', sep = '\t', index=False)

