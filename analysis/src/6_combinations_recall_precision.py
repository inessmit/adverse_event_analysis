#!/usr/bin/env python
# coding: utf-8

# For the cases wich multiple targets/bioactivities available, try all combinations of associated targets and find highest overall recall and precision.

# In[ ]:


import pandas as pd
import pickle
import analysis_functions
import itertools
import urllib.parse
import urllib.request
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# In[ ]:


basedir = '/scratch/ias41/ae_code'


# In[ ]:


with open(basedir + '/analysis/data/dirs_info.pkl', 'rb') as f:
    dirs = pickle.load(f)


# In[ ]:


faers_data = dirs['20200110_faers_unbound_margin_pred_005_PRR2']
sider_data = dirs['20200110_sider_unbound_margin_pred']
combined_destination_dir = 'unbound_margin_pred_faers_vs_sider'


# In[ ]:


os.mkdir(basedir + '/ae_target_links/output/' + faers_data['dir'] + '/combinations')
os.mkdir(basedir + '/ae_target_links/output/' + sider_data['dir'] + '/combinations')


# In[ ]:


faers_pickle = basedir + '/faers_aes/results/20200108_PSM_molregno2aes_PRR2_chi4_faers_min5drugs_all_random_controls.pkl'
sider_pickle = basedir + '/sider/results/20191215_molregno2aes_sider_min5drugs.pkl'


# In[ ]:


target_info_file = basedir + '/ae_target_links/data/target_names.txt'
bioact_data_file = basedir + '/integration_bioact_plasma_conc/results/unbound_median_margin_min5active_added_preds_measured_proteins.txt'


# In[ ]:


# Open drug2ae dictionaries
with open(faers_pickle, 'rb') as f:
    molregno2aes_all = pickle.load(f)
    
# Find compounds with no significantly associated AEs
no_info_molregnos = [molregno for molregno in molregno2aes_all if len(molregno2aes_all[molregno])==0]
nr_compounds_without_aes = len(no_info_molregnos)
# Restrict to drugs with at least one AE
molregno2aes = molregno2aes_all.copy()
for molregno in no_info_molregnos:
    del(molregno2aes[molregno])
assert nr_compounds_without_aes == len(molregno2aes_all) - len(molregno2aes)

# Reverse dictionary
ae2molregno = {}
for molregno in molregno2aes:
    for AE in molregno2aes[molregno]:
        try:
            ae2molregno[AE].add(molregno)
        except KeyError:
            ae2molregno[AE] = {molregno}

with open(sider_pickle, 'rb') as f:
    molregno2aes_all_sider = pickle.load(f)
    
# Reverse dictionary
ae2molregno_sider = {}
for molregno in molregno2aes_all_sider:
    for AE in molregno2aes_all_sider[molregno]:
        try:
            ae2molregno_sider[AE].add(molregno)
        except KeyError:
            ae2molregno_sider[AE] = {molregno}


# In[ ]:


# Target information
target_info = pd.read_csv(target_info_file, sep='\t')
target_info = target_info.loc[target_info['accession_organism']=='Homo sapiens',:]


# In[ ]:


# Open bioactivity data
bioact_data = pd.read_csv(bioact_data_file, sep='\t')

# Restrict species
#bioact_data = bioact_data.loc[bioact_data['accession'].isin(set(target_info['accession'])),:]

# Restrict bioactivity to dataset to those compounds in the current dataset (sider/faers) and at least one AE
bioact_df = bioact_data.loc[bioact_data['parent_molregno'].isin(molregno2aes.keys()),:]
bioact_df_sider = bioact_data.loc[bioact_data['parent_molregno'].isin(ae2molregno_sider.keys()),:]


# In[ ]:


# MedDRA hierchy
meddra_hier = pd.read_excel(basedir + '/analysis/data/all_faers_and_sider_aes_hier_output.xlsx', skiprows=4)
meddra_hier_selection = meddra_hier.loc[meddra_hier['Primary SOC']=='Y',[' Term','HLT','SOC','PT']].drop_duplicates()
meddra_hier_selection['HLT'] = meddra_hier_selection['HLT'].apply(lambda x: x.upper())

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


known_tuples_terms = set([(x[1]['Accession'], x[1][' Term']) for x in known_merged.iterrows()])


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


faers_merged = faers_sign.merge(meddra_hier_selection, left_on='Adverse Event', right_on=' Term')
sider_merged = sider_sign.merge(meddra_hier_selection, left_on='Adverse Event', right_on=' Term')

def find_known(row):
    if ((row['accession'],row['Adverse Event'])) in known_tuples_terms:
        return 1
    else:
        return 0
    
faers_merged['known'] = faers_merged.apply(find_known, axis=1)
sider_merged['known'] = sider_merged.apply(find_known, axis=1)


# In[ ]:


# Get gene symbols for targets
accession_string = ' '.join(list(set(faers_merged['accession']) | set(sider_merged['accession'])))
url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'UNIPROTKB AC/ID',
'to': 'GENENAME',
'format': 'tab',
'query': f'{accession_string}'
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
    response = f.read().decode('utf-8')


# In[ ]:


uniprot2gene_id = dict()
for item in response.split('\n')[1:-1]:
    uniprot2gene_id[item.split('\t')[0]] = item.split('\t')[1]


# In[ ]:


# Save for later use
with open(basedir + '/analysis/data/uniprot2gene_id.pkl', 'wb') as f:
    pickle.dump(uniprot2gene_id, f)


# In[ ]:


# FAERS and SIDER files molregnos are strings, need to turn back into lists to analyse drug overlap
def make_into_lists(row):
    return [int(float(i)) for i in row.strip('[]').split(', ')]


# In[ ]:


faers_merged['ae_vector'] = faers_merged['ae_vector'].apply(lambda x: make_into_lists(x))
faers_merged['activity_vector'] = faers_merged['activity_vector'].apply(lambda x: make_into_lists(x))
sider_merged['ae_vector'] = sider_merged['ae_vector'].apply(lambda x: make_into_lists(x))
sider_merged['activity_vector'] = sider_merged['activity_vector'].apply(lambda x: make_into_lists(x))


# In[ ]:


faers_merged['molregnos'] = faers_merged['molregnos'].apply(lambda x: make_into_lists(x))
sider_merged['molregnos'] = sider_merged['molregnos'].apply(lambda x: make_into_lists(x))


# In[ ]:


faers_merged['active_molregnos'] = faers_merged['active_molregnos'].apply(lambda x: make_into_lists(x))
sider_merged['active_molregnos'] = sider_merged['active_molregnos'].apply(lambda x: make_into_lists(x))


# In[ ]:


faers_merged['kind'] = 'FAERS'
sider_merged['kind'] = 'SIDER'


# ## Iteratively try combinations

# In[ ]:


def try_all_combinations(AE, assoc_df, dataset, targets=None):
    
    if dataset=='FAERS':
        ae_molregnos = ae2molregno[AE]
    elif dataset=='SIDER':
        ae_molregnos = ae2molregno_sider[AE]

    relevant_df = assoc_df.loc[assoc_df['Adverse Event']==AE,:]
    if targets != None:
        relevant_df = assoc_df.loc[(assoc_df['Adverse Event']==AE)&(assoc_df['accession'].isin(targets)),:]
    
    # Find targets associated with AE
    target_list = list(relevant_df.loc[relevant_df['Adverse Event']==AE,:].sort_values(by='PPV', ascending=False)['accession'].drop_duplicates())
    all_compound_lists = [row[1]['molregnos'] for row in relevant_df.loc[relevant_df['Adverse Event']==AE].iterrows()]
    if len(all_compound_lists) == 0:
        print(AE + ': no overlap')
        return
    # Find compounds overlapping between targets
    overlapping_compounds = set.intersection(*map(set, all_compound_lists))
    if len(overlapping_compounds) == 0:
        print(AE + ': no overlap')
        return
      
    nr_compounds_ae_and_measured = len(overlapping_compounds&ae_molregnos)

    # Define all combinations
    combinations_to_test = set()
    len(target_list)
    
    for i in range(1, len(target_list)+1):
        combs = itertools.combinations(target_list, i)
        for comb in combs:
            combinations_to_test.add(comb)
    
    # Find performance of that combination
    combination_performances = dict()
    
    for target_set in combinations_to_test:
        target_df = bioact_df.loc[bioact_df['accession'].isin(target_set),:]
        active_molregnos = set(target_df.loc[target_df['integrated_plasma_activity']==1,'parent_molregno'])
        #molregnos = set(target_df['parent_molregno'])
        
        compounds_found = (active_molregnos & overlapping_compounds).intersection(ae_molregnos)
        false_positives = (active_molregnos & overlapping_compounds) - compounds_found
        assert len(compounds_found) + len(false_positives) == len(active_molregnos & overlapping_compounds)
    
        overall_recall = len(compounds_found) / nr_compounds_ae_and_measured
        overall_PPV = len(compounds_found) / len(active_molregnos & overlapping_compounds)
        
        combination_performances[', '.join(list(target_set))] = {'Targets': set(target_set), 'nr_targets': len(target_set), 'Adverse Event': AE, 'Overall_PPV': overall_PPV, 'Overall_recall': overall_recall, 'Compounds screened': len(overlapping_compounds), 'Compounds screened with AE': nr_compounds_ae_and_measured, 'Compounds found': len(compounds_found), 'False positives': len(false_positives)}
        
    # Find top row based on highest retrieval
    performances_df = pd.DataFrame.from_dict(combination_performances, orient='index').sort_values(by=['Overall_recall', 'Overall_PPV'], ascending=[False,False]).reset_index(drop=True)
    top_recall = performances_df.head(1)['Overall_recall'][0]
    corresponding_PPV = performances_df.head(1)['Overall_PPV'][0]
    top_performances = performances_df.loc[(performances_df['Overall_recall']==top_recall)&(performances_df['Overall_PPV']==corresponding_PPV)]

    # Find minimum sets with highest performance
    min_targets = top_performances['nr_targets'].min()
    combinations_of_interest = top_performances.loc[top_performances['nr_targets']==min_targets]
    
    # Return minimal set with highest performance
    return performances_df, combinations_of_interest


# In[ ]:


aes_w_multiple_targets_df = pd.DataFrame(faers_merged.groupby('Adverse Event')[['accession']].count().sort_values(by='accession', ascending=False).query('accession > 1'))
aes_w_multiple_targets_df_sider = pd.DataFrame(sider_merged.groupby('Adverse Event')[['accession']].count().sort_values(by='accession', ascending=False).query('accession > 1'))


# In[ ]:


sider_cumulative_all = []
sider_cumulative_best = []

for ae in aes_w_multiple_targets_df_sider.index:
    results = try_all_combinations(ae, assoc_df = sider_merged, dataset='SIDER')
    sider_cumulative_all.append(results[0])
    sider_cumulative_best.append(results[1])


# In[ ]:


faers_cumulative_all = []
faers_cumulative_best = []

for ae in aes_w_multiple_targets_df.index:
    results = try_all_combinations(ae, assoc_df = faers_merged, dataset='FAERS')
    faers_cumulative_all.append(results[0])
    faers_cumulative_best.append(results[1])


# In[ ]:


sider_all = pd.concat(sider_cumulative_all, ignore_index=True)
sider_best = pd.concat(sider_cumulative_best, ignore_index=True)
faers_all = pd.concat(faers_cumulative_all, ignore_index=True)
faers_best = pd.concat(faers_cumulative_best, ignore_index=True)


# In[ ]:


sider_destination_dir = sider_data['dir']
sider_all.to_csv(basedir + f'/ae_target_links/output/{sider_destination_dir}/combinations/all_target_combinations.txt', sep='\t', index=False)
sider_best.to_csv(basedir + f'/ae_target_links/output/{sider_destination_dir}/combinations/best_target_combinations.txt', sep='\t', index=False)

faers_destination_dir = faers_data['dir']
faers_all.to_csv(basedir + f'/ae_target_links/output/{faers_destination_dir}/combinations/all_target_combinations.txt', sep='\t', index=False)
faers_best.to_csv(basedir + f'/ae_target_links/output/{faers_destination_dir}/combinations/best_target_combinations.txt', sep='\t', index=False)


# In[ ]:


faers_all = pd.read_csv(basedir + f'/ae_target_links/output/{faers_destination_dir}/combinations/all_target_combinations.txt', sep='\t')

faers_best = pd.read_csv(basedir + f'/ae_target_links/output/{faers_destination_dir}/combinations/best_target_combinations.txt', sep='\t')

sider_all = pd.read_csv(basedir + f'/ae_target_links/output/{sider_destination_dir}/combinations/all_target_combinations.txt', sep='\t')

sider_best = pd.read_csv(basedir + f'/ae_target_links/output/{sider_destination_dir}/combinations/best_target_combinations.txt', sep='\t')


# In[ ]:


faers_best['Target set with highest recall'] = 1
faers_all['Target set with highest recall'] = 0
sider_best['Target set with highest recall'] = 1
sider_all['Target set with highest recall'] = 0


# In[ ]:


faers_best.head()


# In[ ]:


# The 'keep='first' only works because I re-imported the file and then the Targets is a string instead of a set... 
# Combine best and all other combinations FAERS into one file
faers_integrated = pd.concat([faers_best, faers_all], sort=False).sort_values(by='Target set with highest recall', ascending=False)
faers_integrated2 = faers_integrated.drop_duplicates(subset=['Targets', 'nr_targets', 'Adverse Event', 'Overall_PPV','Overall_recall', 'Compounds screened', 'Compounds screened with AE','Compounds found', 'False positives'], keep='first')

faers_integrated2.fillna({'Target set with highest recall': 0}, inplace=True)
faers_integrated2.sort_values(by=['Adverse Event','Target set with highest recall'], ascending=[True,False], inplace=True)
faers_integrated2['Target set with highest recall'] = faers_integrated2['Target set with highest recall'].astype('int')

# Save file
faers_integrated2.to_csv(basedir + '/ae_target_links/output/' + faers_data['dir'] + '/combinations/all_combinations_incl_best.txt', sep='\t', index=False)


# In[ ]:


# Combine best and all other combinations SIDER into one file
sider_integrated = pd.concat([sider_best, sider_all], sort=False).sort_values(by='Target set with highest recall')
sider_integrated2 = sider_integrated.drop_duplicates(subset=['Targets', 'nr_targets', 'Adverse Event', 'Overall_PPV','Overall_recall', 'Compounds screened', 'Compounds screened with AE','Compounds found', 'False positives'], keep='first')
sider_integrated2.fillna({'Target set with highest recall': 0}, inplace=True)
sider_integrated2.sort_values(by=['Adverse Event','Target set with highest recall'], ascending=[True,False], inplace=True)
sider_integrated2['Target set with highest recall'] = sider_integrated2['Target set with highest recall'].astype('int')


# In[ ]:


# Save file
sider_integrated2.to_csv(basedir + '/ae_target_links/output/' + sider_data['dir'] + '/combinations/all_combinations_incl_best.txt', sep='\t', index=False)


# In[ ]:


# Save 10 best combinations across datasets

sider_integrated3 = sider_best.copy()
faers_integrated3 = faers_best.copy()
sider_integrated3['Dataset'] = 'S'
faers_integrated3['Dataset'] = 'F'

both_datasets = pd.concat([faers_integrated3, sider_integrated3], sort=False)

my_selection = both_datasets.loc[both_datasets['nr_targets']>1].sort_values(by=['Overall_recall','Overall_PPV'], ascending=False).drop_duplicates(subset=['Adverse Event'], keep='first')[:10]
my_selection['Target set'] = my_selection['Targets'].apply(lambda x: ', '.join([uniprot2gene_id[i] for i in eval(x)]))
my_selection.rename(columns={'Overall_PPV': 'Overall PPV', 'Overall_recall': 'Overall recall'}, inplace=True)
my_selection['Adverse Event'] = my_selection['Adverse Event'].apply(lambda x: x.lower().capitalize())
my_selection['Overall recall'] = my_selection['Overall recall'].apply(lambda x: '{:.2f}'.format(float(x)))
my_selection['Overall PPV'] = my_selection['Overall PPV'].apply(lambda x: '{:.2f}'.format(float(x)))

my_selection[['Target set', 'Adverse Event', 'Overall PPV',
       'Overall recall','Dataset']].to_csv(basedir + f'/analysis/results/{combined_destination_dir}/top10_recall_combinations.txt', sep='\t', index=False)


# ### Do plots

# In[ ]:


def make_combs_df(ae, best_df, all_df, legend_loc='best'):
    cumulative_performances = list()
    
    best_combination = list(best_df.loc[best_df['Adverse Event']==ae,'Targets'])[0]
    
    combinations_to_test = set()
    for i in range(1, len(best_combination)):
        combs = itertools.combinations(best_combination, i)
        for comb in combs:
            combinations_to_test.add(comb)
    
    for comb in combinations_to_test:
        comb_df = all_df.loc[(all_df['Adverse Event']==ae)&(all_df['Targets']==set(comb))]
        cumulative_performances.append(comb_df)
    
    cumulative_performances.append(best_df.loc[best_df['Adverse Event']==ae])
    cumulative_performances_df = pd.concat(cumulative_performances, ignore_index=True).sort_values(by=['Overall_recall', 'Overall_PPV']).reset_index(drop=True).reset_index(drop=False)
    
    mystyles = iter(['solid', 'dashed', 'dotted', (0, (3, 1, 1, 1)), (0, (3, 5, 1, 5)), (0, (1, 10))])
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 13})
    plt.rc('axes', labelsize=16)
    
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111)
    all_texts = []
        
    my_df = cumulative_performances_df.copy()
    if len(my_df) ==1:
        return
    
    def find_superscript_symbol(accession, ae):
        gene_id = uniprot2gene_id[accession]
        
        ae_hlt = term2hlt[ae]
        superscript1 = ''
        superscript2 = ''
        if (accession, ae_hlt) in known_tuples:
            superscript1 = r'$^\blacktriangledown$'
        if accession in safety_targets:
            superscript2 = r'$^\triangledown$'
        return gene_id + superscript1 + superscript2 
            
    my_df['gene_ids'] = my_df['Targets'].apply(lambda x: ', '.join([find_superscript_symbol(y,ae) for y in x]))
    
    compounds_screened_with_ae = list(my_df['Compounds screened with AE'].drop_duplicates())[0]
    compounds_screened = list(my_df['Compounds screened'].drop_duplicates())[0]

    plt.plot(my_df['index'], my_df['Overall_recall'], '-o', markersize=3,linestyle=next(mystyles), color='black', label=f'Recall (out of n = {compounds_screened_with_ae})')
    plt.plot(my_df['index'], my_df['Overall_PPV'], '-o', markersize=3,linestyle=next(mystyles), color='black', label='Positive Predictive Value')

    gene_ids = list(my_df['gene_ids'])#[uniprot2gene_id[accession] for accession in list(my_df['Targets'])]

    #colors = ['blue' if x==1 else 'black' if y==1 else 'red' for x,y in zip(list(my_df['known association']), list(my_df['known target']))]
#     texts = [plt.text(x,y,z, ha='center', va='center', size=12) for x,y,z in zip(my_df['index'], my_df['Overall_recall'],gene_ids)]
#     adjust_text(texts, arrowprops=dict(arrowstyle='-', color='darkgrey'))

#     for text in texts:
#         all_texts.append(text)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    if legend_loc ==None:
        plt.legend(bbox_to_anchor=(1,1))
    else:
        plt.legend(loc=legend_loc)
    ax.set_ylabel('Performance')
    ax.set_xlabel(f'Target contributions\nn = {compounds_screened} compounds overlapping')
    
    ax.set_xticks(list(my_df['index']))
    ax.set_xticklabels(gene_ids, rotation=45,horizontalalignment='right')

    plt.tick_params(axis='both', which='major', labelsize=15)

    return plt.gcf()


# In[ ]:


# Term to HLT dict
term2hlt = dict()
for row in meddra_hier_selection.iterrows():
    term2hlt[row[1][' Term']] = row[1]['HLT']


# In[ ]:


safety_targets = set(known_merged['Accession'])


# In[ ]:


sider_best['Targets'] = sider_best['Targets'].apply(lambda x: eval(x))
faers_best['Targets'] = faers_best['Targets'].apply(lambda x: eval(x))
sider_all['Targets'] = sider_all['Targets'].apply(lambda x: eval(x))
faers_all['Targets'] = faers_all['Targets'].apply(lambda x: eval(x))


# In[ ]:


for ae in list(sider_best['Adverse Event'].drop_duplicates()):
    myplot = make_combs_df(ae, sider_best, sider_all)
    ae_formatted = ae.replace(' ', '_').replace('/', '-')
    if myplot:
        myplot.savefig(basedir + f'/ae_target_links/output/{sider_destination_dir}/combinations/{ae_formatted}.jpg', bbox_inches='tight')
        plt.close()
    plt.close('all')


# In[ ]:


for ae in list(faers_best['Adverse Event'].drop_duplicates()):
    myplot = make_combs_df(ae, faers_best, faers_all)
    ae_formatted = ae.replace(' ', '_').replace('/', '-')
    if myplot:
        myplot.savefig(basedir + f'/ae_target_links/output/{faers_destination_dir}/combinations/{ae_formatted}.jpg', bbox_inches='tight')
    plt.close('all')


# In[ ]:




