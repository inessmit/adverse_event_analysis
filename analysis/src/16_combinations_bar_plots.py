#!/usr/bin/env python
# coding: utf-8

# In[77]:


import pandas as pd
import pickle
import itertools
import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:


basedir = '/scratch/ias41/ae_code'


# In[5]:


with open(basedir + '/analysis/data/dirs_info.pkl', 'rb') as f:
    dirs = pickle.load(f)


# In[100]:


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


# In[105]:


known_tuples_terms = set([(x[1]['Accession'], x[1][' Term']) for x in known_merged.iterrows()])
safety_targets = set(known_merged['Accession'])


# In[6]:


faers_data = dirs['20200110_faers_unbound_margin_pred_005_PRR2']
sider_data = dirs['20200110_sider_unbound_margin_pred']
combined_destination_dir = 'unbound_margin_pred_faers_vs_sider'


# In[9]:


# Save for later use
with open(basedir + '/analysis/data/uniprot2gene_id.pkl', 'rb') as f:
    uniprot2gene_id = pickle.load(f)


# In[31]:


faers_destination_dir = faers_data['dir']
faers_all = pd.read_csv(basedir + f'/ae_target_links/output/{faers_destination_dir}/combinations/all_target_combinations.txt', sep='\t')
faers_best = pd.read_csv(basedir + f'/ae_target_links/output/{faers_destination_dir}/combinations/best_target_combinations.txt', sep='\t')

sider_destination_dir = sider_data['dir']
sider_all = pd.read_csv(basedir + f'/ae_target_links/output/{sider_destination_dir}/combinations/all_target_combinations.txt', sep='\t')
sider_best = pd.read_csv(basedir + f'/ae_target_links/output/{sider_destination_dir}/combinations/best_target_combinations.txt', sep='\t')


# In[32]:


faers_best['Target set with highest recall'] = 1
faers_all['Target set with highest recall'] = 0
sider_best['Target set with highest recall'] = 1
sider_all['Target set with highest recall'] = 0


# In[33]:


sider_best['Targets'] = sider_best['Targets'].apply(lambda x: eval(x))
faers_best['Targets'] = faers_best['Targets'].apply(lambda x: eval(x))
sider_all['Targets'] = sider_all['Targets'].apply(lambda x: eval(x))
faers_all['Targets'] = faers_all['Targets'].apply(lambda x: eval(x))


# In[102]:


# Term to HLT dict
term2hlt = dict()
for row in meddra_hier_selection.iterrows():
    term2hlt[row[1][' Term']] = row[1]['HLT']


# In[180]:


plt.rcdefaults()
sns.set_style("whitegrid", {'axes.grid' : False, 'xtick.bottom': True, 'ytick.left': True})
plt.rcParams.update({'font.size': 14})
plt.rc('axes', labelsize=16)


# In[184]:


def make_combs_df(ae, best_df, all_df):
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
    
    #plt.rcdefaults()
    #plt.rcParams.update({'font.size': 13})
    #plt.rc('axes', labelsize=16)
    
    #fig = plt.figure(figsize=(7,7))
    #ax = fig.add_subplot(111)
    all_texts = []
        
    my_df = cumulative_performances_df.copy()
    if len(my_df)==1:
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
    
    fig = plt.figure(figsize=(7,7))

    compounds_screened_with_ae = list(my_df['Compounds screened with AE'].drop_duplicates())[0]
    compounds_screened = list(my_df['Compounds screened'].drop_duplicates())[0]

    my_df.rename(columns = {'Overall_recall': f"Fraction of AE-associated drugs that are active\n(out of n = {compounds_screened_with_ae})", 'Overall_PPV': 'Overall PPV'}, inplace=True)

    my_df.index = my_df['gene_ids']
    ax = my_df[['Overall PPV', f'Fraction of AE-associated drugs that are active\n(out of n = {compounds_screened_with_ae})']].plot.bar(rot=0, fontsize=12, color=['#99d8c9','#2ca25f'])

    # Define hatches
    for i,thisbar in enumerate(ax.patches):
        # Set a different hatch for each bar
        if i < (len(ax.patches)/2):
            thisbar.set_hatch('/')
        else:
            thisbar.set_hatch('')

    for p in ax.patches: 
        ax.annotate('{:.2}'.format(p.get_height()), (p.get_x()+p.get_width()/2., p.get_height()), ha='center', va='center', xytext=(0, 5), textcoords='offset points', fontsize=11)

    ax.legend(loc='upper right', bbox_to_anchor=(1, 1.3), ncol=1)

    ax.set_xlabel(f'Target (combinations)\n\nn = {compounds_screened} compounds overlapping')
    plt.tick_params(axis='both', which='major', labelsize=15)

    return plt.gcf()


# In[ ]:


for ae in list(faers_best['Adverse Event'].drop_duplicates()):
    myplot = make_combs_df(ae, faers_best, faers_all)
    ae_formatted = ae.replace(' ', '_').replace('/', '-')
    if myplot:
        myplot.savefig(basedir + f'/ae_target_links/output/{faers_destination_dir}/combinations/{ae_formatted}_bar.jpg', bbox_inches='tight')
    plt.close('all')


# In[ ]:


for ae in list(sider_best['Adverse Event'].drop_duplicates()):
    myplot = make_combs_df(ae, sider_best, sider_all)
    ae_formatted = ae.replace(' ', '_').replace('/', '-')
    if myplot:
        myplot.savefig(basedir + f'/ae_target_links/output/{sider_destination_dir}/combinations/{ae_formatted}_bar.jpg', bbox_inches='tight')
        plt.close()
    plt.close('all')


# In[ ]:




