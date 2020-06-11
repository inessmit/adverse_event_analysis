#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""After MedDRA hierarchy analysis on indication terms, extract HLT and retrieve other PTs/LLT in that HLT. Save dictionary with molregno to MedDRA PT/LLT indications, expanded via HLT. """


# In[ ]:


import pandas as pd
import argparse
import sqlite3 as sqlite
import pickle
import itertools


# In[ ]:


# Command line argparse arguments

parser = argparse.ArgumentParser()
parser.add_argument("--project_dir", help="directory")
parser.add_argument("--mapped_compounds_db", type=str, help="database with mapped compounds (from compound_mapping dir)")
parser.add_argument("--mesh2meddra_df")
parser.add_argument("--meddra_hier_output")
parser.add_argument("--all_aes_meddra_hier_file", type=str, help='meddra output file from hierarchy analysis on all adverse events in dataset')
parser.add_argument("--may_prevent_dict")
parser.add_argument('--may_treat_dict')

args = parser.parse_args()
project_dir = args.project_dir


# In[ ]:


# Open pickled indications

# Pickle dictionaries

with open(args.may_prevent_dict, 'rb') as f:
    may_prevent_rxnorm2disease = pickle.load(f)
with open(args.may_treat_dict, 'rb') as f:
    may_treat_rxnorm2disease = pickle.load(f)


# In[ ]:


# Open previous mesh2meddra file
mesh2meddra_df = pd.read_csv(args.mesh2meddra_df, sep='\t')


# In[ ]:


# Need to make Mesh ID to Meddra dictionary
# Checked that there is a meddra per mesh

mesh2meddra = dict()
for row in mesh2meddra_df[['MEDDRA','MESH_ID']].drop_duplicates().iterrows():
    meddra = row[1]['MEDDRA']
    mesh_id = row[1]['MESH_ID']
    mesh2meddra[mesh_id] = meddra


# In[ ]:


# Open meddra hierarchy analysis output
inds_hlgt = pd.read_excel(args.meddra_hier_output,skiprows=4)
inds_hlgt.columns = ['Row ID', 'Term', ' Code', 'Level', 'PT', 'PT Code', 'HLT', 'HLT Code',
       'HLGT', 'HLGT Code', 'SOC', 'SOC Code', 'Primary SOC']
inds_hlgt['Term'] = [i.upper() for i in inds_hlgt['Term']]


# In[ ]:


# Open sqlite database with rxnorm to molregno mappings
conn = sqlite.connect(args.mapped_compounds_db)
cur = conn.cursor()

# Find all RxNorm IDs that have a mapped molregno

compound_query = """select mapped_parent_molregno, rxnorm_concept
from compound_structures
where mapped_parent_molregno is not null
"""
molregno_rxnorm = cur.execute(compound_query).fetchall()

conn.close()


# In[ ]:


# There is one molregno for each Rxnorm id

rxnorm2molregno = dict()
for i in molregno_rxnorm:
    molregno = i[0]
    rxnorm = i[1]
    rxnorm2molregno[rxnorm] = molregno


# In[ ]:


# Open file with HLTs for all AE meddra terms
meddra_hlgt = pd.read_excel(args.all_aes_meddra_hier_file,skiprows=4)
meddra_hlgt.columns = ['Row ID', 'Term', ' Code', 'Level', 'PT', 'PT Code', 'HLT', 'HLT Code',
       'HLGT', 'HLGT Code', 'SOC', 'SOC Code', 'Primary SOC']
meddra_hlgt['Term'] = [i.upper() for i in meddra_hlgt['Term']]
meddra_hlgt['PT'] = [i.upper() for i in meddra_hlgt['PT']]


# In[ ]:


all_hlgt = pd.concat([inds_hlgt, meddra_hlgt])

term2hlt = dict()
hlt2terms = dict()
for hlt in set(all_hlgt['HLT'].drop_duplicates()):
    hlt2terms[hlt] = set()


# In[ ]:


# If there is a PT, add to dict
for row in all_hlgt.loc[~all_hlgt['PT'].isnull(),['Term','HLT']].iterrows():
    term = row[1]['Term']
    hlt = row[1]['HLT']
    term2hlt[term] = hlt
    hlt2terms[hlt].add(term)


# In[ ]:


# Sometimes the ind is a HLT, not PT
# In those cases add the occurring terms for that HLT from all_se file instead
for row in all_hlgt.loc[all_hlgt['PT'].isnull(),['Term','HLT']].iterrows():
    term = row[1]['Term']
    hlt = row[1]['HLT']
    term2hlt[term] = hlt
    
    occuring_terms = list(meddra_hlgt.loc[meddra_hlgt['HLT']==hlt,'Term'].drop_duplicates())
    for se_term in occuring_terms:
        hlt2terms[hlt].add(se_term)


# In[ ]:


# Start adding meddra inds


# In[ ]:


# Make empty molregno dictionary based on RxNorm ids with indications

rxnorm_ids_w_indications = may_prevent_rxnorm2disease.keys() | may_treat_rxnorm2disease.keys()
molregnos_w_indications = [rxnorm2molregno[rxnorm_id] for rxnorm_id in rxnorm_ids_w_indications]

molregno2meddra_inds_via_hlt = dict()
for molregno in molregnos_w_indications:
    molregno2meddra_inds_via_hlt[molregno] = set()


# In[ ]:


# Loop over indications and add to dictionary
def add_indications(my_dict):

    for rxnorm_id in my_dict.keys():
        molregno = rxnorm2molregno[rxnorm_id]
        mesh_indications = my_dict[rxnorm_id]

        meddra_indications = set()
        for ind in mesh_indications:
            try:
                meddra_ind = mesh2meddra[ind[1]]
                meddra_indications.add(meddra_ind)
            except KeyError:
                continue
        meddra_indications_hlt = [term2hlt[term] for term in meddra_indications]
        all_meddra_indication_terms = set([i for i in itertools.chain(*[hlt2terms[hlt] for hlt in meddra_indications_hlt])])
    
        for meddra_ind in all_meddra_indication_terms:
            molregno2meddra_inds_via_hlt[molregno].add(meddra_ind)

add_indications(my_dict=may_prevent_rxnorm2disease)
add_indications(my_dict=may_treat_rxnorm2disease)


# In[ ]:


with open(project_dir + '/results/molregno2meddra_inds_via_hlt.pkl', 'wb') as f:
    pickle.dump(molregno2meddra_inds_via_hlt, f)

