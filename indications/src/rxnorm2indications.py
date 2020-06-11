#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import sqlite3 as sqlite
import argparse
import requests
import time
import pickle
import datetime
import itertools
from pandas import ExcelWriter


# In[ ]:


# Command line argparse arguments

parser = argparse.ArgumentParser()
parser.add_argument("--project_dir", help="directory")
parser.add_argument("--mapped_compounds_db", type=str, help="database with mapped compounds (from compound_mapping dir)")
parser.add_argument("--mrconso_file", type=str, help="MRCONSO.RRF file downloaded from UMLS")

args = parser.parse_args()
project_dir = args.project_dir


# In[ ]:


# Open sqlite database with rxnorm to molregno mappings
conn = sqlite.connect(args.mapped_compounds_db)
cur = conn.cursor()

# Find all RxNorm IDs that have a mapped molregno
compound_query = """select distinct rxnorm_concept
from compound_structures
where mapped_parent_molregno is not null
"""
rxnorm_ids = [i[0] for i in cur.execute(compound_query).fetchall()]
print('Nr of RxNorm ids: {}'.format(len(rxnorm_ids)))
conn.close()


# In[ ]:


# RxNorm version in the current RxNorm API # Ideally this is the same as the RxNorn downloaded files
response = requests.get('https://rxnav.nlm.nih.gov/REST/version.json')

with open(project_dir + '/results/RxNorm_API_version.txt', 'w') as f:
    f.write(response.content.decode())


# In[ ]:


base = 'https://rxnav.nlm.nih.gov/REST/rxclass'
def do_request(cui, rel):
    r = requests.get(base + f'/class/byRxcui.json?rxcui={cui}&relaSource=MEDRT&relas={rel}')
    if r.status_code == 200:
        return r
    else:
        print(cui, r.status_code)


# In[ ]:


may_treat_rxnorm2disease = dict()

for cui in rxnorm_ids:
    r = do_request(cui, rel='may_treat')
    time.sleep(0.1)

    indications = set()

    try:
        for item in r.json()['rxclassDrugInfoList']['rxclassDrugInfo']:
            disease = item['rxclassMinConceptItem']['className']
            mesh_id = item['rxclassMinConceptItem']['classId']
            indications.add((disease, mesh_id))
        if len(indications) > 0:
            may_treat_rxnorm2disease[cui] = indications
    except KeyError:
        continue


# In[ ]:


print(datetime.datetime.now().strftime("%H:%M:%S") + ' Finished retrieving "may_treat" relationships')


# In[ ]:


may_prevent_rxnorm2disease = dict()

for cui in rxnorm_ids:
    r = do_request(cui, rel='may_prevent')
    time.sleep(0.1)

    indications = set()

    try:
        for item in r.json()['rxclassDrugInfoList']['rxclassDrugInfo']:
            disease = item['rxclassMinConceptItem']['className']
            mesh_id = item['rxclassMinConceptItem']['classId']
            indications.add((disease, mesh_id))
        if len(indications) > 0:
            may_prevent_rxnorm2disease[cui] = indications
    except KeyError:
        continue


# In[ ]:


print(datetime.datetime.now().strftime("%H:%M:%S") + ' Finished retrieving "may_prevent" relationships')


# In[ ]:


# Pickle dictionaries

with open(project_dir + '/data/rxnorm2may_treat.pkl', 'wb') as f:
    pickle.dump(may_treat_rxnorm2disease, f)
with open(project_dir + '/data/rxnorm2may_prevent.pkl', 'wb') as f:
    pickle.dump(may_prevent_rxnorm2disease, f)
print('Pickled dictionaries in /data...')


# In[ ]:


print(datetime.datetime.now().strftime("%H:%M:%S") + ' Starting to load MRCONSO file... ')
# Open MRCONSO file from UMLS, to map between vocabularies
umls_conso_loc = args.mrconso_file
umls_conso = pd.read_csv(umls_conso_loc, sep='|', header=None, names=['CUI','LAT', 'TS', 'LUI', 'STT', 'SUI', 'ISPREF', 'AUI', 'SAUI', 'SCUI', 'SDUI', 'SAB', 'TTY', 'CODE', 'STR', 'SRL', 'SUPPRESS', 'CVF'], index_col=False)
umls_conso['STR'] = umls_conso['STR'].apply(lambda x: str(x).upper())


# In[ ]:


# Locate the rows with indications for the compounds above
all_indications = set([i[1] for i in itertools.chain(*may_prevent_rxnorm2disease.values())]) | set([i[1] for i in itertools.chain(*may_treat_rxnorm2disease.values())])
umls_indications = umls_conso.loc[(umls_conso['SAB']=='MSH')&(umls_conso['SDUI'].isin(all_indications)),['CUI','STR', 'CODE','TTY']].drop_duplicates()
umls_indications.columns = ['CUI', 'MESH', 'MESH_ID', 'MESH TYPE']


# In[ ]:


# Find MedDRA terms corresponding to downloaded indications in MeSH terminology
umls_indications_meddra = umls_conso.loc[(umls_conso['SAB']=='MDR')&(umls_conso['CUI'].isin(umls_indications['CUI'].drop_duplicates())),['CUI','STR','SDUI','TTY']].drop_duplicates()
umls_indications_meddra.columns = ['CUI', 'MEDDRA','MEDDRA_ID', 'MEDDRA TYPE']
conso_merged = umls_indications_meddra.merge(umls_indications, on='CUI', suffixes = ['',''])


# In[ ]:


conso_merged.to_csv(project_dir + '/data/mesh2meddra_df.txt', sep='\t', index=False)
print(datetime.datetime.now().strftime("%H:%M:%S") + ' Wrote file with MeSH and MedDRA terms to /data.')


# In[ ]:


# Need to do MedDRA hierarchy analysis for all these terms
inds_df = pd.DataFrame({'Term': list(conso_merged['MEDDRA'].drop_duplicates())})
inds_df.reset_index(inplace=True)
inds_df.columns = ['Row ID', 'Term']


# In[ ]:


writer = ExcelWriter(project_dir + '/data/meddra_terms_indications_hierarchy_input.xlsx')
inds_df.to_excel(writer,'Sheet1', index=False)
writer.save()
print(datetime.datetime.now().strftime("%H:%M:%S") + ' Wrote input file with MedDRA terms (total {}) to /data for hierarchy analysis.'.format(len(inds_df)))

