#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""Get all MedDRA terms occurring in FAERS AEOLUS and save to excel file for hierarchy analysis in MedDRA web-based browser."""


# In[ ]:


import argparse
import pymysql
import pandas as pd
from pandas import ExcelWriter


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument("--mysql_defaults", type=str, help="file to mysql details")
parser.add_argument("--project_dir", help="directory")

args = parser.parse_args()
project_dir = args.project_dir


# In[ ]:


aeolus_conn = pymysql.connect(read_default_file=args.mysql_defaults, unix_socket='/var/run/mysqld/mysqld.sock')
cursor = aeolus_conn.cursor()
my_query = """select distinct upper(pt) from standard_case_outcome"""
cursor.execute(my_query)
result = [i[0] for i in cursor.fetchall()]
aeolus_conn.close()


# In[ ]:


# Need to do MedDRA hierarchy analysis for all these terms
aes_df = pd.DataFrame(result)
aes_df.reset_index(inplace=True)
aes_df.columns = ['Row ID', 'Term']


# In[ ]:


writer = ExcelWriter(project_dir + '/data/all_aeolus_aes_hierarchy_input.xlsx')
aes_df.to_excel(writer,'Sheet1', index=False)
writer.save()
print('Wrote input file with MedDRA terms (total {}) to /data for hierarchy analysis.'.format(len(aes_df)))

