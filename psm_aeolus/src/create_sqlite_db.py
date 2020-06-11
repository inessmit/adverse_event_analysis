"""Create and populate an sqlite database with two tables for storing drug-indicationa and drug-drug correlations calculated from raw data (drug-isr, drug-primaryid, indication-isr, indication-primaryid pairs)"""

import scipy.stats as stats
import sqlite3
import statsmodels.sandbox.stats.multicomp as statsmodels
import os
import logging
import itertools
import argparse
import time
import datetime

# Command line argparse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="directory")
args = parser.parse_args()

basedir = args.basedir

# Function for loading text files with raw data (drug-report_id and indication_report_id pairs) and populating a dictionary with this data
def populate_dict_from_file(primaryid_file, isr_file, dict_name):
    """create dictionary with indication or drug as key and dictionary containing the sets of primaryids or isrs as value
    kwargs: primaryid_file (assumed to have a header)
            isr_file (assumed to have a header)
            dict_name"""

    with open(primaryid_file, 'r') as f_pi:
        primaryid_data = f_pi.readlines()[1:]
    with open(isr_file, 'r') as f_isr:
        isr_data = f_isr.readlines()[1:]

    for row_pi in primaryid_data:

        my_data = row_pi.strip().split('\t')
        item0 = int(my_data[0])  # the primaryid or isr
        item1 = int(my_data[1])  # the indication or drug

        if item1 not in dict_name.keys():
            dict_name[item1] = {'primaryids': set(), 'isrs': set()}
            dict_name[item1]['primaryids'].add(item0)

        elif item1 in dict_name.keys():
            dict_name[item1]['primaryids'].add(item0)

    for row_isr in isr_data:

        my_data = row_isr.strip().split('\t')
        item0 = int(my_data[0])  # the primaryid or isr
        item1 = int(my_data[1])  # the indication or drug

        if item1 not in dict_name.keys():
            dict_name[item1] = {'primaryids': set(), 'isrs': set()}
            dict_name[item1]['isrs'].add(item0)

        elif item1 in dict_name.keys():
            dict_name[item1]['isrs'].add(item0)

# Create dictionaries to be populated
drug_reports = dict()
indication_reports = dict()

# Populate dictionaries from text files
populate_dict_from_file(primaryid_file= basedir + '/data/primaryid_indications.txt', isr_file=basedir + '/data/isr_indications.txt', dict_name=indication_reports)
populate_dict_from_file(primaryid_file= basedir + '/data/primaryid_drugs.txt', isr_file=basedir + '/data/isr_drugs.txt', dict_name=drug_reports)

print(datetime.datetime.now().strftime("%H:%M:%S") + ' populated dictionaries')

# Make a lookup dictionary for drug and indication names
drug_names = dict()
indication_names = dict()

with open(basedir + '/data/indication_names.txt') as f_ind:
    data = f_ind.readlines()[1:]
for row in data:
    items = row.strip('\n').split('\t')
    indication_names[int(items[0])] = items[1]

with open(basedir + '/data/drug_names.txt') as f_drug:
    data = f_drug.readlines()[1:]
for row in data:
    items = row.strip('\n').split('\t')
    drug_names[int(items[0])] = items[1]

print(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") + ' made lookup dict')
print(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") + ' now starting to compute correlations, may take over 24h...')

# Define function for calculating drug-indication correlations
def indications_calculate_abcd(my_drug, my_indication, total_reports):
    """calculate a, b, c, and d for two-by-to table.
    kwargs: my_drug (int) for standard_concept_id
            my_indication (int) for standard_concept_id"""
    drug_name = drug_names[my_drug]
    indication_name = indication_names[my_indication]

    # Count reports in both drug and indication tables
    # note that some reports are only in the drug table, no indications known. at the moment I don't think these are excluded. Possibly should exclude?
    #  e.g. if a patient is taking some drugs, but no indications are known for this person, this would currently be in count c. However, it's missing data.

    a = len(drug_reports[my_drug]['primaryids'] & indication_reports[my_indication]['primaryids']) + len(
        drug_reports[my_drug]['isrs'] & indication_reports[my_indication]['isrs'])
    b = len(indication_reports[my_indication]['primaryids'] - drug_reports[my_drug]['primaryids']) + len(
        indication_reports[my_indication]['isrs'] - drug_reports[my_drug]['isrs'])
    c = len(drug_reports[my_drug]['primaryids'] - indication_reports[my_indication]['primaryids']) + len(
        drug_reports[my_drug]['isrs'] - indication_reports[my_indication]['isrs'])
    d = total_reports - len(drug_reports[my_drug]['primaryids'] | indication_reports[my_indication]['primaryids']) - len(
        drug_reports[my_drug]['isrs'] | indication_reports[my_indication]['isrs'])

    assert a + b == len(indication_reports[my_indication]['primaryids']) + len(indication_reports[my_indication]['isrs'])
    assert a + c == len(drug_reports[my_drug]['primaryids']) + len(drug_reports[my_drug]['isrs'])

    e = stats.fisher_exact([[a, b], [c, d]])

    return ({'drug': my_drug, 'drug_name': drug_name, 'indication': my_indication
        , 'indication_name': indication_name, 'tbt_a': a, 'tbt_b': b, 'tbt_c': c, 'tbt_d': d, 'fishers_odds': e[0],
             'fishers_pvalue': e[1]})

# Define function for calculating drug-drug correlations
def drugs_calculate_abcd(my_drug1, my_drug2, total_reports):
    """calculate a, b, c, and d for two-by-to table.
    kwargs: drug1 (int) for standard_concept_id of drug1
            drug2 (int) for standard_concept_id of drug2"""
    drug1_name = drug_names[my_drug1]
    drug2_name = drug_names[my_drug2]

    # reports in both drug and indication tables
    # note that some reports are only in the drug table, no indications known. at the moment I don't think these are excluded. Possibly should exclude?
    #  e.g. if a patient is taking some drugs, but no indications are known for this person, this would currently be in count c. However, it's missing data.

    a = len(drug_reports[my_drug1]['primaryids'] & drug_reports[my_drug2]['primaryids']) + len(
        drug_reports[my_drug1]['isrs'] & drug_reports[my_drug2]['isrs'])
    b = len(drug_reports[my_drug2]['primaryids'] - drug_reports[my_drug1]['primaryids']) + len(
        drug_reports[my_drug2]['isrs'] - drug_reports[my_drug1]['isrs'])
    c = len(drug_reports[my_drug1]['primaryids'] - drug_reports[my_drug2]['primaryids']) + len(
        drug_reports[my_drug1]['isrs'] - drug_reports[my_drug2]['isrs'])
    d = total_reports - len(drug_reports[my_drug1]['primaryids'] | drug_reports[my_drug2]['primaryids']) - len(
        drug_reports[my_drug1]['isrs'] | drug_reports[my_drug2]['isrs'])

    assert a + b == len(drug_reports[my_drug2]['primaryids']) + len(drug_reports[my_drug2]['isrs'])
    assert a + c == len(drug_reports[my_drug1]['primaryids']) + len(drug_reports[my_drug1]['isrs'])

    e = stats.fisher_exact([[a, b], [c, d]])

    return ({'drug1': my_drug1, 'drug1_name': drug1_name, 'drug2': my_drug2
        , 'drug2_name': drug2_name, 'tbt_a': a, 'tbt_b': b, 'tbt_c': c, 'tbt_d': d, 'fishers_odds': e[0],
             'fishers_pvalue': e[1]})

# Calculate drug-indications correlations
results_indications = []

# Calculate number of total reports from lists of isr-ids and primary-ids
# what follows flattens the list of report-id lists
ind_common_reports_primaryid = set(itertools.chain(*[indication_reports[key]['primaryids'] for key in indication_reports.keys()])) & set(itertools.chain(*[drug_reports[key]['primaryids'] for key in drug_reports.keys()]))
ind_common_reports_isr = set(itertools.chain(*[indication_reports[key]['isrs'] for key in indication_reports.keys()])) & set(itertools.chain(*[drug_reports[key]['isrs'] for key in drug_reports.keys()]))
ind_total_reports = len(ind_common_reports_primaryid) + len(ind_common_reports_isr)

# Calculate abcd numbers of 2-by-2 table and Fisher's exact test for each drug-indication pair
for indication in list(indication_reports.keys()):
    for drug in list(drug_reports.keys()):
        results_indications.append(indications_calculate_abcd(drug, indication, ind_total_reports))

print(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") + ' calculated drug-indication correlations')

# Calculate drug-drug correlations
results_drugs = []

# Calculate number of total reports from lists of isr-ids and primary-ids
# what follows flattens the list of report-id lists
drug_common_reports_primaryid = set(itertools.chain(*[drug_reports[key]['primaryids'] for key in drug_reports.keys()]))
drug_common_reports_isr = set(itertools.chain(*[drug_reports[key]['isrs'] for key in drug_reports.keys()]))
drug_total_reports = len(drug_common_reports_primaryid) + len(drug_common_reports_isr)

# Calculate abcd numbers of 2-by-2 table and Fisher's exact test for each drug-drug pair
for drug1 in list(drug_names.keys()):
    for drug2 in list(drug_names.keys()):
        if drug1 == drug2:
            continue
        results_drugs.append(drugs_calculate_abcd(drug1, drug2, drug_total_reports))

print(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") + ' calculated drug-drug correlations')

# Do Benjamini-Hochberg multiple testing correction for drug-indication correlations
indications_fishers_pvalues = [i['fishers_pvalue'] for i in results_indications]
bh_results = statsmodels.multipletests(indications_fishers_pvalues, alpha=0.25, method='fdr_bh')

for result, bh_bool, bh_pvalue in zip(results_indications, bh_results[0], bh_results[1]):
    result['bh_bool'] = bh_bool
    result['bh_pvalue'] = bh_pvalue

# Do Benjamini-Hochberg multiple testing correction for drug-drug correlations
drugs_fishers_pvalues = [i['fishers_pvalue'] for i in results_drugs]
bh_results = statsmodels.multipletests(drugs_fishers_pvalues, alpha=0.25, method='fdr_bh')

for result, bh_bool, bh_pvalue in zip(results_drugs, bh_results[0], bh_results[1]):
    result['bh_bool'] = bh_bool
    result['bh_pvalue'] = bh_pvalue

print(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") + ' did multiple testing correction')

# Remove any previous version of the sqlite database if it exists:
if os.path.exists(basedir + '/data/correlations.db'):
    os.remove(basedir + '/data/correlations.db')

conn = sqlite3.connect(basedir + '/data/correlations.db')
cursor = conn.cursor()

# Create database tables
sql_indications = '''create table indications 
(drug_id int
, drug_name text
, indication_id integer
, indication_name text
, tbt_a integer
, tbt_b integer
, tbt_c integer
, tbt_d integer
, fishers_odds real
, fishers_pvalue real
, bh_bool text
, bh_pvalue real)'''

sql_drugs = '''create table drugs
(drug1_id int
, drug1_name text
, drug2_id integer
, drug2_name text
, tbt_a integer
, tbt_b integer
, tbt_c integer
, tbt_d integer
, fishers_odds real
, fishers_pvalue real
, bh_bool text
, bh_pvalue real)'''

cursor.execute(sql_indications)
conn.commit()
cursor.execute(sql_drugs)
conn.commit()

print(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") + ' created db tables')

# Inserts results in sqlite db table for drug-indications correlations
for item in results_indications:
    sql = '''insert into indications (drug_id, drug_name, indication_id, indication_name, tbt_a, tbt_b, tbt_c, tbt_d, fishers_odds, fishers_pvalue, bh_bool, bh_pvalue) values (?,?,?,?,?,?,?,?,?,?,?,?)'''
    cursor.execute(sql, (item['drug'], item['drug_name'], item['indication'], item['indication_name'], item['tbt_a'], item['tbt_b'], item['tbt_c']
                         , item['tbt_d'], item['fishers_odds'], item['fishers_pvalue'], str(item['bh_bool']), item['bh_pvalue']))

# Inserts results in sqlite db table for drug-drug correlations
for item in results_drugs:
    sql = '''insert into drugs 
    (drug1_id, drug1_name, drug2_id, drug2_name, tbt_a, tbt_b, tbt_c, tbt_d, fishers_odds
    , fishers_pvalue, bh_bool, bh_pvalue) values (?,?,?,?,?,?,?,?,?,?,?,?)'''
    cursor.execute(sql, (item['drug1'], item['drug1_name'], item['drug2']
                         , item['drug2_name'], item['tbt_a'], item['tbt_b'], item['tbt_c']
                         , item['tbt_d'], item['fishers_odds'], item['fishers_pvalue']
                         , str(item['bh_bool']), item['bh_pvalue']))

print(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") + ' inserted results')

conn.commit()
conn.close()