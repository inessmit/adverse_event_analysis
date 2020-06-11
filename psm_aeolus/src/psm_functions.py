"""Functions for preparing dictionaries, preparing and doing the propensity score (PS) matching as well as making various plots of PS distributions on FAERS data from FAERS AEOLUS"""

import matplotlib
matplotlib.use('Agg')
import sqlite3
import random
from sklearn import linear_model
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import pandas as pd
import logging
sns.set_style("whitegrid", {"legend.frameon": True})

import pymysql
import scipy.stats as stats
import statsmodels.sandbox.stats.multicomp as statsmodels
import numpy as np
import collections
import itertools


def load_drug_and_indication_dicts(indications_dict_name, drugs_dict_name, indication_names_file, drug_names_file):
    """Populate two dictionaries with drug and indication ids as keys, respectively, and the objects' names in string as values. Don't return anything.
    indications_dict_name -- str name of dictionary for indication_name
    drugs_dict_name -- str name of dictionary for drug names
    indication_names_file -- str path to file with indication names having concept_id (column 1), concept_name (column 2), expects a header
    drug_names_file -- str path to file with drug names having concept_id (column 1), concept_name (column 2), expects a header"""

    with open(indication_names_file) as f:
        data = f.readlines()[1:]
    for row in data:
        items = row.strip('\n').split('\t')
        indications_dict_name[int(items[0])] = items[1]

    with open(drug_names_file) as f:
        data = f.readlines()[1:]
    for row in data:
        items = row.strip('\n').split('\t')
        drugs_dict_name[int(items[0])] = items[1]

def populate_drugind_2reports(primaryid_file, isr_file, dict_name):
    """Create dictionary with indication or drug as key and dictionary containing the sets of primaryids or isrs
    as value from inputfiles
    kwargs: primaryid_file -- str filepath for file with primaryid (column 1) and attribute id (drug or indication) (column 2), expects a header
            isr_file -- str filepath for file with isr (column 1) and attribute id (drug or indication) (column 2), expects a header
            dict_name -- str name of existing empty dict to be populated"""

    with open(primaryid_file, 'r') as f:
        primaryid_data = f.readlines()[1:]
    with open(isr_file, 'r') as f:
        isr_data = f.readlines()[1:]

    for row in primaryid_data:

        mydata = row.strip().split('\t')
        item0 = int(mydata[0])  # the primaryid or isr
        item1 = int(mydata[1])  # the indication or drug

        if item1 not in dict_name.keys():
            dict_name[item1] = {'primaryids': set(), 'isrs': set()}
            dict_name[item1]['primaryids'].add(item0)

        elif item1 in dict_name.keys():
            dict_name[item1]['primaryids'].add(item0)

    for row in isr_data:

        mydata = row.strip().split('\t')
        item0 = int(mydata[0])  # the primaryid or isr
        item1 = int(mydata[1])  # the indication or drug

        if item1 not in dict_name.keys():
            dict_name[item1] = {'primaryids': set(), 'isrs': set()}
            dict_name[item1]['isrs'].add(item0)

        elif item1 in dict_name.keys():
            dict_name[item1]['isrs'].add(item0)


def populate_reports2attributes(indications_file, drugs_file, dict_name):
    """Populate the report2attributes dictionary, which holds the drugs and indications as attributes for each report id.
    Run for either isr or primaryid files. If an id in either file does not have occurence in other file, it is still added to the set with an empty attribute set
    kwargs: indications_file -- str name of file with report id (column 1) and indication_concept_id (column2), include header in file
            drugs_file -- str name of file with report id (column 1) and drug_concept_id (column 2), include header in file
            dict_name -- str name of dictionary that will be populated"""

    with open(indications_file, 'r') as f:
        indications_data = f.readlines()[1:]
    with open(drugs_file, 'r') as f:
        drugs_data = f.readlines()[1:]

    for row in indications_data:
        mydata = row.strip().split('\t')
        item0 = int(mydata[0])  # the primaryid or isr
        item1 = int(mydata[1])  # the indication

        if item0 not in dict_name.keys():
            dict_name[item0] = {'indications': set(), 'drugs': set()}
            dict_name[item0]['indications'].add(item1)
        else:
            dict_name[item0]['indications'].add(item1)

    for row in drugs_data:
        mydata = row.strip().split('\t')
        item0 = int(mydata[0])  # the primaryid or isr
        item1 = int(mydata[1])  # the drug

        if item0 not in dict_name.keys():
            dict_name[item0] = {'indications': set(), 'drugs': set()}
            dict_name[item0]['drugs'].add(item1)
        else:
            dict_name[item0]['drugs'].add(item1)


def load_all_dictionaries(drug_names_file, indication_names_file, primaryid_indications, isr_indications, primaryid_drugs, isr_drugs):
    """Return the following populated dictionaries: drug_names, indication_names, drug2reports, indication2reports, report2attributes_primaryid, report2attributes_isr
    keyword arguments:
    drug_names_file -- str path to file with drug names having concept_id (column 1), concept_name (column 2), expects a header
    indication_names_file -- str path to file with indication names having concept_id (column 1), concept_name (column 2), expects a header
    primaryid_indications -- str filepath for file with primaryid (column 1) and indication id (column 2), expects a header
    isr_indications -- str filepath for file with isr (column 1) and indication id (column 2), expects a header
    primaryid_drugs -- str filepath for file with primaryid (column 1) and drug id (column 2), expects a header
    isr_drugs -- str filepath for file with isr (column 1) and drug id (column 2), expects a header
    """

    # prepare the dictionaries with drug and indication ids and corresponding names
    drug_names = dict()
    indication_names = dict()
    load_drug_and_indication_dicts(indications_dict_name = indication_names, drugs_dict_name = drug_names
                                   , indication_names_file = indication_names_file, drug_names_file = drug_names_file)

    # prepare the two separate dictionaries with report ids for the indications and drugs, respectively
    drug2reports = dict()
    indication2reports = dict()
    populate_drugind_2reports(primaryid_file=primaryid_indications, isr_file=isr_indications,
                                              dict_name=indication2reports)
    populate_drugind_2reports(primaryid_file=primaryid_drugs, isr_file=isr_drugs,
                                              dict_name=drug2reports)

    # prepare the dictionaries holding the attributes {drugs, indications} for each report {a patient}
    report2attributes_isr = dict()
    report2attributes_primaryid = dict()
    populate_reports2attributes(indications_file=primaryid_indications, drugs_file=primaryid_drugs,
                                                dict_name=report2attributes_primaryid)
    populate_reports2attributes(indications_file=isr_indications, drugs_file=isr_drugs,
                                                dict_name=report2attributes_isr)

    return drug_names, indication_names, drug2reports, indication2reports, report2attributes_primaryid, report2attributes_isr



def retrieve_features(drug_id, inds_corrs_db_path, drugs_corrs_db_path):
    """Return two lists with correlated drugs and correlated indications respectively for given drug.
    drug -- int concept_id of the drug
    inds_corrs_db_path -- str path to the sqlite db containing correlations with indications
    drugs_corrs_db_path -- str path to the sqlite db containing correlations with other drugs"""

    # retrieve previously calculated correlated indications
    conn = sqlite3.connect(inds_corrs_db_path)
    cursor = conn.cursor()
    correlated_indications = [(i[0],i[1]) for i in cursor.execute('select indication_id, fishers_odds from indications where bh_pvalue < 0.05 and fishers_odds > 1 and drug_id = {} order by bh_pvalue limit 200'.format(str(drug_id))).fetchall()]
    conn.close()

    # retrieve previously calculated correlated indications
    conn = sqlite3.connect(drugs_corrs_db_path)
    cursor = conn.cursor()
    correlated_drugs = [(i[0],i[1]) for i in cursor.execute('select drug2_id, fishers_odds from drugs where bh_pvalue < 0.05 and fishers_odds > 1 and drug1_id = {} order by bh_pvalue limit 200'.format(str(drug_id))).fetchall()]

    conn.close()

    ind_df = pd.DataFrame(correlated_indications, columns = ['id', 'fishers_odds'])
    ind_df['type'] = 'indication'
    
    drug_df = pd.DataFrame(correlated_drugs, columns = ['id', 'fishers_odds'])
    drug_df['type'] = 'drug'
    
    all_df = ind_df.append(drug_df)
    all_df.sort_values(by = 'fishers_odds', ascending=False, inplace=True)
    
    top_200 = all_df.head(200)

    top_correlated_indications = list(top_200.loc[top_200['type']=='indication',:]['id'])
    top_correlated_drugs = list(top_200.loc[top_200['type']=='drug',:]['id'])

    return top_correlated_drugs, top_correlated_indications

def define_case_control_reports(drug_id, indication_reports, drug_reports, correlated_indications, correlated_drugs):
    """Assemble case and control reports from dictionaries and return two dictionaries with case report ids and control report ids,
    each of which has a set of both isr ids and primaryids.
    drug_id -- int concept_id of the drug
    indication_reports -- name of dict that holds report ids per indication
    drug_reports -- name of dict that holds report ids per drug
    correlated_indications -- str list of correlated indications for the drug
    correlated_drugs -- str list of correlated other drugs for the drug"""
    case_reports = dict()
    case_reports['primaryids'] = drug_reports[drug_id]['primaryids']
    case_reports['isrs'] = drug_reports[drug_id]['isrs']

    control_reports = {'primaryids': set(), 'isrs': set()}

    for feature in correlated_indications:
        control_reports['primaryids'] = control_reports['primaryids'].union(
            indication_reports[feature]['primaryids'] - case_reports['primaryids'])
        control_reports['isrs'] = control_reports['isrs'].union(
            indication_reports[feature]['isrs'] - case_reports['isrs'])

    for feature in correlated_drugs:
        control_reports['primaryids'] = control_reports['primaryids'].union(
            drug_reports[feature]['primaryids'] - case_reports['primaryids'])
        control_reports['isrs'] = control_reports['isrs'].union(drug_reports[feature]['isrs'] - case_reports['isrs'])
    
    return case_reports, control_reports


def limit_control_nr(control_reports):
    """If there are more than 100,000 control reports, change list of control reports by limiting total number of
    samples to 100,000 by random sampling without replacement (as in Tatonetti's paper) and while doing this
    take a proportional sample of isr and primaryid reports.
    control_reports -- dict of control reports"""
    control_n = len(control_reports['isrs']) + len(control_reports['primaryids'])
    limited_control_reports = dict()
    if control_n < 1e5:
        return control_reports
    
    if control_n > 1e5:
        isrs_proportion = len(control_reports['isrs']) / control_n
        primaryid_proportion = len(control_reports['primaryids']) / control_n

        random.seed(42)
        limited_control_reports['isrs'] = random.sample(control_reports['isrs'], round(1e5 * isrs_proportion))
        limited_control_reports['primaryids'] = random.sample(control_reports['primaryids'], round(1e5 * primaryid_proportion))

        assert len(limited_control_reports['isrs']) + len(limited_control_reports['primaryids']) <= 1e5
        
        return limited_control_reports


def prepare_regression_input(drug, case_reports, control_reports, report2attributes_primaryid, report2attributes_isr, correlated_drugs, correlated_indications):
    """Order reports and retrieve feature vector for each and append to the report_attributes_vector (first argument for logistic regression)
    and also return the label_vector which concatenates '1' for cases and then '0' for non-cases.
    case_reports -- dict of case reports
    control_reports -- dict of control reports
    report2attributes_primaryid -- str dictionary previously defined containing drug and indication attributes per report
    report2attributes_isr -- str name of dictionary previously defined defined containing drug and indication attributes per report
    correlated_indications -- str list of correlated indications for the drug
    correlated_drugs -- str list of correlated other drugs for the drug"""
    case_reports['primaryids'] = sorted(case_reports['primaryids'])
    case_reports['isrs'] = sorted(case_reports['isrs'])
    control_reports['primaryids'] = sorted(control_reports['primaryids'])
    control_reports['isrs'] = sorted(control_reports['isrs'])

    report_attributes_vectors = []
    valid_case_reports = {'isrs': set(), 'primaryids': set()}
    
    for case in case_reports['primaryids']:
        attributes_vector = [1 if i in report2attributes_primaryid[case]['drugs'] else 0 for i in correlated_drugs] + [
            1 if i in report2attributes_primaryid[case]['indications'] else 0 for i in correlated_indications]
        if sum(attributes_vector) > 0:
            report_attributes_vectors.append(attributes_vector)
            valid_case_reports['primaryids'].add(case)

    for case in case_reports['isrs']:
        attributes_vector = [1 if i in report2attributes_isr[case]['drugs'] else 0 for i in correlated_drugs] + [
            1 if i in report2attributes_isr[case]['indications'] else 0 for i in correlated_indications]
        if sum(attributes_vector) > 0:
            report_attributes_vectors.append(attributes_vector)
            valid_case_reports['isrs'].add(case)

    for control in control_reports['primaryids']:
        attributes_vector = [1 if i in report2attributes_primaryid[control]['drugs'] else 0 for i in correlated_drugs] + [
            1 if i in report2attributes_primaryid[control]['indications'] else 0 for i in correlated_indications]
        report_attributes_vectors.append(attributes_vector)

    for control in control_reports['isrs']:
        attributes_vector = [1 if i in report2attributes_isr[control]['drugs'] else 0 for i in correlated_drugs] + [
            1 if i in report2attributes_isr[control]['indications'] else 0 for i in correlated_indications]
        report_attributes_vectors.append(attributes_vector)

    nr_valid_cases = len(valid_case_reports['primaryids']) + len(valid_case_reports['isrs'])
    nr_original_cases = len(case_reports['primaryids']) + len(valid_case_reports['isrs'])
    nr_controls = len(control_reports['isrs']) + len(control_reports['primaryids'])
    
    logging.info('{} nr original cases: {}, nr valid cases: {}'.format(str(drug), str(nr_original_cases), str(nr_valid_cases)))

    label_vector = [1] * nr_valid_cases + [0] * nr_controls

    return (valid_case_reports, control_reports, report_attributes_vectors, label_vector)


def return_predicted_by_model(report_attributes_vectors, label_vector):
    """build model and predict on same report, return list of propensity scores"""

    model = linear_model.LogisticRegression()
    model = model.fit(report_attributes_vectors, label_vector)

    predicted = []
    for j in model.predict_proba(report_attributes_vectors):
        predicted.append(j[1])

    return predicted


def plot_scores(data_cases, data_controls, drug, savefig_path, drug_names):
    '''save image of plot of distribution of propensity scores among cases and controls
    kwargs: data_cases: list of propensity scores in cases
            data_controls: list of propensity scores in controls
            drug: drug_id
            savefig_path: path to dir for saving the images
            drug_names: name of dict holding drug id to name'''

    # plot cases


    plt.hist(data_cases, bins=20)

    axes = plt.gca()
    start, end = axes.get_ylim()

    plt.title('Propensity scores among exposed (n = {}) ({})'.format(str(len(data_cases)), drug_names[drug]))
    plt.ylabel('Number of reports')
    plt.xlabel('Propensity score (binned)')

    drug_name = drug_names[drug].replace('/', '-').replace(' ', '_')
    if len(drug_name) > 200:
        drug_name = ''

    plt.savefig(savefig_path + str(drug) + '_' + drug_name + '_exposed.png', dpi = 199)

    plt.close()

    # plot controls

    plt.hist(data_controls, bins=20)

    axes = plt.gca()
    axes.set_ylim([start, end])  # use same y-lim as in the cases graph

    plt.title('Propensity scores among non-exposed (n = {}) ({})'.format(str(len(data_controls)), drug_names[drug]))
    plt.ylabel('Number of reports (cut at {})'.format(str(end)))
    plt.xlabel('Propensity score (binned)')

    plt.savefig(savefig_path + str(drug) + '_' + drug_name + '_non-exposed.png', dpi = 199)
    plt.close()

def plot_cases_per_bin(drug, cases, controls, predicted, label_vector, data_path, savefig_path, drug_names):
    """Save image of plot of number of exposed cases in each of the PS bins (general bins) and associated data as text file,
    as wel as a text file of bins specific to case PS distribution, on which sampling would actually be based.
    kwargs: drug: current drug id
            cases: dictionary of reports for cases (as created by define_case_control_reports)
            controls: dictionary of reports for controls (as created by define_case_control_reports)
            predicted: vector of predicted propensity scores (order is sorted cases primaryids, sorted case isrs, sorted control primaryids, sorted control isrs, as created by prepare_regression_input )
            label_vector: vector of true labels (as created by prepare_regression_input)
            data_path: path to dir for text file with data
            savefig_path: path to dir for image file
            drug_names: name of dict holding drug id to name"""

    def f(x):
        if len(x) != 0:
            return sum(x['true_label']) / len(x)
        else:
            return 0

    id_vector = sorted(cases['primaryids']) + sorted(cases['isrs']) + sorted(controls['primaryids']) + sorted(
        controls['isrs'])
    type_vector = ['primaryid'] * len(cases['primaryids']) + ['isr'] * len(cases['isrs']) + ['primaryid'] * len(
        controls['primaryids']) + ['isr'] * len(controls['isrs'])

    score_label_pairs = [i for i in zip(predicted, label_vector, id_vector, type_vector)]

    label_df = pd.DataFrame(score_label_pairs, columns=['propensity_score', 'true_label', 'report_id', 'id_type'])

    results = {}
    for name, group in label_df.groupby(pd.cut(label_df['propensity_score'],
                                               bins=[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                                                     0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.1], right=False)):
        results[name] = f(group), sum(group['true_label']), len(group) - sum(group['true_label']), len(group)

    results_df = pd.DataFrame.from_dict(results, orient='index')
    results_df.sort_index(inplace=True)
    results_df.columns = ['fraction_cases', 'nr_cases', 'nr_controls', 'cases+controls']

    # save data as text file
    drug_name = drug_names[drug].replace('/', '-').replace(' ', '_')
    if len(drug_name) > 200:
        drug_name = ''

    results_df.to_csv(path_or_buf=data_path + str(drug) + '_' + drug_name + 'FractionPerBin.txt', sep="\t",
                      index_label='Propensity score bin')

    # plot the results
    ax = results_df['fraction_cases'].plot(kind='bar', legend=False, rot=80, ylim=(0, 1.1), xlim=(0, 1),
                                             title='Fraction of cases ({}) per propensity score bin'.format(
                                                 drug_name))

    for p in ax.patches:
        ax.annotate(str(round(p.get_height(), 2)), (p.get_x() * 1.007, p.get_height() * 1.007))

    ax.set_xlabel('Propensity score bin')
    ax.set_ylabel('Fraction of cases')

    # save image
    # plt.show()
    plt.tight_layout()
    plt.savefig(savefig_path + str(drug) + '_' + drug_name + 'FractionPerBin.png', dpi = 199)
    plt.close()

    # divide both cases and controls in 20 bins based on cases' propensity scores, save results as text file

    cases_df = label_df.loc[label_df['true_label'] == 1, :]
    controls_df = label_df.loc[label_df['true_label'] == 0, :]

    cases_bins = pd.cut(cases_df['propensity_score'], bins=20, include_lowest = True, precision=2, retbins=True)[1]

    grouped_cases = cases_df.groupby(pd.cut(cases_df['propensity_score'], bins=20, precision=2, include_lowest = True))
    grouped_controls = controls_df.groupby(
        pd.cut(controls_df['propensity_score'], bins=cases_bins, precision=2, include_lowest = True))

    summary = []
    for group, group2 in zip(grouped_cases, grouped_controls):
        summary.append([group[0], sum(group[1]['true_label']),
                        len(group2[1])])  # group[0] first element is name of the group, second element is group/df

    summary_df = pd.DataFrame(summary)
    summary_df.columns = ['Propensity score bin', 'nr_cases', 'nr_controls']

    summary_df.to_csv(path_or_buf=data_path + str(drug) + '_' + drug_name + 'FractionPerCasesBin.txt',
                      sep="\t",
                      index=False)


def plot_scores_conditional(data_cases, data_controls, drug, savefig_path, drug_names):
    '''save image of plot of distribution of propensity scores among cases and controls
    kwargs: data_cases: list of propensity scores in cases
            data_controls: list of propensity scores in controls
            drug: drug_id
            savefig_path: path to dir for saving the images
            drug_names: name of dict holding drug id to name'''

    # plot cases


    plt.hist(data_cases, bins=20)

    axes = plt.gca()
    start, end = axes.get_ylim()

    plt.title('Propensity scores among exposed (n = {}) ({})'.format(str(len(data_cases)), drug_names[drug]))
    plt.ylabel('Number of reports')
    plt.xlabel('Propensity score (binned)')

    drug_name = drug_names[drug].replace('/', '-').replace(' ', '_')

    if len(drug_name) > 200:
        drug_name = ''

    plt.savefig(savefig_path + str(drug) + '_' + drug_name + '_exposed.png', dpi=199)
    # plt.show()
    plt.close()

    # plot controls

    n, bins, patches = plt.hist(data_controls, 20)

    if max(n) >= end:
        plt.close()

        f = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
        ax = plt.subplot(gs[0])

        n, bins, patches = ax.hist(data_controls, 20)
        start2, end2 = ax.get_ylim()
        new_start = sorted(n, reverse=True)[2]
        ax.set_ylim([new_start, end2])
        # tick_spacing = (end - (new_start))/3
        # tick_spacing = 40000
        # ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

        # ax.set_ylim([5000,100000])

        ax2 = plt.subplot(gs[1])
        ax2.hist(data_controls, 20)
        ax2.set_ylim([start, end])

        ax.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop='off')  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()

        d = .01  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
        ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d / 3, 1 + d / 3), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d / 3, 1 + d / 3), **kwargs)  # bottom-right diagonal

        my_title = 'Number of reports'
        f.text(0.01, 0.5, my_title, va='center', rotation='vertical')
        plt.xlabel('Propensity score (binned)')
        #
        #
        ax.set_title(
            'Propensity scores among non-exposed (n = {}) ({})'.format(str(len(data_controls)), drug_names[drug]))

        plt.tight_layout()
        #plt.show()
        plt.savefig(savefig_path + str(drug) + '_' + drug_name + '_non-exposed.png', dpi = 199)
        plt.close()

    if max(n) < end:
        # plot controls
        plt.close()

        plt.hist(data_controls, bins=20)

        plt.title('Propensity scores among non-exposed (n = {}) ({})'.format(str(len(data_controls)), drug_names[drug]))
        plt.ylabel('Number of reports')
        plt.xlabel('Propensity score (binned)')
        plt.tight_layout()
        #plt.show()
        plt.savefig(savefig_path + str(drug) + '_' + drug_name + '_non-exposed.png', dpi = 199)
        plt.close()


def do_ps_plots(drug, predicted, attributes_vector, label_vector, cases, controls, drug_names, savefig_path, data_path):

    # define some lists needed for plotting ps functions
    assert len(predicted) == len(attributes_vector) == len(label_vector)

    cases_predicted_scores = predicted[:len(cases['primaryids']) + len(cases['isrs'])]
    controls_predicted_scores = predicted[len(cases['primaryids']) + len(cases['isrs']):]

    # save plots
    plot_scores_conditional(data_cases=cases_predicted_scores, data_controls=controls_predicted_scores, drug=drug,
                                drug_names=drug_names, savefig_path=savefig_path)
    plot_cases_per_bin(drug=drug, cases=cases, controls=controls, predicted=predicted,
                                       label_vector=label_vector, data_path=data_path
                                       , savefig_path=savefig_path, drug_names=drug_names)


def sample_controls_and_cases(cases, controls, predicted, label_vector):
    """Divide cases and controls based on PS scores for cases. Group by PS and sample from controls 10x the number of reports in cases bin. Discard bins for which no controls available.
    Return dictionaries with sampled control and case reports.
    cases -- as prepared by prepare_objects
    controls -- as prepared by prepare_objects
    predicted -- as prepared by prepare_objects
    label_vector -- as prepared by prepare_objects
    """

    # Zip all vectors and make a df with propensity score, true label, report id, report type to dataframe

    id_vector = sorted(cases['primaryids']) + sorted(cases['isrs']) + sorted(controls['primaryids']) + sorted(
        controls['isrs'])
    type_vector = ['primaryid'] * len(cases['primaryids']) + ['isr'] * len(cases['isrs']) + ['primaryid'] * len(
        controls['primaryids']) + ['isr'] * len(controls['isrs'])
    score_label_pairs = [i for i in zip(predicted, label_vector, id_vector, type_vector)]

    label_df = pd.DataFrame(score_label_pairs, columns=['propensity_score', 'true_label', 'report_id', 'id_type'])

    # separate cases and controls as separate dfs
    cases_df = label_df.loc[label_df['true_label'] == 1, :]
    controls_df = label_df.loc[label_df['true_label'] == 0, :]

    # divide cases into 20 equally spaced bins (this depends on the range of scores observed)
    cases_bins = pd.cut(cases_df['propensity_score'], bins=20, right=False, retbins=True)[
        1]  # the first element here is the binning, obtained because of retbins == True

    # now divide controls using the same binning (cases_bins derived from previous line)
    grouped_cases = cases_df.groupby(pd.cut(cases_df['propensity_score'], bins=20, right=False))
    grouped_controls = controls_df.groupby(pd.cut(controls_df['propensity_score'], bins=cases_bins, right=False))

    # the following samples control reports for bins based on cases, with replacement.
    # then adds report ids to the sampled_controls dict. then for those bins for which controls are available, all cases are added to the sampled_cases dict.

    sampled_controls = {'isrs': [], 'primaryids': []}
    sampled_cases = {'isrs': [], 'primaryids': []}

    for name, group in grouped_cases:
        try:
            nr_samples_needed = len(grouped_cases.get_group(name)['report_id'])
        except KeyError:
            continue

        try:
            len(grouped_controls.get_group(name)['report_id'])
        except KeyError:
            continue

        def g(x, dictionary):
            if x[1] == 'isr':
                dictionary['isrs'].append(x[0])
            if x[1] == 'primaryid':
                dictionary['primaryids'].append(x[0])

        grouped_controls.get_group(name).sample(nr_samples_needed * 10, replace=True, random_state=42)[
            ['report_id', 'id_type']].apply(lambda x: g(x, sampled_controls), axis=1)
        grouped_cases.get_group(name)[['report_id', 'id_type']].apply(lambda x: g(x, sampled_cases),
                                                                      axis=1)  # this only includes bins for which both cases and controls were available

    return (sampled_controls, sampled_cases)


def calculate_differences(correlated_list, item, names_dict, top_correlated, attributes_vector, cases, controls,
                          sampled_cases, sampled_controls, report2attributes_isr, report2attributes_primaryid):
    """
    item -- 'drugs' or 'indications'
    names_dict -- drug_names or indication_names dictionary as loaded by load_all_dictionaries
    correlated_list -- either correlated_drugs list or correlated_indications list
    top_correlated -- number of top correlated features to calculate the data for e.g. top 10 correlated indications"""


    differences = []

    for num, attribute in enumerate(correlated_list[:top_correlated]):  

        cases_before_have = 0
        control_before_have = 0
        control_after_have = 0
        cases_after_have = 0

        for report_id in controls['primaryids']:
            if attribute in report2attributes_primaryid[report_id][item]:
                control_before_have += 1

        for report_id in controls['isrs']:
            if attribute in report2attributes_isr[report_id][item]:
                control_before_have += 1

        for report_id in cases['primaryids']:
            if attribute in report2attributes_primaryid[report_id][item]:
                cases_before_have += 1

        for report_id in cases['isrs']:
            if attribute in report2attributes_isr[report_id][item]:
                cases_before_have += 1

        for report_id in sampled_controls['primaryids']:
            if attribute in report2attributes_primaryid[report_id][item]:
                control_after_have += 1

        for report_id in sampled_controls['isrs']:
            if attribute in report2attributes_isr[report_id][item]:
                control_after_have += 1

        for report_id in sampled_cases['primaryids']:
            if attribute in report2attributes_primaryid[report_id][item]:
                cases_after_have += 1

        for report_id in sampled_cases['isrs']:
            if attribute in report2attributes_isr[report_id][item]:
                cases_after_have += 1

        attribute_name = names_dict[attribute]
        percentage_cases_before = (cases_before_have / (len(cases['primaryids']) + len(cases['isrs']))) * 100
        percentage_control_before = (control_before_have / (len(controls['primaryids']) + len(controls['isrs']))) * 100
        percentage_cases_after = (cases_after_have / (len(sampled_cases['primaryids']) + len(sampled_cases['isrs']))) * 100
        percentage_control_after = (control_after_have / (len(sampled_controls['primaryids']) + len(sampled_controls['isrs']))) * 100

        differences.append([attribute_name, percentage_cases_before, percentage_control_before, percentage_cases_after,
                            percentage_control_after])

    return (differences)


def plot_differences(drug, differences_list, item, data_path, savefig_path, drug_names):
    """Take differences in concomitant drugs or indications among pre-post matching samples as calculated by calculate_differences and save plot of data.
    kwargs:
    drug -- drug concept id
    differences_list -- output data from calculate_differences function
    item -- 'drugs' or 'indications'
    path_differences -- path to directory for saving plot images
    drug_names -- str name of dictionary for drug names
    """

    drug_name = drug_names[drug].replace('/', '-').replace(' ', '_')
    if len(drug_name) > 200:
        drug_name = ''

    if item == 'drugs':
        column_name = 'drug name'
        topic = 'concomitant drug use'

    if item == 'indications':
        column_name = 'indication name'
        topic = 'disease indications'

    differences_df = pd.DataFrame(differences_list,
                                  columns=[column_name, 'exposed before', 'non-exposed before', 'exposed after',
                                           'non-exposed after'])
    differences_df.to_csv(path_or_buf = data_path + str(drug) + '_' + drug_name + '_' + item + 'Differences.txt', sep = '\t')

    differences_df['difference before'] = abs(differences_df['exposed before'] - differences_df['non-exposed before'])
    differences_df['difference after'] = abs(differences_df['exposed after'] - differences_df['non-exposed after'])

    differences_df['y value'] = [i for i in range(1, len(differences_df) + 1)]

    x = differences_df['difference before']
    x2 = differences_df['difference after']
    y = differences_df['y value']
    labels = [name.lower() + ' (' + str(nr) +')' for name,nr in zip(differences_df[column_name],differences_df['y value'])]

    ax = plt.subplot(111)
    for num in range(1, len(differences_df) + 1):
        plt.axhline(y=num, color='lightGray')
    plt.axvline(x=0, color='lightGray')
    plt.plot(x, y, 'ro', color='MediumBlue', label='Before PSM', marker='D')
    plt.plot(x2, y, 'ro', color='Orange', label='After PSM')
    plt.yticks(y, labels)
    plt.margins(0.05)

    box = ax.get_position()
    ax.legend(numpoints=1, loc=10, bbox_to_anchor=(1.15, box.y0))
    plt.xlabel('Absolute difference in percentage of reports with {}'.format(topic))

    ttl = ax.title
    ttl.set_position([.5, 1.03])

    #plt.tight_layout()
    plt.savefig(savefig_path + str(drug) + '_' + drug_name + '_{}'.format(item) + '_differences.png', dpi = 199, bbox_inches='tight')
    plt.close()


def make_report2ae_dicts(db_defaults):

    faers_conn = pymysql.connect(read_default_file=db_defaults, unix_socket='/var/run/mysqld/mysqld.sock')
    cursor = faers_conn.cursor()

    # retrieve all AEs for all the isrs ids
    cursor.execute('select distinct isr,outcome_concept_id from standard_case_outcome where isr is not null')
    result = [i for i in cursor.fetchall()]

    # retrieve all AEs for all the primaryids
    cursor.execute('select distinct primaryid,outcome_concept_id from standard_case_outcome where primaryid is not null')
    result2 = [i for i in cursor.fetchall()]
    
    cursor.execute('select distinct outcome_concept_id, pt from standard_case_outcome')
    result3 = [i for i in cursor.fetchall()]

    faers_conn.close()

    report2ae_isrs = {}
    
    if result:
        for item in result:
            ae_id = item[1]
            report_id = int(item[0])

            if report_id not in report2ae_isrs.keys():
                report2ae_isrs[report_id] = set()

            report2ae_isrs[report_id].add(ae_id)           
    
    report2ae_primaryids = {}
    
    if result2:       
        for item in result2:
            ae_id = item[1]
            report_id = int(item[0])

            if report_id not in report2ae_primaryids.keys():
                report2ae_primaryids[report_id] = set()

            report2ae_primaryids[report_id].add(ae_id)
    
    ae_concept_id2pt = {}
    
    for item in result3:
        ae_id = item[0]
        pt = item[1]
        ae_concept_id2pt[ae_id] = pt

    return report2ae_isrs, report2ae_primaryids, ae_concept_id2pt

def combine_results(report_dict_for_cases, report_dict_for_controls, report2ae_isrs, report2ae_primaryids, ae_concept_id2pt, drug):
    """Combine the two dictionaries with AE data into one big list with results (i.e. counts and fractions of AE reported in exposed and non-exposed groups.
    Return list of all AE statistics including PRR and significance (as list of lists containing ['concept id', 'adverse event', 'number of cases affected', 'total cases',
                                           'number of controls affected', 'total controls', 'unique controls',
                                           'fraction of cases affected', 'fraction of controls affected', 'PRR', 'chi-square statistic', 'chi-square p-value')]
    kwargs: sampled cases -- dictionary with report ids for exposed (e.g. from sample_cases_controls function)
            sampled controls -- dictionary with report ids for non-exposed (e.g. from sample_cases_controls function)
            ae_dict_for_all_aes -- dictionary with report ids per adverse event from make_ae_dict_all_reports function
            """
    logging.info('{}: Starting total counts..'.format(str(drug)))
    all_aes = {}
    total_nr_controls = len(report_dict_for_controls['isrs']) + len(report_dict_for_controls['primaryids'])
    total_nr_cases = len(report_dict_for_cases['isrs']) + len(report_dict_for_cases['primaryids'])
    unique_controls = len(set(report_dict_for_controls['isrs'])) + len(set(report_dict_for_controls['primaryids']))
    
    logging.info('{}: Start finding AEs in cases and controls'.format(str(drug)))
    # Find which AEs occur among both controls and cases
    cases_aes_isrs = set(itertools.chain(*[report2ae_isrs[report_id] for report_id in report_dict_for_cases['isrs']]))
    cases_aes_primaryids = set(itertools.chain(*[report2ae_primaryids[report_id] for report_id in report_dict_for_cases['primaryids']]))
    cases_aes = cases_aes_isrs | cases_aes_primaryids
    
    controls_aes_isrs = set(itertools.chain(*[report2ae_isrs[report_id] for report_id in report_dict_for_controls['isrs']]))
    controls_aes_primaryids = set(itertools.chain(*[report2ae_primaryids[report_id] for report_id in report_dict_for_controls['primaryids']]))
    controls_aes = controls_aes_isrs | controls_aes_primaryids
    
    logging.info('{}: Start finding overlapping AEs'.format(str(drug)))
    aes_to_do = cases_aes & controls_aes        
    
    # AEs that occur in both exposed and non-exposed, one of the conditions for the PRR as pointed out here:
    # van Puijenbroek EP, Bate A, Leufkens HG, Lindquist M, Orre R, Egberts AC. A comparison of measures of
    # disproportionality for signal detection in spontaneous reporting systems for adverse drug reactions.
    # Pharmacoepidemiol Drug Saf. 2002 Jan-Feb;11(1):3-10. PubMed PMID: 11998548
    logging.info('{}: Start counting per AE'.format(str(drug)))
    for current_ae in aes_to_do:
        
        nr_cases_affected = 0
        nr_controls_affected = 0
        
        for report_id in report_dict_for_cases['isrs']:
            if current_ae in report2ae_isrs[report_id]:
                nr_cases_affected += 1
        
        for report_id in report_dict_for_cases['primaryids']:
            if current_ae in report2ae_primaryids[report_id]:
                nr_cases_affected += 1
        
        if total_nr_controls != unique_controls:
            # Convenience dict isrs
            isrs2ae = dict()
            for report_id in set(report_dict_for_controls['isrs']):
                if current_ae in report2ae_isrs[report_id]:
                    isrs2ae[report_id] = 1
                else:
                    isrs2ae[report_id] = 0

            # Convenience dict primaryids
            primaryids2ae = dict()
            for report_id in set(report_dict_for_controls['primaryids']):
                if current_ae in report2ae_primaryids[report_id]:
                    primaryids2ae[report_id] = 1
                else:
                    primaryids2ae[report_id] = 0
            # Get counts
            for report_id in report_dict_for_controls['isrs']:
                if isrs2ae[report_id] == 1:
                    nr_controls_affected += 1

            for report_id in report_dict_for_controls['primaryids']:
                if primaryids2ae[report_id] == 1:
                    nr_controls_affected += 1
        
        # Don't use convenience dicts if nr unique controls is not different from total controls
        elif total_nr_controls == unique_controls:
            for report_id in report_dict_for_controls['isrs']:
                if current_ae in report2ae_isrs[report_id]:
                    nr_controls_affected += 1
        
            for report_id in report_dict_for_controls['primaryids']:
                if current_ae in report2ae_primaryids[report_id]:
                    nr_controls_affected += 1

        fraction_cases_affected = nr_cases_affected / total_nr_cases
        
        fraction_controls_affected = nr_controls_affected / total_nr_controls
        
        nr_cases_unaffected = total_nr_cases - nr_cases_affected
        nr_controls_unaffected = total_nr_controls - nr_controls_affected
        
        ae_name = ae_concept_id2pt[current_ae]
        comment = None
        
        #print([nr_cases_affected, nr_cases_unaffected, nr_controls_affected, nr_controls_unaffected])
        # only proceed if observed frequency in each cell is at least 5 (necessary for chi-squared test).
        if all([i > 4 for i in [nr_cases_affected, nr_cases_unaffected, nr_controls_affected, nr_controls_unaffected]]) != True:
            PRR, chi2, chi2pvalue = None, None, None
            comment = 'observed frequencies below 5: [{}]'.format(' '.join(str(i) for i in [nr_cases_affected, nr_cases_unaffected, nr_controls_affected, nr_controls_unaffected]))

        else:
            # calculate the Proportional Reporing Ratio (PRR) according to van Puijenbroek et al. (reference above)
            PRR = (nr_cases_affected / (nr_cases_affected + nr_cases_unaffected)) / (nr_controls_affected / (nr_controls_affected + nr_controls_unaffected))

            # significance test on the PRR, this does chi-square test with Yates correction as is commonly used, also according to van Puijenbroek et al. (reference above)
            obs = np.array([[nr_cases_affected, nr_cases_unaffected], [nr_controls_affected, nr_controls_unaffected]])
            chi2_results = stats.chi2_contingency(obs)
            expected = chi2_results[3]
            if all([i > 4 for i in expected.ravel()]) != True: # do not save chi-squared results if expected frequencies in each cell are not at least 5
                PRR, chi2, chi2pvalue = None, None, None
                comment = 'expected frequencies below 5: {}'.format(' '.join(format(i, '.2f') for i in expected.ravel()))
            else:
                chi2 = chi2_results[0]
                chi2pvalue = chi2_results[1]

        ae_info = collections.OrderedDict([('concept id', current_ae)
                   , ('adverse event', ae_name)
                   , ('number of exposed affected', nr_cases_affected)
                   , ('total exposed', total_nr_cases)
                   , ('number of non-exposed affected', nr_controls_affected)
                   , ('total non-exposed', total_nr_controls)
                   , ('unique non-exposed', unique_controls)
                   , ('fraction of exposed affected', fraction_cases_affected)
                   , ('fraction of non-exposed affected', fraction_controls_affected)
                   , ('PRR', PRR)
                   , ('chi-squared statistic', float(format(chi2, '.2f')) if chi2 else None)
                   , ('chi-squared p-value', chi2pvalue)
                   , ('comment', comment)
                   , ('corrected p-value', None)])
        all_aes[current_ae] = ae_info

    # Make ordered list of the concept ids (AEs) of which PRR is not None
    ordered_ae_list = sorted([key for key in all_aes.keys() if all_aes[key]['PRR'] is not None]) # the key is the concept id of the AE

    # Retrieve PRR values for these concepts (AEs) only. Only these are subject to multiple testing correction
    prr_pvalues = [all_aes[concept_id]['chi-squared p-value'] for concept_id in ordered_ae_list]
    
    if len(prr_pvalues) > 1: # Multiple testing correction is only needed if there are at least two AEs with PRRs (large enough sample size)
        # Do multiple testing correction
        multitest_results = statsmodels.multipletests(prr_pvalues, alpha=0.25, method='fdr_bh')

        for concept_id, corrected_pvalue in zip(ordered_ae_list, multitest_results[1]):
            all_aes[concept_id]['corrected p-value'] = corrected_pvalue

    all_aes_list = list(all_aes.values())

    return all_aes_list


def save_all_ae_results_csv(drug, all_aes, data_path, drug_names, file_ending):
    """Save dataframe with results from combine_results function as csv file in specified directory.
    kwargs: drug -- drug id
            all_aes -- list of dicts with results from combine_results function
            data_path -- path to directory for storing csv files
            drug_names -- dictionary of drug_names based on concept id"""

    # define drug name
    drug_name = drug_names[drug].replace('/', '-').replace(' ', '_')
    if len(drug_name) > 200:
        drug_name = ''

    # make dataframe
    ae_df = pd.DataFrame(all_aes)

    # save dataframe
    ae_df.sort_values(by='PRR', ascending=False, inplace=True)
    ae_df.to_csv(path_or_buf=data_path + str(drug) + '_' + drug_name + file_ending, sep="\t", index=False)
