import pandas as pd
import sklearn.metrics
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import logging
import os

def compute_lr(my_target, bioact_df, ae_dict, my_ae, assume_inactive=False, includes_predictions=False):
    """Compute likelihood ratio (LR) and associated counts for AE-target combination, return.
    kwargs: my_ae - str name of adverse event
            my_target - str uniprot ID of target
            ae_dict - name of ad dictionary, faers or sider
            bioact_df - bioactivity df to use, e.g. faers, reaxys"""
    
    target_df = bioact_df.loc[bioact_df['accession']==my_target,:]
    
    if assume_inactive == False:
        active_molregnos = list(target_df.loc[target_df['integrated_plasma_activity']==1,'parent_molregno'])
        molregnos = list(target_df['parent_molregno'])
        activity_vector = list(target_df['integrated_plasma_activity'])
        ae_vector = [1 if my_ae in ae_dict[molregno] else 0 for molregno in molregnos]
        if includes_predictions == True:
            predicted_vector = list(target_df['predicted'])
    
    elif assume_inactive == True:
        active_molregnos = list(target_df['parent_molregno'])
        molregnos = sorted(ae_dict.keys())
        activity_vector = [1 if molregno in active_molregnos else 0 for molregno in molregnos]
        ae_vector = [1 if my_ae in ae_dict[molregno] else 0 for molregno in molregnos]
        if includes_predictions == True:
            predicted_vector = list(target_df['predicted'])
        
    # if the adverse event is not happening for any of the compounds
    if sum(ae_vector) < 1:
        return
        
    y_pred = activity_vector.copy()
    confusion_result = sklearn.metrics.confusion_matrix(ae_vector, y_pred).ravel()
    
    # If the data size is too small, sklearn doesn't return normal result (normal result has length 4)
    if len(confusion_result) < 4:
        logging.info('{}, {}: Quitting because data size is too small for two-class vectors'.format(my_target, my_ae))
        return
    
    tn, fp, fn, tp = confusion_result
    LR = (tp * (fp + tn) ) / ( fp * (tp + fn))
    
    AE_hit_rate = tp / (tp + fn)
    nAE_hit_rate = fp / (fp + tn)
    
    AE_count = ae_vector.count(1)
    nAE_count = ae_vector.count(0)
    active_count = sum(y_pred)
    inactive_count = len(y_pred) - sum(y_pred)
    
    odds, pvalue = stats.fisher_exact([[tp,fn],[fp,tn]])

    results = [my_target
            ,len(ae_vector)
            , AE_count
            , AE_hit_rate
            , nAE_count
            , nAE_hit_rate
            , active_count
            , inactive_count
            , my_ae
            , LR
            , pvalue
            , activity_vector
            , ae_vector
            , molregnos
            , active_molregnos]
    if includes_predictions == True:
        results.append(predicted_vector)
    return results

def bar_plot(my_target, my_ae, activity_vector, ae_vector, my_output_dir, target_name, organism_string):
    """Do bar plot of hit rates (active/inactive) and save in 'plots' directory within my_output_dir
    kwargs: my_target - str uniprot id
            my_ae - str name of adverse event
            activity_vector - binary vector of bioactivity data
            ae_vector - binary vector of ae data
            my_output_dir - str full path to the directory for saving results, must contain 'plots' directory
            organism_string -- str organism name with spaces replaced by underscore"""
    
    
    # Prepare data
    ae_total = ae_vector.count(1)
    nae_total = ae_vector.count(0)
    
    # Bioactivity to boolean 
    plot_data_df = pd.DataFrame({'Bioactivity': activity_vector, 'Adverse event': ae_vector})
    
    # Count boolean actives
    ae_active = plot_data_df.groupby('Adverse event').get_group(1)['Bioactivity'].sum()
    nae_active = plot_data_df.groupby('Adverse event').get_group(0)['Bioactivity'].sum()
    
    # Fractions of active and inactive compounds for plotting on the bar plot
    ae_fractions = (ae_active/ae_total, nae_active/nae_total)
    nae_fractions = ((ae_total-ae_active)/ae_total, (nae_total-nae_active)/nae_total)
    
    # Prepare plots
    n_groups=2
    fig, ax = plt.subplots()
    index = np.arange(n_groups)
    bar_width = 0.35
    opacity = 0.5

    rects1 = ax.bar(index, ae_fractions, bar_width,
                alpha=opacity, color='r'
                , label='Active')

    rects2 = ax.bar(index + bar_width, nae_fractions, bar_width,
                alpha=opacity, color = 'grey'
                , label='Inactive')

    ax.set_ylabel('Fraction')
    ax.set_xticks(index + bar_width / 2)
    ax.set_xticklabels(('Compounds with AE (n={})'.format(str(ae_total)), 'Compounds without AE (n={})'.format(str(nae_total))))
    ax.legend()
    
    ax.set_title(target_name + ' and ' + my_ae.lower())
    
    fig.tight_layout()
    
    fig.savefig(my_output_dir + '/plots/' + my_target + '_' + organism_string + '_' + my_ae.lower().replace(' ', '_').replace('/', '_'), dpi=200, bbox_inches='tight')
    plt.close()

def find_significant_associations(output_dir, LR_threshold, sign_threshold, min_nr_compounds_active, min_nr_compounds_ae, target_info):
    """For a given experiment/output_dir, find significant positive-direction target-AE associations and save to file.
    Also save a file with Nr. of associated AEs per target for significant targets
    kwargs: output_dir -- str name of output_dir
            LR_threshold -- threshold for the likelihood ratio
            sign_thresold -- threshold for the corrected p-value
            min_nr_compounds_active -- minimum number of compounds active at target
            min_nr_compounds_ae -- minimum number of compounds associated with one adverse event
            target_info -- dataframe with target names to merge on uniprot accession etc. 
            """
    
    significant_associations = pd.DataFrame()
    files = [i for i in os.listdir(output_dir + '/data') if 'ratios.txt' in i]
    
    # Loop over all targets
    for file in files:  
        target_df = pd.read_csv(output_dir + '/data/{}'.format(file), sep='\t')
        target_df = target_df.loc[(target_df['nr compounds with AE']>=min_nr_compounds_ae)&(target_df['nr compounds active']>=min_nr_compounds_active)&(~target_df['Likelihood Ratio'].isnull()),:]
        
        if len(target_df) < 1:
            continue
        
        significant = target_df.loc[(target_df['corrected p-value']<=sign_threshold)&(target_df['Likelihood Ratio']>= LR_threshold),:]
        if len(significant) > 0:
            significant_associations = significant_associations.append(significant, sort=False)    
    
    # Save copy of all significant associations
    file_info = 'sign_associations_LR{}_sign{}_minAE{}_minACT{}'.format(str(LR_threshold), str(sign_threshold), str(min_nr_compounds_ae), str(min_nr_compounds_active))
    significant_associations.sort_values(by='Likelihood Ratio', ascending=False).to_csv(output_dir + '/' + file_info + '.txt', sep='\t', index=None)
    
    # Group significant associations by accession to get count of AE per target
    target_counts = pd.DataFrame(significant_associations.groupby('accession')['Adverse Event'].count())
    target_counts.reset_index(inplace=True, drop=False)
    target_counts.columns = ['accession', 'Nr. of associated AEs']    
    target_counts = target_counts.merge(target_info, on='accession', how='left')
    target_counts = target_counts[['accession', 'target_organism', 'pref_name', 'Nr. of associated AEs']]
    
    # Save copy of significant targets and counts
    file_info2 = 'sign_targets_counts_LR{}_sign{}_minAE{}_minACT{}'.format(str(LR_threshold), str(sign_threshold), str(min_nr_compounds_ae), str(min_nr_compounds_active))
    target_counts.sort_values(by='Nr. of associated AEs', ascending=False).to_csv(output_dir + '/' + file_info2 + '.txt', sep='\t', index=None)

def find_all_associations(output_dir, min_n=None):
    plt.rcdefaults()

    all_associations = pd.DataFrame()
    files = [i for i in os.listdir(output_dir + '/data') if 'ratios_all.txt' in i]

    # Loop over all targets
    for file in files:  
        target_df = pd.read_csv(output_dir + '/data/{}'.format(file), sep='\t')

        if len(target_df) < 1:
            continue

        all_associations = all_associations.append(target_df)
    
    if min_n:
        all_associations_min_n = all_associations.loc[(all_associations['nr compounds with AE']>=min_n)&(all_associations['nr compounds active']>=min_n)&(~all_associations['Likelihood Ratio'].isnull()),:]
        return all_associations_min_n
        
    return all_associations

def find_associations(output_dir, min_n, lr, pv, target_info):

    all_associations = pd.DataFrame()
    files = [i for i in os.listdir(output_dir + '/data') if 'ratios.txt' in i]

    # Loop over all targets
    for file in files:  
        target_df = pd.read_csv(output_dir + '/data/{}'.format(file), sep='\t')
        target_df = target_df.loc[(target_df['nr compounds with AE']>=min_n)&(target_df['nr compounds active']>=min_n)&(~target_df['Likelihood Ratio'].isnull()),:]

        if len(target_df) < 1:
            continue

        all_associations = all_associations.append(target_df)
    
    merged = all_associations.merge(target_info, on='accession', how='left')
    significant = merged.loc[(merged['corrected p-value']<= pv)&(merged['Likelihood Ratio']>= lr),:].sort_values(by='Likelihood Ratio', ascending=False)
    
    return merged, significant

def do_dataset_counts(assoc_df, min_n):
    # Restrict df with min n
    assoc_df = assoc_df.loc[(assoc_df['nr compounds with AE']>=min_n)&(assoc_df['nr compounds active']>=min_n)&(~assoc_df['Likelihood Ratio'].isnull()),:]
    assoc_df['Adverse Event'] = assoc_df['Adverse Event'].apply(lambda x: x.upper())
    
    # Compute some counts on the associations dataframe
    # Number of unique adverse events and targets
    unique_aes = len(assoc_df['Adverse Event'].drop_duplicates())
    unique_targets = len(assoc_df['accession'].drop_duplicates())
        
    nr_associations_tested = len(assoc_df[['accession', 'Adverse Event']].drop_duplicates())
    
    # Number of unique compounds need to be extracted from strings/list
    all_molregnos = set()
    for item in list(assoc_df['molregnos'].drop_duplicates()):
        molregnos = [int(i) for i in item.strip('[]').split(', ')]
        for molregno in molregnos:
            all_molregnos.add(molregno)
            
    # Number of unique active compounds need to be extracted from strings/list
    active_molregnos = set()
    for item in list(assoc_df['active_molregnos'].drop_duplicates()):
        stripped = item.strip('[]')
        if len(stripped) != 0:
            molregnos = [int(i) for i in item.strip('[]').split(', ')]
            for molregno in molregnos:
                active_molregnos.add(molregno)
            
    return {'Unique target-AE pairs': nr_associations_tested
         ,'Unique adverse events': unique_aes
         ,'Unique targets': unique_targets
         , 'Unique compounds': len(all_molregnos)
         , 'Of which unique active compounds':  len(active_molregnos)
           }

def calculate_share_of_predictions(bioact_df):
    assert len(bioact_df) == len(bioact_df[['parent_molregno','accession']].drop_duplicates())
    
    # Calculate % of all compound-target pairs derived from measured data or predictions
    measured_df = bioact_df.loc[bioact_df['predicted']==0,:]
    predicted_df = bioact_df.loc[bioact_df['predicted']==1,:]
    perc_pairs_measured = (len(measured_df) / len(bioact_df)) * 100
    perc_pairs_predicted = (len(predicted_df) / len(bioact_df)) * 100
    
    # % of active pairs 
    active_pairs = bioact_df.loc[bioact_df['integrated_plasma_activity']==1,:]
    if len(active_pairs) > 1:
        perc_active_measured = (len(active_pairs.loc[active_pairs['predicted']==0,:]) / len(active_pairs)) * 100
        perc_active_predicted = (len(active_pairs.loc[active_pairs['predicted']==1,:]) / len(active_pairs)) * 100
    else:
        perc_active_measured = 0
        perc_active_predicted = 0
    
    # % of inactive pairs 
    inactive_pairs = bioact_df.loc[bioact_df['integrated_plasma_activity']==0,:]
    if len(inactive_pairs) > 1:
        perc_inactive_measured = (len(inactive_pairs.loc[inactive_pairs['predicted']==0,:]) / len(inactive_pairs)) * 100
        perc_inactive_predicted = (len(inactive_pairs.loc[inactive_pairs['predicted']==1,:]) / len(inactive_pairs)) * 100
    else:
        perc_inactive_measured = 0
        perc_inactive_predicted = 0
        
    return {'Nr compounds': len(set(bioact_df['parent_molregno']))
            ,'% of compounds active': (len(active_pairs) / len(bioact_df)) * 100
            ,'% of compounds inactive': (len(inactive_pairs) / len(bioact_df)) * 100
            ,'% of compound-target pairs from measured data': perc_pairs_measured
            ,'% of compound-target pairs from predicted data': perc_pairs_predicted
            ,'% of active compound-target pairs from measured data': perc_active_measured
            ,'% of active compound-target pairs from predicted data': perc_active_predicted
            ,'% of inactive compound-target pairs from measured data': perc_inactive_measured
            ,'% of inactive compound-target pairs from predicted data': perc_inactive_predicted
           }

def calculate_share_of_predictions_per_target(bioact_df):
    """Return dataframe with percentages of measured and predicted split by active and inactive data points."""
    target_dict = {}
    all_targets = set(bioact_df['accession'])
    for target in all_targets:
        target_dict[target] = calculate_share_of_predictions(bioact_df.loc[bioact_df['accession']==target,:])
    target_df = pd.DataFrame.from_dict(target_dict, orient='index')
    return target_df