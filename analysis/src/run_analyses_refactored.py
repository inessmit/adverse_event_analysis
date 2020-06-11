#!/usr/bin/env python
# coding: utf-8

# In[1]:


import analysis_functions
import pandas as pd
import importlib
import os
import pickle
import argparse
import datetime


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument("--colour_dict_dir", type=str, help="output directory which should be used for generating first plotting colour dictionary (the largest dataset)")
parser.add_argument("--colour_dict_filename", type=str, help="Filename to use for pickle of colour dict in /analysis/data")
parser.add_argument("--list_of_dirs", type=str, help="comma-separated list of directories for which to do the plots and figures")          
parser.add_argument("--dir_info_pickle", type=str, help='Filename to use for pickle of thresholds information dictionary in /analysis/data')
args = parser.parse_args()


# In[ ]:


basedir = '/scratch/ias41/ae_code'


# In[ ]:


directory_list = args.list_of_dirs.split(',')

# Make sure the colour dict is run first
if args.colour_dict_dir and args.colour_dict_dir in directory_list:
    directory_list.remove(args.colour_dict_dir)
    directory_list = [args.colour_dict_dir] + directory_list


# In[ ]:


# Donut sizes meaning:
# 0: The number/count of AEs above which should be normal label display (not line to label)
# 1: The number/count of AEs above which should start putting a line with label
# 2: The number/count of AEs above which to display the count in the square of the donut
# 3: Counting the donut squares with different colors, at which square should start the line labels?


# In[ ]:


dirs = dict({
            '20200110_faers_unbound_margin_pred_005_PRR2': {'dir': '20200110_faers_unbound_margin_pred_005_PRR2'
            , 'min_n': 5
            , 'lr': 2
            , 'pv': 0.05
            , 'top_lr': 2
            , 'top_pv': 0.05
            , 'ppv': 0                                
            , 'hr': 0
            , 'donut_sizes': [9,2,4,6]
            , 'donut_white': [7,8,14,15,16]
            , 'colour_dict': args.colour_dict_filename
            , 'precision_display': 0.3                                 
            },
            '20191229_faers_unbound_margin_005_PRR2': {'dir': '20191229_faers_unbound_margin_005_PRR2'
            , 'min_n': 5
            , 'lr': 2
            , 'pv': 0.05
            , 'top_lr': 2
            , 'top_pv': 0.05
            , 'ppv': 0                                
            , 'hr': 0
            , 'donut_sizes': [4,1,1,5]
            , 'donut_white': [7,8,14,15,16]
            , 'colour_dict': args.colour_dict_filename
            , 'precision_display': 0.3                                 
            },
            '20200110_faers_total_margin_pred_005_PRR2': {'dir': '20200110_faers_total_margin_pred_005_PRR2'
            , 'min_n': 5
            , 'lr': 2
            , 'pv': 0.05
            , 'top_lr': 2
            , 'top_pv': 0.05
            , 'ppv': 0                                
            , 'hr': 0
            , 'donut_sizes': [12,3,3,5]
            , 'donut_white': [7,9,12]
            , 'colour_dict': args.colour_dict_filename
            , 'precision_display': 0.3                                 
            },
            '20200110_faers_cutoff6_pred_005_PRR2': {'dir': '20200110_faers_cutoff6_pred_005_PRR2'
            , 'min_n': 5
            , 'lr': 2
            , 'pv': 0.05
            , 'top_lr': 2
            , 'top_pv': 0.05
            , 'ppv': 0                                
            , 'hr': 0
            , 'donut_sizes': [53,18,18,10]
            , 'donut_white': [1,5,13,17]
            , 'colour_dict': args.colour_dict_filename
            , 'precision_display': 0.3                                 
            },
        '20200110_sider_unbound_margin_pred': {'dir': '20200110_sider_unbound_margin_pred'
            , 'min_n': 5
            , 'lr': 2
            , 'pv': 0.05
            , 'top_lr': 2
            , 'top_pv': 0.05
            , 'ppv': 0                                
            , 'hr': 0
            , 'donut_sizes': [19,3,5,7]
            , 'donut_white': [2,3,5,8,9,11]
            , 'colour_dict': args.colour_dict_filename
            , 'precision_display': 0.3                                 
            },
        '20191229_sider_unbound_margin': {'dir': '20191229_sider_unbound_margin'
            , 'min_n': 5
            , 'lr': 2
            , 'pv': 0.05
            , 'top_lr': 2
            , 'top_pv': 0.05
            , 'ppv': 0                                
            , 'hr': 0
            , 'donut_sizes': [12,2,2,4]
            , 'donut_white': [1,3,6,7]
            , 'colour_dict': args.colour_dict_filename
            , 'precision_display': 0.3                                 
            },
            '20200110_sider_total_margin_pred': {'dir': '20200110_sider_total_margin_pred'
            , 'min_n': 5
            , 'lr': 2
            , 'pv': 0.05
            , 'top_lr': 2
            , 'top_pv': 0.05
            , 'ppv': 0                                
            , 'hr': 0
            , 'donut_sizes': [22,6,6,7]
            , 'donut_white': [2,6,9]
            , 'colour_dict': args.colour_dict_filename
            , 'precision_display': 0.3                                 
            },
            '20200110_sider_cutoff6_pred': {'dir': '20200110_sider_cutoff6_pred'
            , 'min_n': 5
            , 'lr': 2
            , 'pv': 0.05
            , 'top_lr': 2
            , 'top_pv': 0.05
            , 'ppv': 0                                
            , 'hr': 0
            , 'donut_sizes': [98,25,25,10]
            , 'donut_white': [3,5,8,9,11]
            , 'colour_dict': args.colour_dict_filename
            , 'precision_display': 0.3                                 
            }
       }
        )


# In[ ]:


with open(basedir + '/analysis/data/dirs_info.pkl', 'wb') as f:
    pickle.dump(dirs, f)


# In[ ]:


for directory in dirs:
    dirs[directory]['top_conditions_string'] = "lr{}_pv{}_ppv{}_hr{}".format(dirs[directory]['lr'], dirs[directory]['pv'], dirs[directory]['ppv'], dirs[directory]['hr'])


# In[6]:


# Files and info needed for the analysis

# Target information
target_info = pd.read_csv(basedir + '/ae_target_links/data/target_names.txt', sep='\t')
target_info = target_info.loc[target_info['accession_organism']=='Homo sapiens',:]

# MedDRA hierchy
meddra_hier = pd.read_excel(basedir + '/analysis/data/all_faers_and_sider_aes_hier_output.xlsx', skiprows=4)
meddra_hier_selection = meddra_hier.loc[meddra_hier['Primary SOC']=='Y',[' Term','HLT','SOC']].drop_duplicates()
meddra_hier_selection['HLT'] = meddra_hier_selection['HLT'].apply(lambda x: x.upper())

# ATC codes
all_atc_codes_loc = basedir + '/faers_sider_comparison/data/atc_all.txt'
small_molecule_atc_codes_loc = basedir + '/faers_sider_comparison/data/atc_small_molecules.txt'

# Target classification
chembl_target_classification = pd.read_csv(basedir + '/analysis/data/target_classification_all_levels_r.txt', sep = '\t')

# Previously reported associations
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


# In[ ]:


for directory in directory_list:
         
    data = dirs[directory]
    print('now starting {}'.format(directory))
    
    # Find associations
    all_associations = analysis_functions.find_all_associations(basedir + '/ae_target_links/output/' + data['dir'])
    assoc_df, sign_df = analysis_functions.find_associations(basedir + '/ae_target_links/output/' + data['dir'], min_n = data['min_n'], target_info=target_info, lr=data['lr'], pv=data['pv'])
    print('found associations')
        
    # Promiscuity plots
    analysis_functions.do_promiscuity_plots(sign_df, output_dir=basedir + '/ae_target_links/output/' + data['dir'])
    print('done promiscuity plots')
    
    # Calculate Positive predictive value
    sign_df['PPV'] = sign_df.apply(lambda x: analysis_functions.calculate_ppv(x), axis=1)
    
    # Find top associations
    top_df = sign_df.loc[(sign_df['Likelihood Ratio']>=data['top_lr'])&(sign_df['corrected p-value']<=data['top_pv'])&(sign_df['PPV']>=data['ppv'])&(sign_df['ae_hit_rate']>=data['hr']),:]
    print('found top df')
    
    # Make colour dict only for specified directory, keeps colours consistent across plots
    if directory == args.colour_dict_dir:
        if not os.path.exists(basedir + '/analysis/data/{}'.format(args.colour_dict_filename)):
            analysis_functions.do_heatmap_dict(assoc_df = sign_df, meddra_hier= meddra_hier_selection, colour_dict_loc = basedir + '/analysis/data/{}'.format(args.colour_dict_filename))
            print('done new colour dict')
        else:
            print('not done new colour dict')
    
    # Recall and positive predictive value
    analysis_functions.do_pr_plot_and_txt(all_associations_df=assoc_df, meddra_hier=meddra_hier_selection, known_merged=known_merged, pv_cutoffs=[0.1,0.05,0.01,0.001], y_lim=0.15, x_lim=data['precision_display'], output_loc=basedir + '/ae_target_links/output/' + data['dir'])
    print('done pr plots')
    # Plot distributions of pvalues and LRs
    analysis_functions.plot_pv_lr_dist(associations_df = assoc_df, meddra_hier=meddra_hier_selection, known_merged=known_merged, output_dir=basedir + '/ae_target_links/output/' + data['dir'])
    print('done pv lr distributions')

    # Heatmaps
    analysis_functions.do_heatmap(assoc_df=top_df, meddra_hier=meddra_hier_selection, colour_dict_loc=basedir + '/analysis/data/{}'.format(data['colour_dict']), output_loc=basedir + '/ae_target_links/output/'+ data['dir'], output_filename_conditions=data['top_conditions_string'], clustering_method='complete')
    analysis_functions.do_soc_donut(assoc_df=top_df, meddra_hier=meddra_hier_selection, colour_dict_loc=basedir + '/analysis/data/{}'.format(data['colour_dict']), output_loc=basedir + '/ae_target_links/output/'+ data['dir'], output_filename_conditions=data['top_conditions_string'], sizes=data['donut_sizes'], white_txt_nrs=data['donut_white'])
    analysis_functions.do_atc_bar_plot(assoc_df=all_associations, all_atc_codes_loc=all_atc_codes_loc, small_molecule_atc_codes_loc=small_molecule_atc_codes_loc, output_loc=basedir + '/ae_target_links/output/'+ data['dir'])
    print('done heatmaps and atc')
    
    # Overall Recall
    recall_info = [analysis_functions.overall_recall(associations_df=assoc_df, lr=thresholds[0], pv=thresholds[1], known_merged=known_merged, meddra_hier=meddra_hier_selection) for thresholds in [(2,0.01), (2,0.05), (1,1.1)]]
        
    current_date = datetime.date.today().strftime("%Y%m%d")
    with open(basedir + '/ae_target_links/output/' + data['dir'] + '/{}_overall_recall.txt'.format(current_date), 'w') as f:
        f.write('\n'.join(recall_info))
    print('done overall recall')

