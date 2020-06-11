#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import datetime
import os
import itertools
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from sklearn import metrics
from sklearn.metrics import precision_recall_curve, average_precision_score, recall_score
from adjustText import adjust_text


# In[ ]:


plt.rcdefaults()
sns.set_style("whitegrid")


# In[ ]:


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
        all_associations = all_associations.loc[(all_associations['nr compounds with AE']>=min_n)&(all_associations['nr compounds active']>=min_n)&(~all_associations['Likelihood Ratio'].isnull()),:]
    
    return all_associations


# In[ ]:


def find_associations(output_dir, min_n, lr, pv, target_info):
    plt.rcdefaults()

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
    significant = merged.loc[(merged['corrected p-value']<= pv)&(merged['Likelihood Ratio']>= lr),:]
    
    return merged, significant


# In[72]:


def do_promiscuity_plots(assoc_df, output_dir):
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 15})

    target_compounds = assoc_df[['pref_name','molregnos','activity_vector']].drop_duplicates()
    target2compounds = dict()

    for item in target_compounds.iterrows():
        target_name = item[1][0]
        molregnos = [int(float(i)) for i in item[1][1].strip('[]').split(', ')]
        activities = [int(float(i)) for i in item[1][2].strip('[]').split(', ')]
        active_compounds = [molregno for molregno, activity in zip(molregnos, activities) if activity==1]
        target_name_ext = target_name + ' ({})'.format(len(set(active_compounds)))
        target2compounds[target_name_ext] = set(active_compounds)
    
    sorted_targets = sorted(target2compounds.keys())
    overlap_data = dict()

    for target_1 in sorted_targets:
        overlap_vector = []
        for target_2 in sorted_targets:
            overlap_perc = len(target2compounds[target_1] & target2compounds[target_2])/len(target2compounds[target_1])
            overlap_vector.append(overlap_perc)
        overlap_data[target_1] = overlap_vector
    overlap_matrix_target1 = pd.DataFrame.from_dict(overlap_data, orient='index', columns=sorted_targets)

    plt.rc('axes', labelsize=20)
    plt.rcParams.update({'font.size': 20})
    
    current_date = datetime.date.today().strftime("%Y%m%d")
    fig = plt.figure(figsize=(16,16))
    g = sns.heatmap(overlap_matrix_target1, square=True,cmap='afmhot_r',xticklabels=True,yticklabels=True,  cbar_kws={'shrink':0.3, 'label':'Fraction of active compounds overlapping'})#annot=True,annot_kws={'size':9},
    g.set_ylabel('Proteins for which fraction is displayed (Nr. active compounds)', size=20) # overlapping fraction is shown
    g.set_xlabel('Proteins (Nr. active compounds)', size=20)
    fig.savefig(output_dir + '/{}_promiscuity.png'.format(current_date), dpi=200, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    


# In[ ]:


def calculate_ppv(row):
    ae_active = int(row['nr compounds with AE'] * row['ae_hit_rate'])
    ppv = ae_active / row['nr compounds active']
    return ppv


# In[1]:


def do_hit_rate_ppv_plot(df1, df1_name, df1_color, df2, df2_name, df2_color, output_dir, annotations=None, target_abbreviations=None, additional_filename='', alpha=0.6):
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 15})

    markers = iter(['x','x'])
    colours = iter([df2_color, df1_color])
    distances = [0.1,0.2,0.3,0.4]

    for df in [df1, df2]:
        df['PPV'] = df.apply(lambda x: calculate_ppv(x), axis=1)
    
    df1['kind'] = df1_name
    df2['kind'] = df2_name
    concat_df = pd.concat([df1, df2], sort=False)
    concat_df.rename(columns={'ae_hit_rate': "Fraction of AE-associated drugs that are active"}, inplace=True)
    
    def multivariateGrid(output_dir, col_x, col_y, col_k, df, k_is_color=False, scatter_alpha=alpha):
        plt.rc('axes', labelsize=14)
        plt.rcParams['xtick.labelsize'] = 14
        plt.rcParams['ytick.labelsize'] = 14

        def colored_scatter(x,y,c=None):
            def scatter(*args, **kwargs):
                args = (x, y)
                if c is not None:
                    kwargs['c'] = c
                kwargs['alpha'] = scatter_alpha
                kwargs['s'] = 7
                kwargs['marker'] = next(markers)
                plt.scatter(*args, **kwargs)
            
            return scatter

        g = sns.JointGrid(
            x=col_x,
            y=col_y,
            data=df
        )

        legends=[]
        for name, df_group in df.groupby(col_k):
            color = next(colours)
            legends.append(name)
            if k_is_color:
                color=name
            g.plot_joint(
                colored_scatter(df_group[col_x],df_group[col_y],color),
            )
            
            sns.distplot(
                df_group[col_x].values,
                ax=g.ax_marg_x,
                color=color,
                kde=False,
                norm_hist=True,
                bins=np.arange(0, 1.05, 0.05)
            )
            sns.distplot(
                df_group[col_y].values,
                ax=g.ax_marg_y,
                color=color,
                vertical=True,
                kde=False,
                norm_hist=True,
                bins=np.arange(0, 1.05, 0.05)
            )
        plt.legend(legends, loc=1, fontsize=12) # bbox_to_anchor=(1.2, 1)
        
        if annotations:
            for annotation in annotations:
                plt.annotate(annotation[0], xy=annotation[1], xytext=annotation[2], arrowprops=dict(arrowstyle="->",connectionstyle="arc3", alpha=0.7, color='darkgrey'), fontsize=8, verticalalignment='bottom')

        if target_abbreviations:
            plt.text(x=0, y=-0.2, s=target_abbreviations, fontsize=10, alpha=0.8)
        
        plt.xlim([0, 1.1])
        current_date = datetime.date.today().strftime("%Y%m%d")
        plt.savefig(output_dir + '/{}_hit_rates_va-PPV{}.png'.format(current_date, additional_filename), dpi=200, bbox_inches='tight')
        
        
    multivariateGrid(output_dir, col_y='Value-added PPV', col_x='Fraction of AE-associated drugs that are active', df = concat_df, col_k='kind')
    plt.clf()
    plt.cla()
    plt.close()


# In[ ]:


def do_heatmap_dict(assoc_df, meddra_hier, colour_dict_loc):
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 15})

    # Make an SOC colour dict that can stay consistent between plots

    colour_dict = dict()
    common_socs = assoc_df.merge(meddra_hier, left_on='Adverse Event', right_on=' Term').groupby('SOC').count().sort_values(by='Adverse Event', ascending=False).index
    colour_list = ['#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'] + ["#9970ab", "#c2a5cf", "#e7d4e8"] + ["#c7eae5", "#80cdc1", "#35978f", "#01665e"] + ["#8c510a","#bf812d", "#dfc27d", "#f6e8c3", "#543005", "#f6e8c3"]

    for soc, colour in zip(common_socs, colour_list[:len(common_socs)]):
        colour_dict[soc] = colour

    # Add rest of SOCS
    rest_socs = set(meddra_hier['SOC'].drop_duplicates()) - set(colour_dict.keys())
    for soc, colour in zip(rest_socs, sns.color_palette("dark",len(rest_socs))):
        colour_dict[soc] = colour
    
    with open(colour_dict_loc, 'wb') as f:
        pickle.dump(colour_dict, f)


# In[ ]:


def do_heatmap(assoc_df, clustering_method, meddra_hier, colour_dict_loc, output_loc, output_filename_conditions):
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 20})
    plt.rc('axes', labelsize=18)
    #sns.set(font_scale=1.5)
    
    with open(colour_dict_loc, 'rb') as f:
        imported_colour_dict = pickle.load(f)
    
    pt2colour = dict()
    for item in meddra_hier[[' Term', 'SOC']].drop_duplicates().iterrows():
        pt2colour[item[1][0]] = imported_colour_dict[item[1][1]]
    
    def find_colour(x):
        try:
            colour = pt2colour[x]
            return colour
        except KeyError:
            return 'w'

    unique_aes = sorted(list(assoc_df['Adverse Event'].drop_duplicates()))
    all_targets = list(assoc_df['pref_name'].drop_duplicates())

    # Make matrix
    my_dict = dict()

    for target in all_targets:
        target_df = assoc_df.loc[assoc_df['pref_name']==target,['Likelihood Ratio','Adverse Event']]
        target_aes = list(target_df['Adverse Event'].drop_duplicates())
        lr_vector = []
        for ae in unique_aes:
            if ae in target_aes:
                my_lr = list(target_df.loc[target_df['Adverse Event']==ae,'Likelihood Ratio'])[0]
                if my_lr == np.inf:
                    my_lr = max(assoc_df.loc[assoc_df['Likelihood Ratio']!=np.inf,'Likelihood Ratio'])
            else:
                my_lr = 1
            lr_vector.append(my_lr)
        my_dict[target] = lr_vector   
    ae_matrix = pd.DataFrame.from_dict(my_dict, orient='index', columns=[i.upper() for i in unique_aes])

    # Do plot
    my_colours = pd.Series(ae_matrix.columns.map(find_colour))
    my_colours.name=('MedDRA System Organ Class')
    my_colours.index=ae_matrix.columns

    plt.figure(figsize=(150, 150))

    g = sns.clustermap(ae_matrix, method=clustering_method, yticklabels=1, metric='euclidean',cmap='gist_heat_r',xticklabels=0, col_colors=my_colours, cbar_kws={'label':'Likelihood Ratio (LR)'}) # 'ticks':ticks

    #g.ax_row_dendrogram.set_visible(False)
    #g.ax_col_dendrogram.set_visible(False)

    g.cax.set_position([0, .2, .03, .45])

    #g.cax.set_yticklabels(ticklabels)

    plt.text(0,min(g.cax.get_yticks()-g.cax.get_yticks()[1]),'Not strongly\nassociated\nor no data', wrap=True, size=14)
    plt.text(0,(max(g.cax.get_yticks())+g.cax.get_yticks()[1])/1.2,'Strongly\nassociated', wrap=True, size=14)
    
    ax = g.ax_heatmap
    ax.set_xlabel('{} adverse events'.format(len(unique_aes)), size=18)
    ax.set_ylabel('{} proteins'.format(len(all_targets)), size=18)
    ax.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 12)
    
    current_date = datetime.date.today().strftime("%Y%m%d")
    plt.savefig(output_loc + '/{}_ae_heatmap_{}.png'.format(current_date, output_filename_conditions), dpi=200, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


# In[2]:


def do_pm_heatmap(assoc_df, meddra_hier, pm_enriched, colour_dict_loc, output_loc, output_filename_conditions, clustering_method='complete'):
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 15})
    
    # Import colour dict
    with open(colour_dict_loc, 'rb') as f:
        imported_colour_dict = pickle.load(f)
    
    pt2colour = dict()
    for item in meddra_hier[[' Term', 'SOC']].drop_duplicates().iterrows():
        pt2colour[item[1][0]] = imported_colour_dict[item[1][1]]
    
    def find_colour(x):
        try:
            colour = pt2colour[x]
            return colour
        except KeyError:
            return 'w'

    # Make matrix
    pm_assoc = assoc_df.loc[assoc_df['Adverse Event'].isin([i.upper() for i in pm_enriched.index]),:]
    pm_aes = sorted(list(pm_assoc['Adverse Event'].drop_duplicates()))
    if len(pm_aes) < 2:
        return
    
    pm_targets = list(pm_assoc['pref_name'].drop_duplicates())    
    
    my_dict = dict()

    for target in pm_targets:
        target_df = pm_assoc.loc[pm_assoc['pref_name']==target,['Likelihood Ratio','Adverse Event']]
        target_aes = list(target_df['Adverse Event'].drop_duplicates())

        lr_vector = []

        for ae in pm_aes:
            if ae in target_aes:
                my_lr = list(target_df.loc[target_df['Adverse Event']==ae,'Likelihood Ratio'])[0]
                if my_lr == np.inf:
                    my_lr = max(pm_assoc.loc[pm_assoc['Likelihood Ratio']!=np.inf,'Likelihood Ratio'])
            else:
                my_lr = 1

            lr_vector.append(my_lr)

        my_dict[target] = lr_vector   
    ae_matrix = pd.DataFrame.from_dict(my_dict, orient='index', columns=[i.lower().capitalize() for i in pm_aes])
    
    if len(my_dict) < 2:
        return
    
    # Do plot
    my_colours = pd.Series([i.upper() for i in ae_matrix.columns]).map(find_colour)
    my_colours.name=('MedDRA System Organ Class')
    my_colours.index=ae_matrix.columns

    plt.figure(figsize=(150, 150))
    sns.set(font_scale=1.5)
    g = sns.clustermap(ae_matrix, method=clustering_method, col_colors=my_colours, yticklabels=1, cmap='gist_yarg', metric='euclidean',xticklabels=True, cbar_kws={'label':'Likelihood Ratio (LR)'})#'ticks':[1, 4,8,12,16,20]

    g.cax.set_position([-.10, .2, .03, .45])

    plt.text(0.1,min(g.cax.get_yticks()-g.cax.get_yticks()[1]),'Not strongly\nassociated\nor no data', wrap=True)
    plt.text(0.1,(max(g.cax.get_yticks())+g.cax.get_yticks()[1])/1.2,'Strongly\nassociated', wrap=True)

    current_socs = list(pm_assoc.merge(meddra_hier, left_on='Adverse Event', right_on=' Term', suffixes=['x','']).groupby('SOC').count().sort_values(by='Adverse Event', ascending=False).index)
    #current_socs.append('Not classified')

    for label in current_socs:
         g.ax_col_dendrogram.bar(0, 0, color=imported_colour_dict[label],label=label, linewidth=0)
    g.ax_col_dendrogram.legend(loc='right', ncol=2, bbox_to_anchor=(1.3, 1.8), title='MedDRA System Organ Class',framealpha=0)

    ax = g.ax_heatmap
    ax.set_xlabel('Post-marketing relevant adverse events', size=20)
    ax.set_ylabel('Proteins', size=20)
    plt.rc('axes', labelsize=12)

    current_date = datetime.date.today().strftime("%Y%m%d")
    plt.savefig(output_loc + '/{}_ae_heatmap_pm{}.png'.format(current_date, output_filename_conditions), dpi=200, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


# In[1]:


def do_soc_donut(assoc_df, meddra_hier, colour_dict_loc, output_loc, output_filename_conditions, sizes=[1,1,1,1], white_txt_nrs=[]):
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 15})
    
    with open(colour_dict_loc, 'rb') as f:
        imported_colour_dict = pickle.load(f)

    plt.rcParams['font.size'] = 12
    my_data = assoc_df.merge(meddra_hier, how='left',left_on='Adverse Event', right_on=' Term')
    
    data_values = list(my_data.groupby('SOC').count().sort_values(by='Adverse Event', ascending=False)['accession'].values)
    data_index = list(my_data.groupby('SOC').count().sort_values(by='Adverse Event', ascending=False)['accession'].index)

    data_index_large_only = [i if j > sizes[0]  else '' for i,j in zip(data_index, data_values)]
    data_index_count = [True if j < sizes[1] else False for j in data_values]
    data_index_medium = [i if j > sizes[1]  else '' for i,j in zip(data_index, data_values)]

    def make_autopct(values):
        def my_autopct(pct):
            total = sum(values)
            val = int(round(pct*total/100.0))
            return '{v:d}'.format(v=val)
        return my_autopct
    def absolute_value(val):
        a  = int(np.round(val/100.*sum(data_values), 0))
        if a > sizes[2]:
            return a
        else:
            return ''

    fig, ax = plt.subplots()
    wedges, texts, autotexts = ax.pie(data_values, wedgeprops=dict(width=0.3,edgecolor='white'), pctdistance=0.85, startangle=0, labels=data_index_large_only, colors=[imported_colour_dict[x] for x in data_index], autopct=absolute_value)
    for nr, autotext in enumerate(autotexts):
        if nr in white_txt_nrs:
            autotext.set_color('white')

    bbox_props = dict(boxstyle="square,pad=0", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),zorder=0, va="center") # bbox=bbox_props, 

    for i, p in enumerate(wedges[sizes[3]:len([i for i in data_index_medium if i!=''])], start=sizes[3]):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle, "color": 'black'})
        ax.annotate(data_index_medium[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.6*y),
                    horizontalalignment=horizontalalignment, **kw)

    plt.text(-0.5,0,'      {}\nprotein-AE pairs'.format(len(assoc_df)))

    current_date = datetime.date.today().strftime("%Y%m%d")
    plt.savefig(output_loc + '/{}_ae_soc_donut_{}.png'.format(current_date, output_filename_conditions), dpi=200, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


# In[ ]:


def do_atc_bar_plot(assoc_df, all_atc_codes_loc, small_molecule_atc_codes_loc, output_loc):
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 15})
    plt.figure(figsize=(8, 6)) 
    sns.set_style("white")
    plt.rc('axes', labelsize=15)
    
    all_molregnos = set([int(x) for i in list(assoc_df['molregnos'].drop_duplicates()) for x in i.strip('[]').split(', ')])

    all_atc = pd.read_csv(all_atc_codes_loc, sep='\t')
    project_atc = all_atc.loc[all_atc['molregno'].isin(all_molregnos),:]
    background_atc = pd.read_csv(small_molecule_atc_codes_loc, sep='\t')

    project_len = len(project_atc)
    background_len = len(background_atc)

    background_atc_counts = background_atc.groupby('level1_description').count().sort_values(by='molregno',ascending=False)[['molregno']]
    background_atc_counts.index = [i.lower().capitalize() for i in background_atc_counts.index]
    background_atc_counts.reset_index(inplace=True)

    project_atc_counts = project_atc.groupby('level1_description').count().sort_values(by='molregno',ascending=False)[['molregno']]
    project_atc_counts.index = [i.lower().capitalize() for i in project_atc_counts.index]
    project_atc_counts.reset_index(inplace=True)

    comparison = background_atc_counts.merge(project_atc_counts, how='outer', on='index')
    comparison.columns = ['ATC level 1 description', 'Small molecule approved drugs count', 'Drugs in dataset analysed count']
    comparison.set_index('ATC level 1 description', drop=True,inplace=True)

    comparison['Small molecule approved drugs %'] = comparison['Small molecule approved drugs count'].apply(lambda x: x/sum(comparison['Small molecule approved drugs count'])*100)
    comparison['Dataset drugs %'] = comparison['Drugs in dataset analysed count'].apply(lambda x: x/sum(comparison['Drugs in dataset analysed count'])*100)

    comparison_r = comparison.sort_values(by='Small molecule approved drugs %')
    ind = np.arange(len(comparison_r))
    
    fig, ax = plt.subplots()
    width = 0.35

    ax.barh(ind, comparison_r['Small molecule approved drugs %'], width, label='Approved small molecule drugs (ChEMBL) ({} ATC codes)'.format(background_len), color='darkgrey')
    ax.barh(ind + width, comparison_r['Dataset drugs %'], width, label='Dataset drugs ({} ATC codes)'.format(project_len), color='#67a9cf')

    ax.set(yticks=ind + width, yticklabels=comparison_r.index, ylim=[2*width - 1, len(comparison_r)])
    ax.set_xlabel('Percentage of ATC labels')

    handles, labels = plt.gca().get_legend_handles_labels()
    order = [1,0]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc=(0,1.05))

    #ax.set_xticks([i for i in range(1,21,2)])

    current_date = datetime.date.today().strftime("%Y%m%d")
    plt.savefig(output_loc + '/{}_ATC_bar.png'.format(current_date) , dpi=200, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


# In[ ]:


def do_atc_bar_plot_three_datasets(assoc_df1, df1_name, df1_color, assoc_df2, df2_name, df2_color, all_atc_codes_loc, small_molecule_atc_codes_loc, output_loc):
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 15})
    plt.figure(figsize=(8, 6)) 
    sns.set_style("white")
    plt.rc('axes', labelsize=15)
    
    df1_molregnos = set([int(x) for i in list(assoc_df1['molregnos'].drop_duplicates()) for x in i.strip('[]').split(', ')])
    df2_molregnos = set([int(x) for i in list(assoc_df2['molregnos'].drop_duplicates()) for x in i.strip('[]').split(', ')])

    all_atc = pd.read_csv(all_atc_codes_loc, sep='\t')
    df1_atc = all_atc.loc[all_atc['molregno'].isin(df1_molregnos),:]
    df2_atc = all_atc.loc[all_atc['molregno'].isin(df2_molregnos),:]
    background_atc = pd.read_csv(small_molecule_atc_codes_loc, sep='\t')

    df1_len = len(df1_atc)
    df2_len = len(df2_atc)
    background_len = len(background_atc)

    background_atc_counts = background_atc.groupby('level1_description').count().sort_values(by='molregno',ascending=False)[['molregno']]
    background_atc_counts.index = [i.lower().capitalize() for i in background_atc_counts.index]
    background_atc_counts.reset_index(inplace=True)

    df1_atc_counts = df1_atc.groupby('level1_description').count().sort_values(by='molregno',ascending=False)[['molregno']]
    df1_atc_counts.index = [i.lower().capitalize() for i in df1_atc_counts.index]
    df1_atc_counts.reset_index(inplace=True)
    
    df2_atc_counts = df2_atc.groupby('level1_description').count().sort_values(by='molregno',ascending=False)[['molregno']]
    df2_atc_counts.index = [i.lower().capitalize() for i in df2_atc_counts.index]
    df2_atc_counts.reset_index(inplace=True)

    comparison = background_atc_counts.merge(df1_atc_counts, how='outer', on='index').merge(df2_atc_counts, how='outer', on='index')
    comparison.columns = ['ATC level 1 description', 'Small molecule approved drugs count', f'Drugs in {df1_name} dataset count', f'Drugs in {df2_name} dataset count']
    comparison.set_index('ATC level 1 description', drop=True,inplace=True)

    comparison['Small molecule approved drugs %'] = comparison['Small molecule approved drugs count'].apply(lambda x: x/sum(comparison['Small molecule approved drugs count'])*100)
    comparison[f'{df1_name} drugs %'] = comparison[f'Drugs in {df1_name} dataset count'].apply(lambda x: x/sum(comparison[f'Drugs in {df1_name} dataset count'])*100)
    comparison[f'{df2_name} drugs %'] = comparison[f'Drugs in {df2_name} dataset count'].apply(lambda x: x/sum(comparison[f'Drugs in {df2_name} dataset count'])*100)

    comparison_r = comparison.sort_values(by='Small molecule approved drugs %')
    ind = np.arange(len(comparison_r))
    
    fig, ax = plt.subplots()
    width = 0.3

    ax.barh(ind, comparison_r['Small molecule approved drugs %'], width, label='Approved small molecule drugs (ChEMBL) ({} ATC codes)'.format(background_len), color='darkgrey')
    ax.barh(ind + width, comparison_r[f'{df1_name} drugs %'], width, label='{} drugs ({} ATC codes)'.format(df1_name, df1_len), color=df1_color)
    ax.barh(ind + width + width, comparison_r[f'{df2_name} drugs %'], width, label='{} drugs ({} ATC codes)'.format(df2_name, df2_len), color=df2_color)
    
    ax.set(yticks=ind + width, yticklabels=comparison_r.index, ylim=[2*width - 1, len(comparison_r)])
    ax.set_xlabel('Percentage of ATC labels')

    handles, labels = plt.gca().get_legend_handles_labels()
    order = [2,1,0]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc=(0,1.05))

    #ax.set_xticks([i for i in range(1,21,2)])

    current_date = datetime.date.today().strftime("%Y%m%d")
    plt.savefig(output_loc + '/{}_ATC_bar_combined.png'.format(current_date) , dpi=200, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


# In[ ]:


def do_target_class_bar_plot_without_safety(df1, df1_name, df1_color, df2, df2_name, df2_color, output_loc, chembl_target_classification):
    """chembl_target_classification -- df with chembl target classification, must contain columns accession, level_1, level_2
    """
    
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 16})
    
    sns.set_style("white")
    sns.set_style('ticks')
    plt.rc('axes', labelsize=16)
    
    def find_integrated(x):
        if x['level_2'] == 'Not available':
            return x['level_1']
        else:
            return x['level_2']
    
    df1_init = pd.DataFrame(chembl_target_classification.loc[chembl_target_classification['accession'].isin(set(list(df1['accession']))),:])
    df2_init = pd.DataFrame(chembl_target_classification.loc[chembl_target_classification['accession'].isin(set(list(df2['accession']))),:])
    
    for df in [df1_init, df2_init]:
        df['integrated_level'] = df.apply(find_integrated, axis=1)
    
    df1_counts = df1_init.groupby(['integrated_level'])['accession'].count().reset_index()
    df2_counts = df2_init.groupby(['integrated_level'])['accession'].count().reset_index()

    for df, name in zip([df1_counts, df2_counts], [df1_name, df2_name]):
        df.columns = ['Target Class', name]
        df[name + ' %'] = df[name].apply(lambda x: x/sum(df[name])*100)
           
    comparison = df1_counts.merge(df2_counts, how='outer')
    comparison.fillna(0, inplace=True)
    comparison_r = comparison.sort_values(by=df2_name + ' %')
    comparison_r.set_index('Target Class', inplace=True)

    ind = np.arange(len(comparison_r))

    fig, ax = plt.subplots(figsize=(8,11))
    width = 0.35
    ax.barh(ind, comparison_r[df1_name + ' %'], width, label=df1_name, color=df1_color)
    ax.barh(ind + width, comparison_r[df2_name + ' %'], width, label=df2_name, color=df2_color)
    
    ax.set(yticks=ind + width/2, yticklabels=comparison_r.index, ylim=[2*width - 1, len(comparison_r)])
    ax.set_xlabel('Unique targets in dataset (%)')

    handles, labels = plt.gca().get_legend_handles_labels()
    order = [1,0]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='lower right')
    
    current_date = datetime.date.today().strftime("%Y%m%d")
    plt.savefig(output_loc + '/{}_target_class_bar_without_safety.png'.format(current_date) , dpi=200, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


# In[1]:


def do_target_class_bar_plot(df1, df1_name, df1_color, df2, df2_name, df2_color, safety_targets, df3_name, df3_color, output_loc, chembl_target_classification):
    """chembl_target_classification -- df with chembl target classification, must contain columns accession, level_1, level_2
    safety_targets -- set of known safety targets
    df3_name -- str, name of safety target set"""
    
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 16})
    
    sns.set_style("white")
    sns.set_style('ticks')
    plt.rc('axes', labelsize=16)
    
    def find_integrated(x):
        if x['level_2'] == 'Not available':
            return x['level_1']
        else:
            return x['level_2']
    
    df1_init = pd.DataFrame(chembl_target_classification.loc[chembl_target_classification['accession'].isin(set(list(df1['accession']))),:])
    df2_init = pd.DataFrame(chembl_target_classification.loc[chembl_target_classification['accession'].isin(set(list(df2['accession']))),:])
    df3_init = pd.DataFrame(chembl_target_classification.loc[chembl_target_classification['accession'].isin(set(safety_targets)),:])
    
    for df in [df1_init, df2_init, df3_init]:
        df['integrated_level'] = df.apply(find_integrated, axis=1)
    
    df1_counts = df1_init.groupby(['integrated_level'])['accession'].count().reset_index()
    df2_counts = df2_init.groupby(['integrated_level'])['accession'].count().reset_index()
    df3_counts = df3_init.groupby(['integrated_level'])['accession'].count().reset_index()

    for df, name in zip([df1_counts, df2_counts, df3_counts], [df1_name, df2_name, df3_name]):
        df.columns = ['Target Class', name]
        df[name + ' %'] = df[name].apply(lambda x: x/sum(df[name])*100)
           
    comparison = df1_counts.merge(df2_counts, how='outer').merge(df3_counts, how='outer')
    comparison.fillna(0, inplace=True)
    comparison_r = comparison.sort_values(by=df2_name + ' %')
    comparison_r.set_index('Target Class', inplace=True)

    ind = np.arange(len(comparison_r))

    fig, ax = plt.subplots(figsize=(8,11))
    width = 0.3
    ax.barh(ind, comparison_r[df1_name + ' %'], width, label=df1_name, color=df1_color)
    ax.barh(ind + width, comparison_r[df2_name + ' %'], width, label=df2_name, color=df2_color)
    ax.barh(ind + width + width, comparison_r[df3_name + ' %'], width, label=df3_name, color=df3_color)
    
    ax.set(yticks=ind + width, yticklabels=comparison_r.index, ylim=[2*width - 1, len(comparison_r)])
    ax.set_xlabel('Unique targets in dataset (%)')

    handles, labels = plt.gca().get_legend_handles_labels()
    order = [2,1,0]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
    
    current_date = datetime.date.today().strftime("%Y%m%d")
    plt.savefig(output_loc + '/{}_target_class_bar.png'.format(current_date) , dpi=200, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


# In[2]:


def cumulative_curve(associations_df, known_merged, sort_by_columns, ascending=False):
    
    associations_df_upper = associations_df.copy()
    associations_df_upper['Adverse Event'] = associations_df_upper['Adverse Event'].apply(lambda x: x.upper())
    
    known_tuples = set([(x[1]['Accession'],x[1]['HLT'].upper()) for x in known_merged.iterrows()])
    
    def find_known(row):
        if ((row['accession'],row['HLT'])) in known_tuples:
            return 1
        else:
            return 0

    associations_df_upper['known'] = associations_df_upper.apply(find_known, axis=1)
    associations_df_upper_known = associations_df_upper.loc[associations_df_upper['known']==1,:]
    
    known_retrieved_total_tuples = set([(x[1]['accession'],x[1]['HLT'].upper()) for x in associations_df_upper_known.iterrows()])
                
    associations_df_sorted = associations_df_upper.sort_values(by=sort_by_columns, ascending=ascending)
    
    associations_unique = associations_df_sorted.drop_duplicates(subset=['accession', 'Adverse Event'], keep='first')
    nr_associations_tested = len(associations_unique)
    
    frac_considered = []
    frac_recalled = []

    for threshold in np.linspace(start=0.0, stop=1.0, num=100):
        frac_considered.append(threshold)

        current_df = associations_unique[:int(threshold*len(associations_unique))]
        retrieved_tuples = set([(i,j) for i, j in zip(current_df['accession'], current_df['HLT'])])

        nr_found = len(retrieved_tuples & known_retrieved_total_tuples)
        frac_recalled.append(nr_found/len(known_retrieved_total_tuples))
    
    nr_targets = len(set(known_merged['Accession']) & set(associations_df['accession']))
    return frac_considered, frac_recalled, nr_targets, nr_associations_tested


# In[ ]:


def do_cumulative_recall(df1, df1_name, df1_color, df2, df2_name, df2_color, known_merged, meddra_hier, sort_by_columns, ascending, output_loc):
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 16})
    plt.figure(figsize=(8, 6)) 
    
    sns.set_style("white")
    plt.rc('axes', labelsize=15)

    df1_merged = df1.merge(meddra_hier, left_on='Adverse Event', right_on=' Term')
    df2_merged = df2.merge(meddra_hier, left_on='Adverse Event', right_on=' Term')
    
    df1_considered, df1_recalled, df1_targets_considered, df1_nr_associations_tested = cumulative_curve(df1_merged, sort_by_columns=sort_by_columns, ascending=ascending, known_merged=known_merged)
    df2_considered, df2_recalled, df2_targets_considered, df2_nr_associations_tested = cumulative_curve(df2_merged, sort_by_columns=sort_by_columns, ascending=ascending, known_merged=known_merged)
    
    df1_auc = metrics.auc(df1_considered, df1_recalled)
    df2_auc = metrics.auc(df2_considered, df2_recalled)
    
    x1,y1=[0,1],[0,max(df1_recalled)]
    x2,y2=[0,1],[0,max(df2_recalled)]
    plt.plot(x1,y1,color='black', alpha=0.5)
    plt.plot(x2,y2,color='black', alpha=0.5)
    
    plt.step(x=df1_considered,y=df1_recalled,alpha=1, where='post',label='{} (AUC = {:.3f})'.format(df1_name, df1_auc), color=df1_color)
    plt.step(x=df2_considered,y=df2_recalled,alpha=1, where='post',linestyle=':',label='{} (AUC = {:.3f})'.format(df2_name, df2_auc), color=df2_color)
    
    max_value = max(df1_recalled + df2_recalled)
    plt.axis([0,1,0,max_value+0.05])
    plt.legend()
    plt.ylabel('Fraction recalled')
    plt.xlabel('Fraction of associations considered')
    plt.text(1.1,max_value/2,'{}:\n{} targets overlap with extracted safety targets\nTotal number of target-AE(HLT) pairs = {}\nTotal number of associations = {}'.format(df1_name, df1_targets_considered, len(df1_merged[['accession','HLT']].drop_duplicates()), df1_nr_associations_tested), size=12)
    
    plt.text(1.1,max_value/4,'{}:\n{} targets overlap with exctracted safety targets\nTotal number of target-AE(HLT) pairs = {}\nTotal number of associations = {}'.format(df2_name, df2_targets_considered, len(df2_merged[['accession','HLT']].drop_duplicates()), df2_nr_associations_tested), size=12)
    current_date = datetime.date.today().strftime("%Y%m%d")
    plt.savefig(output_loc + '/{}_cumulative_recall_hlt_normalised.png'.format(current_date), dpi=200, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


# In[ ]:


def pr_recall_at_pv_cutoff(associations_df, pv_cutoff, known_values, all_targets, all_target_ae_combinations):
    associations_df_upper = associations_df.copy()
    associations_df_upper['Adverse Event'] = associations_df_upper['Adverse Event'].apply(lambda x: x.upper())
    
    associations_unique = associations_df_upper.drop_duplicates(subset=['accession', 'Adverse Event'], keep='first')

    significant = associations_unique.loc[associations_unique['corrected p-value']<=pv_cutoff,:]

    accession2lrs = dict()
    for accession in all_targets:
        accession2lrs[accession] = dict()
    for item in significant.iterrows():
        current_accession = item[1]['accession']
        current_ae = item[1]['HLT']
        current_lr = item[1]['Likelihood Ratio']
        accession2lrs[current_accession][current_ae] = current_lr

    # Make accession 2 ae vector
    ae_values = []
    for combination in all_target_ae_combinations:
        accession = combination[0]
        ae = combination[1]
        try:
            lr = accession2lrs[accession][ae]
            if lr == np.inf:
                ae_values.append(100)
            else:
                ae_values.append(lr)
        except KeyError:
            ae_values.append(0)
            
    precision, recall, thresholds = precision_recall_curve(known_values, ae_values)
    
    return precision, recall, thresholds


# In[ ]:


def overall_recall(associations_df, lr, pv, known_merged, meddra_hier):
    """Calculate recall at specified LR and q-value thresholds. Restrict to targets and AEs present in dataset
    , but do not normalise to 1 based on pairs present."""
    
    associations_df_upper = associations_df.copy()
    associations_df_upper['Adverse Event'] = associations_df_upper['Adverse Event'].apply(lambda x: x.upper())
    associations_df_sorted = associations_df_upper.sort_values(by=['corrected p-value', 'Likelihood Ratio'], ascending=[True, False])
    associations_unique = associations_df_sorted.drop_duplicates(subset=['accession', 'Adverse Event'], keep='first')
    
    associations_merged = associations_unique.merge(meddra_hier, left_on='Adverse Event', right_on=' Term')
    
    occurring_targets = set(list(associations_merged['accession']))
    occurring_aes = set(list(associations_merged['HLT']))
    
    # Restrict known associations to the targets and AEs occurring in associations df
    known_restricted = known_merged.loc[(known_merged['Accession'].isin(occurring_targets))&(known_merged['HLT'].isin(occurring_aes)),:]
    known_tuples = set([(x[1]['Accession'],x[1]['HLT'].upper()) for x in known_restricted.iterrows()])
    
    associations_cutoff = associations_merged.loc[(associations_merged['Likelihood Ratio']>=lr)&(associations_merged['corrected p-value']<=pv),:]
    retrieved_tuples = set([(i,j) for i, j in zip(associations_cutoff['accession'], associations_cutoff['HLT'])])

    nr_found = len(retrieved_tuples & known_tuples)
    recall = nr_found/len(known_tuples)

    info = f"""Recall of known associations (at HLT level) restricted to those targets and AEs (HLT) that occur at least once in the dataset.\nAt thresholds LR>={lr} and q-value<={pv}, recall is {recall}, which is {nr_found} out of {len(known_tuples)} known associations (AE-HLT pairs). Number of unique target-HLT associations considered at this threshold: {len(retrieved_tuples)}\n"""
    
    return info


# In[ ]:


def do_pr_plot_and_txt(all_associations_df, meddra_hier, known_merged, pv_cutoffs, y_lim, x_lim, output_loc):
    
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 16})
    plt.figure(figsize=(8, 6)) 
    
    all_associations_df['Adverse Event'] = all_associations_df['Adverse Event'].apply(lambda x: x.upper())
    all_associations_df = all_associations_df.drop_duplicates(subset=['accession', 'Adverse Event'], keep='first')
    all_associations_df = all_associations_df.merge(meddra_hier, left_on='Adverse Event', right_on=' Term')
    
    # Make set of unique AEs to loop over later
    all_aes = sorted(list(all_associations_df['HLT'].drop_duplicates()))
    all_targets = sorted(list(all_associations_df['accession'].drop_duplicates()))

    # Restrict known associations to targets covered in the dataset
    current_targets = list(set(all_associations_df['accession']))
    current_known_associations = known_merged.loc[(known_merged['Accession'].isin(current_targets))&(known_merged['HLT'].isin(all_aes)),:]
    
    all_target_ae_combinations = sorted(list(itertools.product(sorted(all_targets), sorted(all_aes))))
    
    known_tuples = set([(x[1]['Accession'],x[1]['HLT'].upper()) for x in current_known_associations.iterrows()])
    known_values = [1 if combination in known_tuples else 0 for combination in all_target_ae_combinations]

    colours = iter(['r', 'b', 'darkorange', 'teal'])
    
    cutoff_dfs = []
    
    for cutoff in pv_cutoffs:
        precision, recall, thresholds = pr_recall_at_pv_cutoff(associations_df=all_associations_df, pv_cutoff=cutoff, known_values=known_values, all_targets=all_targets, all_target_ae_combinations=all_target_ae_combinations)

        plt.step(recall, precision, alpha=0.5, color = next(colours), where='post',label='q-value <= {}'.format(cutoff))
        
        selected_thresholds = pd.DataFrame({'Precision': precision[:-1], 'Recall': recall[:-1], 'LR Threshold': thresholds})
        selected_thresholds['q-value'] = cutoff
        cutoff_dfs.append(selected_thresholds[['LR Threshold', 'q-value', 'Precision', 'Recall']].drop_duplicates())
        
    
    plt.legend()
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, y_lim])
    plt.xlim([0.0, x_lim])
    
    current_date = datetime.date.today().strftime("%Y%m%d")
    plt.savefig(output_loc + '/{}_PR_cutoffs.png'.format(current_date), bbox_inches='tight', dpi = 199)
    plt.clf()
    plt.cla()
    plt.close()
    
    concatenated_cutoff_dfs = pd.concat(cutoff_dfs).sort_values(by='Recall', ascending=False)
    concatenated_cutoff_dfs.sort_values(by=['LR Threshold', 'q-value'], ascending=[True, True], inplace=True)
    concatenated_cutoff_dfs['LR Threshold'] = concatenated_cutoff_dfs['LR Threshold'].apply(lambda x: '{:.1f}'.format(x))
    concatenated_cutoff_dfs['Precision'] = concatenated_cutoff_dfs['Precision'].apply(lambda x: '{:.2f}'.format(x))
    concatenated_cutoff_dfs['Recall'] = concatenated_cutoff_dfs['Recall'].apply(lambda x: '{:.2f}'.format(x))
    concatenated_cutoff_dfs.drop_duplicates(inplace=True)
    
    def find_nr_assoc(x):
        current_lr = x['LR Threshold']
        current_pv = x['q-value']
        nr_assoc = len(all_associations_df.loc[(all_associations_df['Likelihood Ratio']>=float(current_lr))&(all_associations_df['corrected p-value']<=float(current_pv)),['accession', 'Adverse Event']].drop_duplicates())
        return nr_assoc

    concatenated_cutoff_dfs['Number of associations at this threshold'] = concatenated_cutoff_dfs.apply(find_nr_assoc, axis=1)
    concatenated_cutoff_dfs.to_csv(output_loc + '/{}_PR_cutoffs.txt'.format(current_date), index=False, sep='\t')
    
    


# In[4]:


def plot_pv_lr_dist(associations_df, meddra_hier, known_merged, output_dir):
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 16})
    plt.rcParams["patch.force_edgecolor"] = False
    plt.figure(figsize=(8, 6)) 
    
    sns.set_style("white")
    
    known_tuples = set([(x[1]['Accession'],x[1]['HLT'].upper()) for x in known_merged.iterrows()])
    
    def find_known(row):
        if ((row['accession'],row['HLT'])) in known_tuples:
            return 1
        else:
            return 0
    
    associations_df_upper = associations_df.copy()
    associations_df_upper['Adverse Event'] = associations_df_upper['Adverse Event'].apply(lambda x: x.upper())
    associations_df_sorted = associations_df_upper.sort_values(by=['corrected p-value', 'Likelihood Ratio'], ascending=[True, False])
    associations_unique = associations_df_sorted.drop_duplicates(subset=['accession', 'Adverse Event'], keep='first')
    
    associations_unique = associations_unique.merge(meddra_hier, left_on='Adverse Event', right_on=' Term')
    
    associations_unique['known'] = associations_unique.apply(find_known, axis=1)
    associations_known = associations_unique.loc[associations_unique['known']==1,:]
    associations_rest = associations_unique.loc[associations_unique['known']==0,:]
    
    current_date = datetime.date.today().strftime("%Y%m%d")
    # LRs
    known_lr = associations_known.loc[associations_known['Likelihood Ratio']<5,'Likelihood Ratio']
    rest_lr = associations_rest.loc[associations_rest['Likelihood Ratio']<5,'Likelihood Ratio']
    
    ax = sns.distplot(known_lr, label='Previously reported associations (n={})'.format(str(len(known_lr))), kde=False, norm_hist=True, color='mediumpurple', hist_kws=dict(edgecolor="mediumpurple", linewidth=0), bins=50)
    sns.distplot(rest_lr, label='Other pairs (n={})'.format(str(len(rest_lr))), kde=False, norm_hist=True, color='grey', hist_kws=dict(edgecolor="grey", linewidth=0), bins=50)

    ax.set(xlabel='Likelihood Ratio (up to 5)')
    ax.set(ylabel='Density')

    plt.legend(prop={'size': 12}, loc=1) # bbox_to_anchor=(0.5, -0.15)
    ax.tick_params(labelsize=12)
    plt.xlim(1,5)

    plt.savefig(output_dir + '/{}_LR_dist_3pub-hlt.png'.format(current_date), bbox_inches='tight', dpi = 199)
    plt.clf()
    plt.cla()
    plt.close()

    # P-values full range
    known_pv = associations_known['corrected p-value']
    rest_pv = associations_rest['corrected p-value']

    ax = sns.distplot(known_pv, label='Previously reported associations (n={})'.format(str(len(known_pv))), kde=False, norm_hist=True, color='mediumpurple', bins=50, hist_kws=dict(edgecolor="mediumpurple", linewidth=0))
    sns.distplot(rest_pv, label='Other pairs (n={})'.format(str(len(rest_pv))), kde=False, norm_hist=True, color='grey', hist_kws=dict(edgecolor="grey", linewidth=0), bins=50)

    ax.set(xlabel='q-value')
    ax.set(ylabel='Density')

    plt.legend(prop={'size': 12}, loc=2) #bbox_to_anchor=(0.5, -0.15)
    ax.tick_params(labelsize=12)
    plt.xlim(0,1.01)
    
    
    plt.savefig(output_dir + '/{}_pvalue_dist_3pub-hlt-fullt.png'.format(current_date), bbox_inches='tight', dpi = 199)
    plt.clf()
    plt.cla()
    plt.close()
    
    # P-values low range
    known_pv = associations_known.loc[associations_known['corrected p-value']<=0.3,'corrected p-value']
    rest_pv = associations_rest.loc[associations_rest['corrected p-value']<=0.3,'corrected p-value']
    
    fig, ax = plt.subplots(figsize=(5,8))
    ax = sns.distplot(known_pv, label='Previously reported associations (n={})'.format(str(len(known_pv))), kde=False, norm_hist=True, color='mediumpurple', hist_kws=dict(edgecolor="mediumpurple", linewidth=0), bins=15)
    sns.distplot(rest_pv, label='Other pairs (n={})'.format(str(len(rest_pv))), kde=False, norm_hist=True, color='grey', hist_kws=dict(edgecolor="grey", linewidth=0), bins=15)

    ax.set(xlabel='q-value')
    ax.set(ylabel='Density')

    plt.legend(prop={'size': 12}, loc=1) # bbox_to_anchor=(0.5, -0.15)
    ax.tick_params(labelsize=12)
    plt.xlim(0,0.3)

    plt.savefig(output_dir + '/{}_pvalue_dist_3pub-hlt.png'.format(current_date), bbox_inches='tight', dpi = 199)
    plt.clf()
    plt.cla()
    plt.close()


# In[ ]:


def find_integrated(x):
    if x['level_2'] == 'Not available':
        return x['level_1']
    else:
        return x['level_2']


# In[ ]:


def make_sign_target_overview_table(df1, df1_name, df2, df2_name, chembl_target_classification, known_merged, output_loc):
    known_targets = set(list(known_merged['Accession'].drop_duplicates()))
    
    all_sign = pd.concat([df1, df2], sort=False)

    both_targets = set(list(df1['accession'])) & set(list(df2['accession']))
    df1_only_targets = set(list(df1['accession'])) - set(list(df2['accession']))
    df2_only_targets = set(list(df2['accession'])) - set(list(df1['accession']))
    
    both_sign = all_sign.loc[all_sign['accession'].isin(both_targets)]
    both_df = both_sign.merge(chembl_target_classification, left_on='accession', right_on='accession')
    both_df['dataset'] = 'Both'

    df1_only = df1.loc[df1['accession'].isin(df1_only_targets),['accession','pref_name']].drop_duplicates().merge(chembl_target_classification, left_on='accession', right_on='accession')
    df1_only['dataset'] = df1_name
    df2_only = df2.loc[df2['accession'].isin(df2_only_targets),['accession','pref_name']].drop_duplicates().merge(chembl_target_classification, left_on='accession', right_on='accession')
    df2_only['dataset'] = 'SIDER'
    
    sign_target_overview = pd.concat([df1_only, both_df, df2_only], sort=False)
    sign_target_overview['integrated_level'] = sign_target_overview.apply(find_integrated, axis=1)
    
    sign_target_selected = sign_target_overview[['accession', 'pref_name', 'integrated_level', 'dataset']]
    sign_target_selected.columns = ['accession', 'Target name', 'Target Class', 'dataset']
    sign_target_selected.drop_duplicates(inplace=True)
            
    sign_target_selected['Previously reported safety target'] = sign_target_selected['accession'].apply(lambda x: 1 if x in known_targets else 0)
    
    sorter = ['FAERS', 'Both', 'SIDER']
    sorter_dict = dict(zip(sorter,range(len(sorter))))
    sign_target_selected['rank'] = sign_target_selected['dataset'].map(sorter_dict)
    
    sign_target_selected.sort_values(by=['Previously reported safety target', 'rank', 'Target Class'], inplace=True)
    sign_target_selected.drop(labels='rank', axis=1, inplace=True)
    
    current_date = datetime.date.today().strftime("%Y%m%d")
    sign_target_selected.to_csv(output_loc + '/{}_sign_target_overview.txt'.format(current_date), index=False, sep='\t')

    


# In[ ]:


def plot_actives_per_target(assoc_df1, dataset1_name, df1_color, assoc_df2, dataset2_name, df2_color, output_loc):
    plt.rcdefaults()
    plt.rc('axes', labelsize=16)
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    plt.rcParams["patch.force_edgecolor"] = False
    
    active_counts1 = assoc_df1[['accession','nr compounds active']].drop_duplicates()
    active_counts1.rename(columns = {'nr compounds active': 'Number of actives per protein'}, inplace=True)

    active_counts2 = assoc_df2[['accession','nr compounds active']].drop_duplicates()
    active_counts2.rename(columns = {'nr compounds active': 'Number of actives per protein'}, inplace=True)
    
    fig = plt.figure()
    plt.hist(active_counts1['Number of actives per protein'].values, bins=30, alpha=0.6, label=dataset1_name + f' ({len(active_counts1)} proteins)', color=df1_color)
    plt.hist(active_counts2['Number of actives per protein'].values, bins=30, alpha=0.6, label=dataset2_name + f' ({len(active_counts2)} proteins)', color=df2_color)

    plt.ylabel('Frequency')
    plt.xlabel('Number of active drugs per protein')
    plt.legend(loc=1, fontsize=14)
    current_date = datetime.date.today().strftime("%Y%m%d")
    plt.savefig(f'{output_loc}/{current_date}_pos_actives_hist.png', dpi=200, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()

