## Analysis of target-adverse event associations

Do various plots and analysis on the data in each of the output directories of 'ae_target_links', including
- distributions of Likelihood Ratios and p-values
- Recall of previously reported associations
- Donut plots of SOCs of AEs
- Try out combinations of targets
- Compute additional performance metrics for single targets
Scripts use info from the prev_reported_safety_associations

Notebook 10_further_analysis is with rainpy environment (detailed in `requirements.txt`) for raincloud plots<sup id="a1">[1](#f1)</sup>, installed according to instructions at from [pog87/PtitPrince](https://github.com/pog87/PtitPrince). In addition, the following packages were installed on top of that:
```
conda install scikit-learn
conda install xlrd
```

### Notebooks

1. Run 1_extract_all_aes_for_hierchy_analysis and do MedDRA hierarchy analysis on result
Gets AEs from FAERS AEOLUS and use MedDRA web-based browser hierarchy analysis 

2. Run 2_target_classification_all_levels_r.sql to get target classification

3. Run 3_run_analyses_refactored.sh
The 3_run_analyses_refactored.ipynb defines a dictionary for each of the folders to run.
Need to manually edit parameters here.
The main directory to be used to generate a colour dictionary with SOC2colour for the plots - then this same dictionary is used for other conditions.
Specify which directories to run in the .sh script.
Especially the ae_donut plots, which will be placed in each of the output directories, have some parameters which need to be manually defined.. basically inspect the plot and adjust parameters in the dictionary in the 3_run_analyses_refactored.ipynb.

4. 4.inspect_data_PPV_annotations
Now starting to prepare the comparison of FAERS and SIDER ae-target links.
First create a directory in 'results' for comparison results. 
This notebook: To annotate a few datapoints on the PPV-hit rate scatterplot to be made, inspect the PPV values and create txt files with the datapoints x and y values and gene abbreviations for the plot. Placed in the comparison results dir and used by next script.
Done interactively for each comparison (sorry.. ).

5. 5_run_analyses_faers_versus_sider.sh
Specify destination directory (should be in 'results') and two directories from ae_target_links/output to compare.

#### To redo the annotations for PPV plots - redo steps 4 and 5.

6. Execute 6_combinations_recall_precision.py
Creates files with target combinations and best combinations in ae_target_links output folder.

7. 7_combination_counts.ipynb
Interactive in notebook. In what % of AE assessed was a combination better in retrieval than a single target?
Distribution of improvement/change in recall and PPV.

8. 8_performance_metrics.py
Additional performance metrics for target-AE associations.

9. 9_target_flow_counts.ipynb
Interactive in notebook. Inspect nr of significant targets, novel targets, total targets assessed etc. 

10. 10_further_analysis.ipynb.  
**Uses the rainpy environment**  
Run notebook to create figures of PPV, LR etc. distributions by SOC and target class. Also plots PPV and LR vs prevalence.

11. 11_cutoff_vs_unbound_safety_known.ipynb
To see effect of plasma concentrations on safety target associations: Compare number of associations in the unbound versus cut-off datasets, also overlap, and quantification of LRs etc.

12. 12_adverse_events_related_to_targets.ipynb

13. Execute 13_withdrawn.sql on ChEMBL_25.
To retrieve drugs that are withdrawn + reason of withdrawal from ChEMBL.

14. 14_SOC_priorities.ipynb
Make tables of target-AE associations within high-priority organ systems (vital organ systems)


##### References
<b id="f1">1</b> Allen M, Poggiali D, Whitaker K, Marshall TR, Kievit RA. Raincloud plots: a multi-platform tool for robust data visualization. Wellcome Open Res. 2019;4:63. Published 2019 Apr 1. doi:[10.12688/wellcomeopenres.15191.1](https://doi.org/10.12688/wellcomeopenres.15191.1) [↩](#a1)

