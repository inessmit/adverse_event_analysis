## Analysis of Adverse Events and experimental and predicted protein target interactions

This repository contains the code and notebooks for my PhD project analysing statistical associations between adverse drug effects and protein targets of drugs. Adverse events are extracted from the Food and Drug Administration Adverse Event Reporting System (FAERS)<sup id="a1">[1](#f1)</sup><sup id="a2">[,2](#f2)</sup> and Side Effect Resource (SIDER)<sup id="a3">[3](#f3)</sup><sup id="a4">[,4](#f4)</sup>.


### Directory structure
Each step of the analysis has its own directory, each containing its own README, data, src, and results directories, or sometimes a set of notebooks.

Generally the 'data' is too large to be included and must be downloaded from original sources or is created using SQL queries. Smaller data files are included. In other cases the input data is the result from another directory, so the directories depend on each other and the particular order they were executed, so the 'basedir' is specificied at the top of scripts/notebooks.

### Conda environments
The environment as specified in `requirements.txt` was used unless otherwise specified. Other environments were used occasionally as specified within README of the individual directory or at the start of notebooks.

The environment in requirements.txt was created with:
```
conda create -n release python=3.6 anaconda statsmodels=0.9.0 pymysql=0.8.1 nb_conda_kernels
pip install multiprocessing-logging==0.2.6 adjusttext=0.7.3.1
conda install matplotlib-venn
```

### DIRECTORIES/SUB-PARTS IN ORDER:

[compound mapping](compound_mapping): Maps drug concepts from FAERS AEOLUS to ChEMBL identifiers/molecular structures via RxNorm, UniChem, and synonyms.
```
compound_mapping
├── data
├── results
└── src
```

[psm_aeolus](psm_aeolus): Applies Propensity Score Matching (PSM) as reported in Tatonetti et al.<sup id="a5">[5](#f5)</sup> on FAERS AEOLUS<sup id="a2">[2](#f2)</sup> database.
```
psm_aeolus
├── data
├── logs
├── results
│   ├── data
│   └── figures
└── src
```

[faers_aes](faers_aes): Extracts the adverse events with proportional reporting ratio (PRR) => 2 or 3 from the results files in the [psm_aeolus](psm_aeolus) and saves them into a single pickle.
```
faers_aes
├── results
└── src
```

[indications](indications): Downloads an indepdendent dataset of drug indications from the RxClass API [6]
```
indications
├── data
├── results
└── src
```

[psm_effect_inds](psm_effect_inds): Analyses the effect of PSM on indication bias using the independent set of drug indications from the RxClass API<sup id="a6">[6](#f6)</sup>.
```
psm_effect_inds
├── figures
└── src
```

[sider](sider): Using the download files from Side Effect Resource (SIDERv4.1)<sup id="a3">[3](#f3)</sup><sup id="a4">[,4](#f4)</sup>, map drugs (excluding biologicals) to ChEMBL compound IDs and select the clinical-trial as opposed to post-marketing effects, save as pickle for further analysis.
Note: SIDER is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA 4.0), thus this same license is applied on the data and results in this directory.
```
sider
├── data
├── results
└── src
```

[faers_sider_comparison](faers_sider_comparison): Compare the significant AEs from FAERS (post-marketing) with those reported in SIDER (clinical).
```
faers_sider_comparison
├── data
├── figures
├── results
└── src
```

[bioactivities](bioactivities): Query experimental bioactivities from ChEMBL<sup id="a7">[7](#f7)</sup><sup id="a8">[,8](#f8)</sup> and use the tool PIDGIN<sup id="a9">[9](#f9)</sup><sup id="a10">[,10](#f10)</sup> to obtain ligand-based target predictions.
```
bioactivities
├── data
├── results
└── src
```

[plasma_concentrations](plasma_concentrations): Compile and process (standardise) drug plasma concentrations, both total and unbound, and integrate with fraction unbound/percentage plasma protein binding, from different sources.
```
plasma_concentrations
├── data
├── figures
├── results
│   └── interim
└── src
```

[integration_bioact_plasma_conc](integration_bioact_plasma_conc): Integrate the experimental and predicted bioactivities with the compiled plasma concentrations using the ratio of pXC50/Cmax to assign activity calls for drug-target pairs.
```integration_bioact_plasma_conc
├── figures
├── results
└── src
```

[prev_reported_safety_associations](prev_reported_safety_associations)
Extract protein targets and associated adverse effects from three reported safety target panels and map them to the MedDRA vocabulary.
```
prev_reported_safety_associations
├── data
└── src
```

[ae_target_links](ae_target_links)
Compute statistical associations between adverse events in FAERS or SIDER and protein targets.

```
ae_target_links
├── data
├── output
│   ├── 20200110_faers_cutoff6_pred_005_PRR2
│   │   ├── data
│   │   ├── logs
│   │   └── plots
│   ├── 20200110_faers_unbound_margin_pred_005_PRR2
│   │   ├── combinations
│   │   ├── data
│   │   ├── logs
│   │   └── plots
│   └── etc.
└── src
```

[analysis](analysis)
Various plots and high-level analyses (distributions, more association metrics) on the derived target-adverse event associations. Also analysis of combinations of targets.
```
analysis
├── data
├── results
│   ├── cutoff_pred_faers_vs_sider
│   ├── overarching
│   └── unbound_margin_pred_faers_vs_sider
└── src
```

Ines Smit  
[Contact](https://www.ch.cam.ac.uk/person/ias41 "Contact")

Released under the MIT license, except for the data in the [sider](sider) directory as specified in [sider/README.md](sider/README.md). 

##### Acknowledgements, references, and links

Supervised by Dr. Andreas Bender, University of Cambridge  
Funding by Lhasa Limited, additional supervision by Thierry Hanser from Lhasa Limited

<b id="f1">1</b> [Food and Drug Administration Adverse Event Reporting System download files](https://www.fda.gov/drugs/questions-and-answers-fdas-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-latest-quarterly-data-files) [↩](#a1)  

<b id="f2">2</b> Banda JM, Evans L, Vanguri RS, Tatonetti NP, Ryan PB, Shah NH. A curated and standardized adverse drug event resource to accelerate drug safety research. Sci Data. 2016;3:160026. Published 2016 May 10. doi: [10.1038/sdata.2016.26](https://doi.org/10.1038/sdata.2016.26) [↩](#a2)  

<b id="f3">3</b> [SIDER website](http://sideeffects.embl.de/) [↩](#a3)  

<b id="f4">4</b> Kuhn M, Letunic I, Jensen LJ, Bork P. The SIDER database of drugs and side effects. Nucleic Acids Res. 2016;44(D1):D1075‐D1079. doi: [10.1093/nar/gkv1075](https://doi.org/10.1093/nar/gkv1075) [↩](#a4)

<b id="f5">5</b> Tatonetti NP, Ye PP, Daneshjou R, Altman RB. Data-driven prediction of drug
effects and interactions. Sci Transl Med. 2012 Mar 14;4(125):125ra31. doi:
[10.1126/scitranslmed.3003377](https://doi.org/10.1126/scitranslmed.3003377) [↩](#a5)  

<b id="f6">6</b> [RxClass Overview](https://rxnav.nlm.nih.gov/RxClassIntro.html) [↩](#a6)  

<b id="f7">7</b> [ChEMBL website](https://www.ebi.ac.uk/chembl/) [↩](#a7)  

<b id="f8">8</b> Mendez D, Gaulton A, Bento AP, et al. ChEMBL: towards direct deposition of bioassay data. Nucleic Acids Res. 2019;47(D1):D930‐D940. doi: [10.1093/nar/gky1075](https://doi.org/10.1093/nar/gky1075) [↩](#a8)  

<b id="f9">9</b> [PIDGIN documentation](https://pidginv3.readthedocs.io/en/latest/) [↩](#a9)  

<b id="f10">10</b> Mervin LH, Bulusu KC, Kalash L, et al. Orthologue chemical space and its influence on target prediction. Bioinformatics. 2018;34(1):72‐79. doi: [10.1093/bioinformatics/btx525](https://doi.org/10.1093/bioinformatics/btx525) [↩](#a10)  
