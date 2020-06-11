These scripts are an implementation of the propensity score matching (PSM) technique described by Tatonetti et al.<sup id="a1">[1](#f1)</sup> on the curated AEOLUS version of FAERS<sup id="a2">[2](#f2)</sup>.
The scripts calculate statistical associations between AEOLUS/RxNorm drugs and adverse events using the Proportional Reporting Ratio (PRR) after using PSM and also without PSM as comparison.

The work was inspired by and builds on the work of Tatonetti et al.<sup id="a1">[1](#f1)</sup><sup id="a2">[3](#f2)</sup> and the PhD thesis by Tatonetti<sup id="a1">[4](#f4)</sup>.
The only intended differences with the original description of PSM by Tatonetti et al.<sup id="a1">[1](#f1)</sup> are that Chi-squared test instead of the Spearman correlation coefficient (rho) is used for identifying features correlated with drug prescription that are used in the PSM. Furthermore, the code is Python only and specific for the AEOLUS database.

Everything was ran using the 'release' environment specified in the README and requirements.txt file in the root directory.

Everything was ran on a server running Ubuntu 18.04 with access to multiple CPU cores and 256 GB of RAM. Running the runall.sh script took over 24h with this setup but also depends on the number of drugs.

The scripts need access to a MySQL installation of the FAERS AEOLUS database, installed according to the scripts by the authors<sup id="a2">[,2](#f2)</sup>.
We noticed the database lacked a few indeces that are needed to execute the code in this repository. See src/aeolus_create_additional_indexes.sql file.

The MySQL login details and parameters were specified in a file with the contents:
```
[client]
user="""
password=""
database='FAERS_AEOLUS'
port=3306
host='localhost'
```

Specify the location of this file in the runall.sh script, and specify the basedir. Then running the runall.sh script does everything.
The drugs that are analysed can be restricted to a subset by specifying a txt file with AEOLUS drug concept ids in the runall.sh script. The list of ids can be prepared depending on overlap with ChEMBL bioactivity data in the aeolus_ids_for_ae_detection.ipynb notebook.


The final results will be placed in 'results' directory. The main output files for the drug-adverse event associatons are the files in results/data beginning with an AEOLUS drug id and ending with ae_associations.txt.
The results/figures contains more intermediate results, i.e. figures about the propensity score matching.


<b id="f1">1</b> Tatonetti NP, Ye PP, Daneshjou R, Altman RB. Data-driven prediction of drug
effects and interactions. Sci Transl Med. 2012 Mar 14;4(125):125ra31. doi:
[10.1126/scitranslmed.3003377](https://doi.org/10.1126/scitranslmed.3003377) [↩](#a1)  
<b id="f2">2</b> Banda JM, Evans L, Vanguri RS, Tatonetti NP, Ryan PB, Shah NH. A curated and
standardized adverse drug event resource to accelerate drug safety research. Sci 
Data. 2016 May 10;3:160026. doi: [10.1038/sdata.2016.26](https://doi.org/10.1038/sdata.2016.26) [↩](#a2)  
<b id="f3">3</b> [http://tatonettilab.org/offsides/](http://tatonettilab.org/offsides/) [↩](#a3)  
<b id="f4">4</b> [https://searchworks.stanford.edu/view/9625324](https://searchworks.stanford.edu/view/9625324) [↩](#a4)  
