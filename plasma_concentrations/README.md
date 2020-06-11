Set of notebooks for processing drug plasma concentrations, both total and unbound, compiled from different sources. Combines data from publications and from ChEMBL<sup id="a1">[1](#f1)</sup><sup id="a2">[,2](#f2)</sup>.

One major source of total plasma concentrations is the publication by Schulz et al.<sup id="a3">[3](#f3)</sup>.

The upper values from the above publication were extracted from its supplementary data (data/Schulz_ea_upper_values_formatted.txt)

The set of notebooks in the 'src' directory combine this with additional data from ChEMBl, and then take the median of various measurements available per drug.

Then unbound concentrations are calculated using Fraction unbound (Fu) and % Plasma Protein Binding (PPB) data from range of sources, and again the median is taken per compound.

Main results in the 'results' directory.

More details in the respective notebooks.

To recreate: 

The ChEMBL queries (using ChEMBL version 24.1) were executed on a local MySQL installation and results saved in 'data' directory.

Execute notebooks in 'src' in order. 

Use environment as in top-level requirements.txt, except for notebook 3 'format_extracted_unbound_datasets' which uses RDKit environment.
```
conda create -c rdkit -n my-rdkit-env rdkit python=3.6 xlrd nb_conda_kernels pymysql
```

<b id="f1">1</b> [ChEMBL website](https://www.ebi.ac.uk/chembl/) [↩](#a1)  

<b id="f2">2</b> Mendez D, Gaulton A, Bento AP, et al. ChEMBL: towards direct deposition of bioassay data. Nucleic Acids Res. 2019;47(D1):D930‐D940. doi: [10.1093/nar/gky1075](https://doi.org/10.1093/nar/gky1075) [↩](#a2)  

<b id="f3">3</b> Schulz M, Iwersen-Bergmann S, Andresen H, Schmoldt A. Therapeutic and toxic
blood concentrations of nearly 1,000 drugs and other xenobiotics. Crit Care. 2012
Jul 26;16(4):R136. doi: [10.1186/cc11441](https://doi.org/10.1186/cc11441) [↩](#a3)  