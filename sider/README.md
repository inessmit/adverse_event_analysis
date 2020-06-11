Notebooks for extracting data from the downloaded files from SIDER v4.1<sup id="a1">[1](#f1)</sup><sup id="a2">[,2](#f2)</sup>. First maps the compounds (excluding biologicals) to ChEMBL ids using UniChem (matching on InChI keys)<sup id="a3">[3](#f3)</sup><sup id="a4">[,4](#f4)</sup>. Then selects the clinical trial (as opposed to post-marketing) adverse events -as far as currently possible- from SIDER and saves these only in a pickle. Thus a selection of the original data is made for further analysis.

SIDER is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA 4.0), thus this same license is applied on the data and result in this directory.

Also see [https://github.com/mkuhn/sider](https://github.com/mkuhn/sider).


<b id="f1">1</b> [SIDER website](http://sideeffects.embl.de/) [↩](#a1)  
<b id="f2">2</b> Kuhn M, Letunic I, Jensen LJ, Bork P. The SIDER database of drugs and side effects. Nucleic Acids Res. 2016;44(D1):D1075‐D1079. doi: [10.1093/nar/gkv1075](https://doi.org/10.1093/nar/gkv1075) [↩](#a2)  
<b id="f3">3</b> [UniChem](https://www.ebi.ac.uk/unichem/) [↩](#a3)  
<b id="f4">4</b> Chambers J, Davies M, Gaulton A, et al. UniChem: a unified chemical structure cross-referencing and identifier tracking system. J Cheminform. 2013;5(1):3. Published 2013 Jan 14. doi: [10.1186/1758-2946-5-3](https://doi.org/10.1186/1758-2946-5-3) [↩](#a4)  
