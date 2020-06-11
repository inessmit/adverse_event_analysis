Scripts to retrieve drug indications from RxClass API<sup id="a1">[1](#f1)</sup>. These will be used in subsequent part as independent dataset of drug indications to check the effect of Propsensity Score Matching technique on indication bias.

Some manual steps required so divided in three parts.

Needs the mapping sqlite db created in other [compound_mapping](../compound_mapping) dir.
Also needs db defaults file (explained in README of [psm_aeolus](../psm_aeolus))

Part 1:

"""Get all MedDRA terms occurring in FAERS AEOLUS and save to excel file for hierarchy analysis in MedDRA web-based browser<sup id="a2">[2](#f2)</sup>."""

Run the run_part1.sh script. Then do Hierarchy Analysis on all adverse event terms (data/all_aeolus_aes_hierarchy_input.xlsx file) in MedDRA web-based browser<sup id="a2">[2](#f2)</sup>. Define the output file in the run_part3.sh script.

Part 2 :

"""For each RxNorm id, get drug indications from RxClass API. Use MRCONSO file downloaded from UMLS<sup id="a3">[3](#f3)</sup> to map to MedDRA terms. Save file with MedDRA terms for hierarchy analysis."""

Run the run_part2.sh script. Then do Hierarchy Analysis on the MedDRA terms for indications (data/meddra_terms_indications_hierarchy_input.xlsx)in MedDRA web-based browser. Define the outputfile in the run_part3.sh script.

Part 3:

"""After MedDRA hierarchy analysis on indication terms and ae terms, extract HLT for each indication and retrieve other PTs in that HLT (from occurring AEs in AEOLUS). Save dictionary with molregno to MedDRA PT/LLT indications, expanded via HLT. """

Run the run_part3.sh script.

<b id="f1">1</b> [RxClass Overview](https://rxnav.nlm.nih.gov/RxClassIntro.html) [↩](#a1)  

<b id="f2">2</b> [https://www.meddra.org/](https://www.meddra.org/) [↩](#a2)  

<b id="f3">3</b> [https://www.nlm.nih.gov/research/umls/index.html](https://www.nlm.nih.gov/research/umls/index.html) [↩](#a3)  