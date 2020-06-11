Notebooks for mapping drug concepts in the FAERS AEOLUS database to ChEMBL identifiers/molecular structures via RxNorm, UniChem, and synonyms. Some of the notebooks required manual inspection/selection of the mapping.

Description of notebooks and what is done in each:

**1. faers_drug_concepts.sql and make_mapping_db.ipynb**
- RxNorm identifiers from FAERS
- prepare mapping DB 'drug_mapping_process.db'

**2. Download RxNorm flat files into data directory**
This is the RxNorm hierarchy. Store in data/raw directory.

**3. remapped_retired**
Check version of RxNorm API and get remapped and retired concepts
- Recent version of RxNorm API record version for changes compared to current RxNorm CONSO file
- Use RELA file to get parents of salt forms (for all rxnorm concepts)
- Use CONSO file to get DrugBank IDS for parents if available, otherwise the rxnorm concept

**4. drugbank_inchi_unichem.ipynb**
- record Unichem, versions of DrugBank and ChEMBL and date
- For the concepts with a Drugbank ID, retrieve inchi and inchikey from UniChem webservices and insert into sqlite mapping db.
- For the 'errors' in retrieving inchi, check if they are biological and flag them as such in the mapping db
- also adds 'drugbank_type' for the rows with drugbank_id

**5. FAERS_compound_mapping.ipynb**
Use ChEMBL stored database get parent ID based on InchI for each of the compounds
- Inspect the retrieved mappings (all done via Drugbank mappings)
if they are not max_phase=4, try to see if mapping of rxnorm name on chembl pref_name, compound_name, or synonyms matches another compound with max_phase 4 
if that is the case, prefer that compound. 
- First without structures, then retrieve structures (not inner join)
- makes new database for mapped compounds: mapped_compounds.db

**6. FAERS_synonym_mapping.ipynb**
- For compounds left over in drug_mapping_process.db (unmapped) (including those with DrugBank ID where inchi was not in ChEMBL, and compounds without Drugbank ID, and retired RxNorm concepts), try to match on pref names, compound_names and synonyms to ChEMBL
- first without structures, then retrieve structures (not inner join)
- Following hierarchy of synonym mappings is done: 
    1. RxNorm name against ChEMBL pref_name where max_phase = 4 (only approved compounds)
    2. RxNorm name against compound_name,  max_phase = 4
    3. RxNorm name against ChEMBL synonyms (excluding TRADE_NAME), max_phase = 4  
    Note: for 1, 2 and 3: if multiple compounds with that name, manually inspect if there was more than one match
    4. RxNorm name against pref_name, not max_phase = 4
    5. RxNorm MTHSPL (structured product label approved by FDA) synonym against ChEMBL pref_name, not max_phase = 4
    6. RxNorm synonyms (excluding NDFRT, SNOMED_US, CVX) against ChEMBL synonyms, max_phase = 4
      - Uses chembl database to get ChEMBL parent ID and InchI for these rest of the compounds
      - Updates apped_compounds.db with additional mappings

Notes:
- 201903_drug_mapping_process is used in the PROCESS only - not final results
- 201903_mapped_compounds has the final mappings and structures
- FAERS_compound mapping was done on local computer, using sql queries and scp
- then FAERS synonym mapping on calculon (server), using pymysql to connect to calculon directly
- ideally should change FAERS_compound_mapping script to use same technique as FAERS synonym mapping.
