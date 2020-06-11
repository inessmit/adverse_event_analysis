use chembl_25;
select distinct
mh.parent_molregno
, md.chembl_id as parent_chembl_id
, cs.accession
, td.tid
, td.organism as target_organism
, td.pref_name as target_name
, td.target_type
, act.standard_type
, act.standard_relation
, act.pchembl_value
, act.standard_flag
, act.activity_comment
, act.data_validity_comment
, assays.description
, assays.chembl_id as assay_chembl_id
, src.src_id
, src.src_description
from activities act
join molecule_hierarchy mh on mh.molregno = act.molregno
join assays on assays.assay_id = act.assay_id
join target_dictionary td on assays.tid = td.tid
join target_components tc on td.tid = tc.tid
join component_sequences cs on cs.component_id = tc.component_id
join molecule_dictionary md on md.molregno = mh.parent_molregno
join source src on src.src_id = assays.src_id
and act.standard_type in ('IC50','EC50','XC50','AC50','Ki','Kd','Potency')
and act.potential_duplicate = 0
and assays.confidence_score in (7,9)
and assays.assay_type in ('B', 'F')
and td.target_type in ('SINGLE PROTEIN', 'PROTEIN COMPLEX')
and td.tax_id = 9606 -- Human;