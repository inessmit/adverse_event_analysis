use chembl_24_1;
select distinct mh.parent_molregno
, md.chembl_id as parent_chembl_id
, md.pref_name
, md.max_phase
, act.standard_type
, act.standard_relation
, act.standard_value
, act.upper_value
, act.standard_flag
, assays.assay_test_type
, assays.description
, assays.assay_organism
, assays.assay_tissue
, act.activity_comment
, act.data_validity_comment
, assays.chembl_id as assay_id
from activities act
join molecule_hierarchy mh on mh.molregno = act.molregno
join assays on assays.assay_id = act.assay_id
join molecule_dictionary md on md.molregno = mh.parent_molregno
where act.standard_type in ('Fu', 'PPB') 
and act.potential_duplicate=0;
