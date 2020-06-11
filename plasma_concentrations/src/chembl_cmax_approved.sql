use chembl_24_1;
select distinct
act.activity_id
, mh.parent_molregno
, md2.chembl_id as parent_chembl_id
, props.mw_freebase
, md2.pref_name as parent_pref_name
, md.molregno as version_molregno
, md.chembl_id as version_chembl_id
, md.pref_name as version_pref_name
, act.published_value
, act.published_units
, act.standard_type
, act.standard_value
, act.standard_upper_value
, act.standard_units
, act.standard_text_value
, act.data_validity_comment
, act.activity_comment
, assays.description
, assays.assay_organism 
, assays.assay_tissue
, assays.assay_cell_type
, assays.chembl_id as assay_chembl_id
, src.src_description
, docs.pubmed_id
, docs.title
from molecule_hierarchy mh
join molecule_dictionary md on md.molregno = mh.molregno
join activities act on act.molregno = md.molregno
join molecule_dictionary md2 on mh.parent_molregno = md2.molregno
join assays assays on assays.assay_id = act.assay_id
join source src on src.src_id = assays.src_id
join docs docs on assays.doc_id = docs.doc_id
left join compound_properties props on props.molregno = md2.molregno
where md2.max_phase = 4
and act.standard_flag=1
and act.standard_type = 'Cmax';

