use chembl_24_1;
select count(distinct parent_molregno)
from molecule_hierarchy hier
join molecule_dictionary md on md.molregno = hier.molregno
join activities act on act.molregno = md.molregno
join assays assays on assays.assay_id = act.assay_id
where act.standard_type = 'Cmax'
and assays.assay_organism = 'Homo sapiens'
and (upper(assays.assay_tissue) in ('BLOOD', 'SERUM', 'PLASMA') OR assays.assay_tissue is NULL)
and act.standard_flag = 1
group by max_phase;
