use chembl_24_1;
select count(distinct parent_molregno)
from molecule_hierarchy hier
join molecule_dictionary md on md.molregno = hier.molregno
join activities act on act.molregno = md.molregno
join assays on assays.assay_id = act.assay_id
join target_dictionary td on assays.tid = td.tid
join target_components tc on td.tid = tc.tid
join component_sequences cs on cs.component_id = tc.component_id
where act.standard_type in ('IC50','EC50','XC50','AC50','Ki','Kd','Potency')
and assays.assay_type in ('B', 'F')
and assays.confidence_score in (7,9)
and td.target_type in ('SINGLE PROTEIN', 'PROTEIN COMPLEX')
and cs.tax_id = 9606
group by max_phase;
