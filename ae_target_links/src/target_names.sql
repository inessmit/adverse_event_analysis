use chembl_24_1;
select distinct
td.tid
, td.pref_name
, td.target_type
, cs.accession
, cs.organism as accession_organism
, td.organism as target_organism
from target_dictionary td
join target_components tc on td.tid = tc.tid
join component_sequences cs on cs.component_id = tc.component_id
where td.target_type in ('SINGLE PROTEIN')
and cs.organism = 'Homo sapiens';
