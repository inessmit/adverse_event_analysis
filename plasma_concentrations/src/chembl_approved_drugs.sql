use chembl_24_1;
select distinct
md.chembl_id
, md.molregno
, md.pref_name
, cr.compound_name
, cr.compound_key
, syn.synonyms
, syn.syn_type
, str.standard_inchi_key
, str.canonical_smiles
, prop.mw_freebase
from molecule_dictionary md
left join compound_structures str on str.molregno = md.molregno
left join compound_records cr on cr.molregno = md.molregno
left join compound_properties prop on prop.molregno = md.molregno
left join molecule_synonyms syn on syn.molregno = md.molregno
where md.max_phase = 4;
