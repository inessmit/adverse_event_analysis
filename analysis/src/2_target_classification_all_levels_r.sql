use chembl_24_1;
select
comp_seqs.accession
, case
when prot_class.class_level = 1 then prot_class.pref_name
when prot_class.class_level = 2 then prot_class2.pref_name
when prot_class.class_level = 3 then prot_class3.pref_name
when prot_class.class_level = 4 then prot_class4.pref_name
when prot_class.class_level = 5 then prot_class5.pref_name
when prot_class.class_level = 6 then prot_class6.pref_name
end as level_1
, case
when prot_class.class_level = 1 then 'Not available'
when prot_class.class_level = 2 then prot_class.pref_name
when prot_class.class_level = 3 then prot_class2.pref_name
when prot_class.class_level = 4 then prot_class3.pref_name
when prot_class.class_level = 5 then prot_class4.pref_name
when prot_class.class_level = 6 then prot_class5.pref_name
end as level_2
, case
when prot_class.class_level = 3 then prot_class.pref_name
when prot_class.class_level = 4 then prot_class2.pref_name
when prot_class.class_level = 5 then prot_class3.pref_name
when prot_class.class_level = 6 then prot_class4.pref_name
else NULL 
end as level_3
, case
when prot_class.class_level = 4 then prot_class.pref_name
when prot_class.class_level = 5 then prot_class2.pref_name
when prot_class.class_level = 6 then prot_class3.pref_name
else NULL
end as level_4
, case
when prot_class.class_level = 5 then prot_class.pref_name
when prot_class.class_level = 6 then prot_class2.pref_name
else NULL
end as level_5
, case
when prot_class.class_level = 6 then prot_class.pref_name
else NULL
end as level_6
from 
protein_classification prot_class
join component_class comp_class on comp_class.protein_class_id = prot_class.protein_class_id
join component_sequences comp_seqs on comp_seqs.component_id = comp_class.component_id
left join protein_classification prot_class2 on prot_class2.protein_class_id = prot_class.parent_id
left join protein_classification prot_class3 on prot_class3.protein_class_id = prot_class2.parent_id
left join protein_classification prot_class4 on prot_class4.protein_class_id = prot_class3.parent_id
left join protein_classification prot_class5 on prot_class5.protein_class_id = prot_class4.parent_id
left join protein_classification prot_class6 on prot_class6.protein_class_id = prot_class5.parent_id;


