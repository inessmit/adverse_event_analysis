select distinct con.concept_id, con.concept_name, con.concept_code 
from FAERS_AEOLUS.concept con
join FAERS_AEOLUS.standard_case_drug drug  on drug.standard_concept_id = con.concept_id
where con.vocabulary_id = 'RxNorm'
and con.concept_class_id = 'Ingredient'
and con.standard_concept = 'S';
