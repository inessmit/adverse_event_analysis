select concept_id, concept_name
from concept
where concept_id in (select distinct indication_concept_id from standard_case_indication);
