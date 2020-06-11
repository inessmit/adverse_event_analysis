select distinct isr, standard_concept_id 
from standard_case_drug 
where isr != '' 
and isr in (select distinct isr from standard_case_outcome where isr is not null);
