select distinct isr, indication_concept_id 
from standard_case_indication 
where isr != '' 
and isr in (select distinct isr from standard_case_outcome where isr is not null);
