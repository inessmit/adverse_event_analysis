select distinct primaryid, standard_concept_id 
from standard_case_drug 
where primaryid != '' 
and primaryid in (select distinct primaryid from standard_case_outcome where primaryid is not null);
