use chembl_25;
select distinct
molregno
, pref_name
, withdrawn_year
, withdrawn_reason
, withdrawn_country
from molecule_dictionary
where withdrawn_flag = 1;