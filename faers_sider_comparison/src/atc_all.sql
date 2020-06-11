use chembl_25;;
select distinct
  atc.molregno
, atc_desc.level5
, level1_description
, level2_description
, level3_description
, level4_description
from atc_classification atc_desc
join molecule_atc_classification atc on atc.level5 = atc_desc.level5;
