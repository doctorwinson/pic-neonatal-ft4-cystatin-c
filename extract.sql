-- PIC source queries. Executed by analysis.py in read-only transactions.

-- query: screening/connection
SELECT current_database(),current_setting('transaction_read_only') read_only;

-- query: screening/cohort
SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'));

-- query: screening/labs
WITH a AS (SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT l.*,d.label,d.label_cn FROM pic_icu.labevents l JOIN a USING(subject_id,hadm_id) JOIN pic_icu.d_labitems d USING(itemid)
WHERE l.itemid IN(6095,6193,6164,6165,6183,6180,6182,5032,5041,6954,5042,5075,6261,5024,6465,5626,5821,5227)
AND l.charttime BETWEEN a.admittime AND a.dischtime ORDER BY l.hadm_id,l.charttime,l.itemid;

-- query: screening/diagnoses
WITH a AS(SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT dx.*,d.icd10_code,d.title_cn,d.title FROM pic_icu.diagnoses_icd dx JOIN a USING(subject_id,hadm_id)
LEFT JOIN pic_icu.d_icd_diagnoses d ON d.icd10_code_cn=dx.icd10_code_cn;

-- query: screening/drugs
WITH a AS(SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT r.* FROM pic_icu.prescriptions r JOIN a USING(subject_id,hadm_id)
WHERE concat_ws(' ',r.drug_name,r.drug_name_cn,r.drug_name_generic) ~* '(thyro|甲状|优甲乐|左甲|dopamin|多巴胺|hydrocortis|氢化可的松|dexameth|地塞米松|methylpred|甲泼|prednis|泼尼|epinephrine|肾上腺素|dobutamin|多巴酚丁胺|milrinone|米力农|vancomy|万古霉素|gentamic|庆大霉素|amikacin|阿米卡星|furosemide|呋塞米|indometh|吲哚美辛|ibuprofen|布洛芬)';

-- query: screening/weights
WITH a AS(SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT c.* FROM pic_icu.chartevents c JOIN a USING(subject_id,hadm_id)
WHERE c.itemid='1014' AND c.charttime BETWEEN a.admittime AND a.admittime+interval '7 days';

-- query: screening/surgery
WITH a AS(SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT s.* FROM pic_icu.surgery_info s JOIN a USING(subject_id,hadm_id);

-- query: supplement/weights_all
WITH a AS (SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT c.* FROM pic_icu.chartevents c JOIN a USING(subject_id,hadm_id) WHERE c.itemid='1014' AND c.charttime BETWEEN a.admittime AND a.dischtime AND c.charttime<a.dob+interval '28 days';

-- query: supplement/rrt_dictionary
SELECT * FROM pic_icu.d_items WHERE concat_ws(' ',label,label_cn) ~* '(dialys|filtrat|透析|血滤|血液滤|肾脏替代|CRRT|ECMO)';

-- query: supplement/rrt_diagnoses
WITH a AS(SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT dx.*,d.title_cn,d.title,d.icd10_code FROM pic_icu.diagnoses_icd dx JOIN a USING(subject_id,hadm_id) LEFT JOIN pic_icu.d_icd_diagnoses d ON d.icd10_code_cn=dx.icd10_code_cn WHERE concat_ws(' ',d.title,d.title_cn) ~* '(dialys|透析|肾脏替代|ECMO)';

-- query: supplement/rrt_surgery
WITH a AS(SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT s.* FROM pic_icu.surgery_info s JOIN a USING(subject_id,hadm_id) WHERE s.surgery_name ~* '(dialys|hemofil|CRRT|ECMO|透析|血滤|血液滤|肾脏替代)';

-- query: metadata/connection
SELECT current_database(),current_setting('transaction_read_only') AS read_only;

-- query: metadata/columns
SELECT table_schema,table_name,column_name,data_type FROM information_schema.columns WHERE table_schema NOT IN ('information_schema','pg_catalog') ORDER BY table_schema,table_name,ordinal_position;

-- query: metadata/items
SELECT * FROM pic_icu.d_items ORDER BY itemid;

-- query: metadata/lab_items
SELECT * FROM pic_icu.d_labitems WHERE concat_ws(' ',label,label_cn) ~* '(lactate|lactic|乳酸|albumin|白蛋白|reactive|反应蛋白|血小板|platelet)';

-- query: metadata/neonatal_diagnosis_dictionary
SELECT * FROM pic_icu.d_icd_diagnoses WHERE icd10_code ~ '^P0[578]' OR concat_ws(' ',title,title_cn) ~* '(gestation|胎龄|孕周|出生体重|birth.?weight|small for gestational)';

-- query: audit/cohort
SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'));

-- query: audit/labs
WITH a AS (SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT l.*,d.label,d.label_cn FROM pic_icu.labevents l JOIN a USING(subject_id,hadm_id) JOIN pic_icu.d_labitems d USING(itemid)
WHERE l.itemid IN(6095,6193,6164,6165,6183,6180,6182,5032,5041,6954,5042,5075,6261,5024,6465,5626,5821,5227)
AND l.charttime BETWEEN a.admittime AND a.dischtime ORDER BY l.hadm_id,l.charttime,l.itemid;

-- query: audit/diagnoses
WITH a AS(SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT dx.*,d.icd10_code,d.title_cn,d.title FROM pic_icu.diagnoses_icd dx JOIN a USING(subject_id,hadm_id)
LEFT JOIN pic_icu.d_icd_diagnoses d ON d.icd10_code_cn=dx.icd10_code_cn;

-- query: audit/drugs
WITH a AS(SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT r.* FROM pic_icu.prescriptions r JOIN a USING(subject_id,hadm_id)
WHERE concat_ws(' ',r.drug_name,r.drug_name_cn,r.drug_name_generic) ~* '(thyro|甲状|优甲乐|左甲|dopamin|多巴胺|hydrocortis|氢化可的松|dexameth|地塞米松|methylpred|甲泼|prednis|泼尼|epinephrine|肾上腺素|dobutamin|多巴酚丁胺|milrinone|米力农|vancomy|万古霉素|gentamic|庆大霉素|amikacin|阿米卡星|furosemide|呋塞米|indometh|吲哚美辛|ibuprofen|布洛芬)';

-- query: audit/surgery
WITH a AS(SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT s.* FROM pic_icu.surgery_info s JOIN a USING(subject_id,hadm_id);

-- query: audit/connection
SELECT current_database(),current_setting('transaction_read_only') AS read_only;

-- query: audit/icustays
WITH a AS(SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT i.* FROM pic_icu.icustays i JOIN a USING(subject_id,hadm_id);

-- query: audit/all_admissions_of_subjects
WITH a AS(SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT ad.subject_id,ad.hadm_id,ad.admittime,ad.dischtime FROM pic_icu.admissions ad WHERE ad.subject_id IN (SELECT subject_id FROM a);

-- query: audit/thyroid_drugs_expanded
WITH a AS(SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT r.* FROM pic_icu.prescriptions r JOIN a USING(subject_id,hadm_id) WHERE concat_ws(' ',drug_name,drug_name_cn,drug_name_generic) ~* '(thyro|甲状|优甲乐|左甲|methimazole|thiamazole|carbimazole|propylthiouracil|甲巯咪唑|丙硫氧嘧啶|卡比马唑)';

-- query: audit/thyroid_renal_dictionary
SELECT * FROM pic_icu.d_icd_diagnoses WHERE trim(icd10_code) ~ '^(E0[0-7]|Q6[0-4])' OR concat_ws(' ',title,title_cn) ~* '(thyroid|甲状腺|肾发育|肾缺如|多囊肾|先天.*肾|肾.*先天)';

-- query: audit/rrt_ecmo_dictionary
SELECT * FROM pic_icu.d_items WHERE concat_ws(' ',label,label_cn) ~* '(dialys|hemofil|ECMO|CRRT|透析|血滤|血液滤|肾脏替代)';

-- query: audit/rrt_ecmo_surgery
WITH a AS(SELECT a.*,p.dob,p.gender FROM pic_icu.admissions a JOIN pic_icu.patients p USING(subject_id)
WHERE a.admittime>=p.dob AND a.admittime<p.dob+interval '28 days'
AND EXISTS(SELECT 1 FROM pic_icu.icustays i WHERE i.hadm_id=a.hadm_id AND (i.first_careunit='NICU' OR i.last_careunit='NICU'))) SELECT s.* FROM pic_icu.surgery_info s JOIN a USING(subject_id,hadm_id) WHERE s.surgery_name ~* '(dialys|hemofil|CRRT|ECMO|透析|血滤|血液滤|肾脏替代)';

-- query: covariates/connection
SELECT current_database(),current_setting('transaction_read_only') AS read_only;

-- query: covariates/symptoms
WITH c AS (SELECT * FROM jsonb_to_recordset(%(paired_tests)s::jsonb) AS pairs(subject_id integer,hadm_id integer,thyroid_time timestamp,renal_time timestamp)) SELECT e.* FROM pic_icu.emr_symptoms e JOIN c USING(subject_id,hadm_id);

-- query: covariates/diagnoses
WITH c AS (SELECT * FROM jsonb_to_recordset(%(paired_tests)s::jsonb) AS pairs(subject_id integer,hadm_id integer,thyroid_time timestamp,renal_time timestamp)) SELECT e.*,d.icd10_code,d.title,d.title_cn FROM pic_icu.diagnoses_icd e JOIN c USING(subject_id,hadm_id) LEFT JOIN pic_icu.d_icd_diagnoses d USING(icd10_code_cn);

-- query: covariates/vitals
WITH c AS (SELECT * FROM jsonb_to_recordset(%(paired_tests)s::jsonb) AS pairs(subject_id integer,hadm_id integer,thyroid_time timestamp,renal_time timestamp)) SELECT e.* FROM pic_icu.chartevents e JOIN c USING(subject_id,hadm_id) JOIN pic_icu.admissions a USING(subject_id,hadm_id) WHERE e.itemid IN ('1001','1002','1003','1004','1006','1015','1016') AND e.charttime BETWEEN greatest(a.admittime,least(c.thyroid_time,c.renal_time)-interval '24 hours') AND least(c.thyroid_time,c.renal_time);

-- query: covariates/weight
WITH c AS (SELECT * FROM jsonb_to_recordset(%(paired_tests)s::jsonb) AS pairs(subject_id integer,hadm_id integer,thyroid_time timestamp,renal_time timestamp)) SELECT e.*,p.dob FROM pic_icu.chartevents e JOIN c USING(subject_id,hadm_id) JOIN pic_icu.patients p USING(subject_id) JOIN pic_icu.admissions a USING(subject_id,hadm_id) WHERE e.itemid='1014' AND e.charttime BETWEEN a.admittime AND least(c.thyroid_time,c.renal_time);

-- query: covariates/labs
WITH c AS (SELECT * FROM jsonb_to_recordset(%(paired_tests)s::jsonb) AS pairs(subject_id integer,hadm_id integer,thyroid_time timestamp,renal_time timestamp)) SELECT e.* FROM pic_icu.labevents e JOIN c USING(subject_id,hadm_id) JOIN pic_icu.admissions a USING(subject_id,hadm_id) WHERE e.itemid IN (5227,5024,6465,5626,5821) AND e.charttime BETWEEN greatest(a.admittime,least(c.thyroid_time,c.renal_time)-interval '24 hours') AND least(c.thyroid_time,c.renal_time);

-- query: covariates/drugs
WITH c AS (SELECT * FROM jsonb_to_recordset(%(paired_tests)s::jsonb) AS pairs(subject_id integer,hadm_id integer,thyroid_time timestamp,renal_time timestamp)) SELECT e.* FROM pic_icu.prescriptions e JOIN c USING(subject_id,hadm_id) WHERE e.startdate<=least(c.thyroid_time,c.renal_time) AND (e.enddate>=least(c.thyroid_time,c.renal_time)-interval '48 hours' OR e.enddate IS NULL) AND concat_ws(' ',e.drug_name,e.drug_name_cn,e.drug_name_generic) ~* '(dopamin|多巴胺|hydrocortis|氢化可的松|dexameth|地塞米松|methylpred|甲泼|prednis|泼尼|epinephrine|肾上腺素|dobutamin|多巴酚丁胺|milrinone|米力农)';
