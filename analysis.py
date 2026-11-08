"""PIC neonatal FT4-cystatin C analyses.

Run with --extract to read PostgreSQL; otherwise use the existing local snapshots.
Primary definitions were frozen after initial results on 15 October 2026.
Supplementary spline, optimizer and selection checks were added on 16 October 2026.
"""
from pathlib import Path
from datetime import datetime
import argparse, hashlib, json, os, re, subprocess, sys

HERE = Path(__file__).resolve().parent
WORK = HERE / 'results'
AGE_KNOTS = [5.13230671301, 10.8030960648, 23.18747222222]
SEED = 20261015

def initialize():
    for name in ['screening', 'initial_analysis', 'covariate_analysis', 'primary_analysis',
                 'supplementary_analysis/evidence', 'selection_analysis/evidence', 'verification/evidence']:
        (WORK / name).mkdir(parents=True, exist_ok=True)

def read_queries():
    text = (HERE / 'extract.sql').read_text(encoding='utf-8')
    chunks = re.split(r'^-- query: ([a-z_]+/[a-z_]+)\s*$', text, flags=re.M)
    return dict(zip(chunks[1::2], [x.strip().rstrip(';') for x in chunks[2::2]]))

def extract_group(group, destination, parameters=None, nested=True):
    import psycopg2
    from psycopg2.extras import RealDictCursor
    connection = psycopg2.connect(
        host=os.getenv('PGHOST', 'localhost'), port=os.getenv('PGPORT', '5432'),
        dbname=os.getenv('PGDATABASE', 'pic'), user=os.getenv('PGUSER', 'postgres'),
        password=os.environ['PGPASSWORD'],
        options='-c default_transaction_read_only=on -c statement_timeout=120000')
    connection.set_session(readonly=True, isolation_level='REPEATABLE READ')
    queries = {}
    try:
        with connection.cursor(cursor_factory=RealDictCursor) as cursor:
            for key, statement in read_queries().items():
                prefix, name = key.split('/')
                if prefix != group:
                    continue
                cursor.execute(statement, parameters) if parameters else cursor.execute(statement)
                queries[name] = {'sql': statement, 'rows': [dict(row) for row in cursor.fetchall()]}
        connection.rollback()
    finally:
        connection.close()
    now = datetime.now().astimezone().isoformat()
    output = {'executed_at': now, 'extracted_at': now, 'queries': queries} if nested else queries
    destination.write_text(json.dumps(output, ensure_ascii=False, default=str, indent=2), encoding='utf-8')

def extract_sources():
    extract_group('screening', WORK/'screening/feasibility_raw.json')
    extract_group('supplement', WORK/'screening/supplement_raw.json', nested=False)
    extract_group('metadata', WORK/'covariate_analysis/metadata_audit.json')
    extract_group('audit', WORK/'primary_analysis/reaudit_source_local.json')

def extract_covariates():
    cohort = json.loads((WORK/'screening/cystatin_primary_6h_patient_audit.json').read_text(encoding='utf-8'))
    pairs = [dict(subject_id=r['subject_id'], hadm_id=r['hadm_id'],
                  thyroid_time=r['charttime'], renal_time=r['renal_time']) for r in cohort]
    extract_group('covariates', WORK/'covariate_analysis/extended_source_snapshot.json',
                  {'paired_tests': json.dumps(pairs)})


def build_pairs(window=6, write_outputs=True):
    import json,re,sys
    from pathlib import Path
    import pandas as pd
    ROOT=WORK/'screening'
    WINDOW=window
    E=json.loads((ROOT/'feasibility_raw.json').read_text(encoding='utf-8'))
    S=json.loads((ROOT/'supplement_raw.json').read_text(encoding='utf-8'))
    A=pd.DataFrame(E['queries']['cohort']['rows'])
    for c in ['admittime','dischtime','dob','deathtime']:A[c]=pd.to_datetime(A[c])
    A=A.sort_values(['admittime','hadm_id']).drop_duplicates('subject_id')
    L=pd.DataFrame(E['queries']['labs']['rows']);L['charttime']=pd.to_datetime(L.charttime)
    L['numeric']=L.valuenum.notna()&L.valuenum.gt(0)&~L.value.fillna('').str.contains(r'[<>]|LESS|GREATER',regex=True,case=False)
    mapping={6095:'TSH',6193:'TSH',6164:'FT3',6165:'FT4',6183:'FT4',5032:'Cr',5041:'Cr',6954:'Cr',6261:'CysC',5024:'albumin',6465:'albumin',5626:'CRP',5821:'CRP'}
    V=L[L.itemid.isin(mapping)&L.numeric].copy();V['test']=V.itemid.map(mapping)
    V=V.groupby(['subject_id','hadm_id','charttime','test']).valuenum.median().unstack('test').reset_index()
    V=V.merge(A[['subject_id','hadm_id','admittime','dischtime','dob','gender']],on=['subject_id','hadm_id'],validate='many_to_one')
    V['age_days']=(V.charttime-V.dob).dt.total_seconds()/86400
    D=pd.DataFrame(E['queries']['diagnoses']['rows']);flags={}
    for h,g in D.groupby('hadm_id'):
     text=' '.join(g.title_cn.fillna(''));codes=' '.join(g.icd10_code.fillna(''))
     flags[h]={'thyroid_dx':bool(re.search(r'甲状腺|甲状腺素',text)),'renal_malformation':bool(re.search(r'Q6[0-4]',codes)),
     'preterm_recorded':bool(re.search(r'P07\.[23]',codes)),'lbw_recorded':bool(re.search(r'P07\.[01]',codes)),'sepsis_recorded':bool(re.search(r'P36|A4[01]',codes))}
    R=pd.DataFrame(E['queries']['drugs']['rows']);R['startdate']=pd.to_datetime(R.startdate);R['enddate']=pd.to_datetime(R.enddate)
    R['text']=R[['drug_name','drug_name_cn','drug_name_generic']].fillna('').agg(' '.join,axis=1)
    RG={h:g for h,g in R.groupby('hadm_id')}
    W=pd.DataFrame(S['weights_all']['rows']);W['charttime']=pd.to_datetime(W.charttime);W=W[W.valuenum.gt(0)]
    WG={h:g for h,g in W.groupby('hadm_id')};VG={h:g for h,g in V.groupby('hadm_id')}
    thy=V[(V.age_days>=3)&(V.age_days<28)].dropna(subset=['TSH','FT4']).copy()
    renal=V[(V.age_days>=3)&(V.age_days<28)].dropna(subset=['CysC','Cr']).copy()
    renal_groups={h:g for h,g in renal.groupby('hadm_id')}
    matched=[]
    for ix,r in thy.iterrows():
     g=renal_groups.get(r.hadm_id)
     if g is None:continue
     g=g.copy();g['lag']=(g.charttime-r.charttime).dt.total_seconds()/3600
     g=g[g.lag.abs()<=WINDOW].sort_values(['charttime']);g=g.iloc[g.lag.abs().argsort(kind='stable')]
     if len(g):
      hit=g.iloc[0];row=r.to_dict();row.update({'CysC':hit.CysC,'Cr':hit.Cr,'renal_time':hit.charttime,'lag_hours':hit.lag});matched.append(row)
    F=pd.DataFrame(matched)
    rows=[]
    for r in F.itertuples():
     flags_h=flags.get(r.hadm_id,{})
     drugs=RG.get(r.hadm_id,R.iloc[:0]);prior=drugs[drugs.startdate<=r.charttime]
     recent=prior[(prior.enddate>=r.charttime-pd.Timedelta(hours=48))|prior.enddate.isna()]
     wg=WG.get(r.hadm_id,W.iloc[:0]);weight=wg[(wg.charttime<=r.charttime)&(wg.charttime>=r.charttime-pd.Timedelta(days=7))]
     vg=VG[r.hadm_id];cov=vg[(vg.charttime<=r.charttime)&(vg.charttime>=r.charttime-pd.Timedelta(hours=24))]
     row={c:getattr(r,c) for c in ['subject_id','hadm_id','charttime','renal_time','lag_hours','age_days','gender','TSH','FT4','FT3','CysC','Cr']}
     row.update(flags_h)
     extracorp=any(s['hadm_id']==r.hadm_id and pd.to_datetime(s['surgery_begin_time'])<=max(r.charttime,r.renal_time) for s in S['rrt_surgery']['rows'])
     row.update({'extracorporeal_before':extracorp,'thyroxine_before':bool(prior.text.str.contains(r'thyro|甲状|优甲乐|左甲',case=False,regex=True).any()),
     'vaso_recent':bool(recent.text.str.contains(r'dopamin|多巴胺|epinephrine|肾上腺素|dobutamin|多巴酚丁胺|milrinone|米力农',case=False,regex=True).any()),
     'steroid_recent':bool(recent.text.str.contains(r'hydrocortis|氢化可的松|dexameth|地塞米松|methylpred|甲泼|prednis|泼尼',case=False,regex=True).any()),
     'weight_record':len(weight)>0,'albumin_24h':cov.albumin.notna().any(),'CRP_24h':cov.CRP.notna().any()})
     rows.append(row)
    F=pd.DataFrame(rows).sort_values(['subject_id','charttime'])
    flows={'neonatal_NICU_admissions':len(E['queries']['cohort']['rows']),'first_admission_subjects':len(A),
     'thyroid_same_timestamp_age3to28_subjects':int(thy.subject_id.nunique()),'paired_within_window_subjects':int(F.subject_id.nunique()),'paired_within_window_records':len(F)}
    for label,mask in [('no_thyroid_dx',~F.thyroid_dx),('no_renal_malformation',~F.renal_malformation),('no_prior_thyroxine',~F.thyroxine_before),('no_prior_extracorporeal',~F.extracorporeal_before)]:
     F=F[mask.reindex(F.index)];flows[label+'_subjects']=int(F.subject_id.nunique())
    if not write_outputs:
     return {'F':F, 'E':E, 'S':S, 'L':L}
    I=F.drop_duplicates('subject_id')
    repeat=[]
    for sub,g in F.groupby('subject_id'):
     last=None;selected=[]
     for ix,r in g.iterrows():
      if last is None or r.charttime>=last+pd.Timedelta(hours=48):selected.append(ix);last=r.charttime
     repeat.extend(selected)
    J=F.loc[repeat];nrep=J.groupby('subject_id').size()
    summary={'pairing_window_hours':WINDOW,'flow':flows,'primary_n':len(I),'all_qualifying_pair_records':len(F),'repeat_48h_subjects':int((nrep>=2).sum()),'repeat_48h_records':len(J),'repeat_only_records':int(nrep[nrep>=2].sum()),'repeat_48h_threeplus_subjects':int((nrep>=3).sum()),
     'primary_age_days_median_iqr':I.age_days.quantile([.25,.5,.75]).to_dict(),'primary_age7plus':int((I.age_days>=7).sum()),
     'primary_with_FT3':int(I.FT3.notna().sum()),'primary_lag_abs_hours_quantiles':I.lag_hours.abs().quantile([0,.25,.5,.75,1]).to_dict(),'primary_lag_within6h':int(I.lag_hours.abs().le(6).sum()),'primary_lag_within12h':int(I.lag_hours.abs().le(12).sum()),'primary_renal_before_thyroid':int(I.lag_hours.le(0).sum()),
     'primary_covariate_counts':{c:int(I[c].sum()) for c in ['preterm_recorded','lbw_recorded','sepsis_recorded','vaso_recent','steroid_recent','weight_record','albumin_24h','CRP_24h']},
     'duplicate_values_same_test_time':int(L[L.numeric].duplicated(['hadm_id','itemid','charttime']).sum()),
     'units':L[L.itemid.isin([5032,5041,6954,6261,6165,6183])].groupby('itemid').valueuom.unique().apply(list).to_dict(),
     'notes':[f'Thyroid TSH/FT4 exact charttime; renal CysC/Cr exact charttime; closest panels within {WINDOW} hours. Exact time is not proof of same specimen.','Weight: positive numeric charted in prior 7 days, extracted throughout neonatal hospitalization.','No exposure-outcome association tests performed.','Available RRT/ECMO labels checked, no complete therapy-status registry; systemic steroid route was not adjudicated in this feasibility-stage screen; the final cohort is defined by the later audit.']}
    I.to_json(ROOT/f'cystatin_primary_{WINDOW}h_patient_audit.json',orient='records',date_format='iso',force_ascii=False,indent=2)
    J.to_json(ROOT/f'cystatin_repeat_{WINDOW}h_patient_audit.json',orient='records',date_format='iso',force_ascii=False,indent=2)
    (ROOT/f'cystatin_feasibility_{WINDOW}h_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,default=str),encoding='utf-8')

    return {'F':F, 'E':E, 'S':S, 'L':L}


def initial_analysis(preparation_only=False):
    from pathlib import Path
    import json,hashlib,platform,warnings
    from datetime import datetime
    from zoneinfo import ZoneInfo
    import numpy as np
    import pandas as pd
    import scipy
    from scipy import stats
    import statsmodels
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests
    from statsmodels.stats.outliers_influence import variance_inflation_factor
    from patsy import build_design_matrices

    BASE=WORK
    SOURCE=BASE/'screening'
    OUT=BASE/'initial_analysis';OUT.mkdir(exist_ok=True)
    SEED=20261015

    def get_cohort(window):
     d=pd.read_json(SOURCE/f'cystatin_primary_{window}h_patient_audit.json')
     for c in ['charttime','renal_time']:d[c]=pd.to_datetime(d[c])
     return d

    primary=get_cohort(6);wide=get_cohort(24)
    assert len(primary)==primary.subject_id.nunique()==760
    assert len(wide)==wide.subject_id.nunique()==1005
    # Reconstruct all eligible 6h pairs using exactly the existing screening code,
    # stopping before any output writing. This is used only to reselect first pairs.
    scope=build_pairs(6, write_outputs=False)
    allpairs=scope['F'].copy()
    rebuilt=allpairs.sort_values(['subject_id','charttime']).drop_duplicates('subject_id').set_index('subject_id')
    assert set(rebuilt.index)==set(primary.subject_id)
    for c in ['FT4','FT3','TSH','Cr','CysC','age_days','lag_hours']:
     assert np.allclose(primary.set_index('subject_id')[c].sort_index(),rebuilt[c].sort_index(),equal_nan=True)

    E=scope['E'];S=scope['S'];L=scope['L'].copy()
    for c in ['subject_id','hadm_id']:assert primary[c].notna().all()
    assert np.isfinite(primary[['FT4','TSH','CysC','Cr','age_days','lag_hours']].to_numpy()).all()
    assert (primary[['FT4','TSH','CysC','Cr']]>0).all().all()
    assert primary.age_days.between(3,28,inclusive='left').all()
    assert primary.lag_hours.abs().le(6).all()
    assert not primary[['thyroid_dx','renal_malformation','thyroxine_before','extracorporeal_before']].any().any()

    # Verify the actual records underlying each frozen paired observation.
    thy_ids=[6095,6193,6164,6165,6183];renal_ids=[6261,5032,5041,6954]
    keys=pd.concat([
     primary[['subject_id','hadm_id','charttime']].assign(panel='thyroid'),
     primary[['subject_id','hadm_id','renal_time']].rename(columns={'renal_time':'charttime'}).assign(panel='renal')])
    selected=L.merge(keys,on=['subject_id','hadm_id','charttime'],how='inner',validate='many_to_one')
    selected=selected[((selected.panel=='thyroid')&selected.itemid.isin(thy_ids))|((selected.panel=='renal')&selected.itemid.isin(renal_ids))]
    selected['nonexact']=selected.value.fillna('').str.contains(r'[<>≤≥]|LESS|GREATER',case=False,regex=True)
    duplicate_groups=selected.groupby(['subject_id','hadm_id','charttime','itemid']).agg(n=('valuenum','size'),distinct_numeric=('valuenum','nunique'))
    discordant=duplicate_groups[(duplicate_groups.n>1)&(duplicate_groups.distinct_numeric>1)]
    bad_subjects=set(discordant.index.get_level_values('subject_id'))
    qual={
     'frozen_n':len(primary),'unique_subjects':int(primary.subject_id.nunique()),'reconstruction_matches':True,
     'gender_counts':primary.gender.value_counts().to_dict(),
     'raw_selected_record_units':{str(k):list(g.valueuom.dropna().unique()) for k,g in selected.groupby('itemid')},
     'selected_duplicate_groups':int((duplicate_groups.n>1).sum()),'discordant_duplicate_groups':len(discordant),
     'discordant_duplicate_patients':len(bad_subjects),
     'censored_records_at_selected_times':selected[selected.nonexact].groupby('itemid').size().to_dict(),
     'all_numeric_FT3_zero_or_negative':int(primary.FT3.le(0).sum()),
     'missing':primary[['FT4','FT3','TSH','Cr','CysC','gender','age_days']].isna().sum().to_dict(),
     'source_files':{f:hashlib.sha256((SOURCE/f).read_bytes()).hexdigest() for f in ['feasibility_raw.json','supplement_raw.json','cystatin_primary_6h_patient_audit.json','cystatin_primary_24h_patient_audit.json']},
     'notes':['No association-dependent exclusion applied to primary analysis.','Frozen screening uses whole-stay thyroid diagnostic labels; timing and coding remain provisional.','No measured GFR, numerical gestational age, full severity score or complete therapy registry.']
    }

    cov_labs=L[L.itemid.isin([5024,6465])&L.valuenum.gt(0)&~L.value.fillna('').str.contains(r'[<>≤≥]|LESS|GREATER',regex=True,case=False)].copy()
    weights=pd.DataFrame(S['weights_all']['rows']);weights.charttime=pd.to_datetime(weights.charttime)
    weights=weights[weights.valuenum.gt(0)&weights.valueuom.eq('kg')]
    def covariates(d):
     d=d.copy();albumins=[];ws=[]
     for r in d.itertuples():
      al=cov_labs[(cov_labs.hadm_id==r.hadm_id)&(cov_labs.charttime<=r.charttime)&(cov_labs.charttime>=r.charttime-pd.Timedelta(hours=24))].sort_values(['charttime','itemid'])
      wg=weights[(weights.hadm_id==r.hadm_id)&(weights.charttime<=r.charttime)&(weights.charttime>=r.charttime-pd.Timedelta(days=7))].sort_values('charttime')
      albumins.append(al.valuenum.iloc[-1] if len(al) else np.nan);ws.append(wg.valuenum.iloc[-1] if len(wg) else np.nan)
     d['albumin']=albumins;d['weight_kg']=ws
     return d
    primary=covariates(primary)
    qual['albumin_units']=list(cov_labs.valueuom.dropna().unique())
    qual['missing_expansion']={c:int(primary[c].isna().sum()) for c in ['albumin','weight_kg']}

    def prepare(d):
     d=d.copy()
     for c in ['FT4','FT3','TSH','CysC','Cr']:d['log_'+c]=np.log(d[c])
     for c in ['preterm_recorded','vaso_recent','steroid_recent']:d[c]=d[c].astype(int)
     d['male']=(d.gender=='M').astype(int)
     return d

    primary=prepare(primary);wide=prepare(wide)
    if preparation_only:
     return {'primary':primary, 'allpairs':allpairs}
    descriptive=primary[['FT4','FT3','TSH','CysC','Cr','age_days','lag_hours','albumin','weight_kg']].describe(percentiles=[.01,.05,.25,.5,.75,.95,.99]).T
    descriptive.to_json(OUT/'descriptive_statistics.json',orient='index',force_ascii=False,indent=2)

    def spearman_boot(x,y,n_boot=2000):
     good=np.isfinite(x)&np.isfinite(y);x=np.asarray(x)[good];y=np.asarray(y)[good]
     rho,p=stats.spearmanr(x,y);rng=np.random.default_rng(SEED)
     boot=np.empty(n_boot)
     for i in range(n_boot):
      idx=rng.integers(0,len(x),len(x));boot[i]=stats.spearmanr(x[idx],y[idx]).statistic
     low,high=np.quantile(boot,[.025,.975]);pr=stats.pearsonr(np.log(x),np.log(y))
     return {'n':len(x),'rho':float(rho),'ci_low':float(low),'ci_high':float(high),'p':float(p),'pearson_logs_r':float(pr.statistic),'pearson_logs_p':float(pr.pvalue)}

    correlations={h:spearman_boot(primary[h].values,primary.CysC.values) for h in ['FT4','FT3','TSH']}
    for h,q in zip(correlations,multipletests([v['p'] for v in correlations.values()],method='fdr_bh')[1]):correlations[h]['BH_q_three_tests']=float(q)

    knots=primary.age_days.quantile([.1,.5,.9]).to_numpy()
    age_term=f"cr(age_days, knots=({knots[1]!r},), lower_bound={knots[0]!r}, upper_bound={knots[2]!r}, constraints='center')"
    base=f'log_Cr + {age_term}'
    adjust=f'{base} + male + preterm_recorded + vaso_recent + lag_hours'
    iqr={h:float(np.log(primary[h].quantile(.75)/primary[h].quantile(.25))) for h in ['FT4','FT3','TSH']}
    models={};results=[]

    def fit(d,h='FT4',label='primary',terms=None):
     terms=adjust if terms is None else terms
     formula='log_CysC ~ log_'+h+(' + '+terms if terms else '')
     model=smf.ols(formula,data=d,missing='drop').fit(cov_type='HC3',use_t=True)
     name='log_'+h;b=float(model.params[name]);lo,hi=model.conf_int().loc[name]
     result={'label':label,'hormone':h,'n':int(model.nobs),'coefficient':b,'ci_low':float(lo),'ci_high':float(hi),'p':float(model.pvalues[name]),
     'percent_difference_IQR':float(100*np.expm1(b*iqr[h])),'percent_ci_low':float(100*np.expm1(lo*iqr[h])),'percent_ci_high':float(100*np.expm1(hi*iqr[h])),
     'IQR_from':float(primary[h].quantile(.25)),'IQR_to':float(primary[h].quantile(.75)),
     'r_squared':float(model.rsquared),'formula':formula,'covariance':'HC3 with t reference distribution'}
     models[label+'_'+h]=model;results.append(result);return model,result

    for h in ['FT4','FT3','TSH']:
     fit(primary,h,'unadjusted','')
     fit(primary,h,'creatinine_age_adjusted',base)
     fit(primary,h,'fully_adjusted')
    main_rows=[r for r in results if r['label']=='fully_adjusted']
    for r,q in zip(main_rows,multipletests([v['p'] for v in main_rows],method='fdr_bh')[1]):r['BH_q_three_hormone_models']=float(q)

    # Pre-declared sensitivity analyses. First eligible pair is reselected where needed.
    short=allpairs[allpairs.lag_hours.abs()<=1].sort_values(['subject_id','charttime']).drop_duplicates('subject_id')
    late=allpairs[(allpairs.age_days>=7)&((allpairs.renal_time-allpairs.charttime).dt.total_seconds()/86400+allpairs.age_days>=7)].sort_values(['subject_id','charttime']).drop_duplicates('subject_id')
    fit(prepare(short),'FT4','sensitivity_1h_reselected')
    fit(wide,'FT4','sensitivity_24h_reselected')
    fit(prepare(late),'FT4','sensitivity_age7d_reselected')
    fit(primary[(primary.vaso_recent==0)&(primary.steroid_recent==0)],'FT4','sensitivity_no_vaso_or_steroid',base+' + male + preterm_recorded + lag_hours')
    if bad_subjects:fit(primary[~primary.subject_id.isin(bad_subjects)],'FT4','sensitivity_no_discordant_duplicates')
    for cov in ['albumin','weight_kg']:
     sub=primary.dropna(subset=[cov]);fit(sub,'FT4','complete_subset_'+cov)
     fit(sub,'FT4','expanded_'+cov,adjust+' + '+cov)
    sub=primary.dropna(subset=['albumin','weight_kg']);fit(sub,'FT4','complete_subset_albumin_weight')
    fit(sub,'FT4','expanded_albumin_weight',adjust+' + albumin + weight_kg')

    m=models['fully_adjusted_FT4'];influence=m.get_influence();cook=influence.cooks_distance[0]
    influential=cook>4/len(primary)
    fit(primary.loc[~influential],'FT4','sensitivity_remove_Cook_gt4_over_n')
    robust=smf.rlm(m.model.formula,data=primary,M=sm.robust.norms.HuberT()).fit(cov='H3')
    b=robust.params['log_FT4'];lo,hi=robust.conf_int().loc['log_FT4']
    robust_result={'label':'Huber_RLM','n':len(primary),'coefficient':float(b),'ci_low':float(lo),'ci_high':float(hi),'p':float(robust.pvalues['log_FT4']),
     'percent_difference_IQR':float(100*np.expm1(b*iqr['FT4'])),'percent_ci_low':float(100*np.expm1(lo*iqr['FT4'])),'percent_ci_high':float(100*np.expm1(hi*iqr['FT4']))}
    vif={name:float(variance_inflation_factor(m.model.exog,i)) for i,name in enumerate(m.model.exog_names) if name!='Intercept'}

    # Within-between mixed model, explicitly restricted to the repeat-measured subgroup.
    repeat=pd.read_json(SOURCE/'cystatin_repeat_6h_patient_audit.json')
    counts=repeat.groupby('subject_id').size();repeat=prepare(repeat[repeat.subject_id.isin(counts[counts>=2].index)])
    assert len(repeat)==162 and repeat.subject_id.nunique()==77
    repeat['charttime']=pd.to_datetime(repeat.charttime)
    assert not repeat.duplicated(['subject_id','renal_time']).any()
    for v in ['FT4','Cr']:
     repeat['between_'+v]=repeat.groupby('subject_id')['log_'+v].transform('mean')
     repeat['within_'+v]=repeat['log_'+v]-repeat['between_'+v]
    mixed_formula='log_CysC ~ within_FT4 + between_FT4 + within_Cr + between_Cr + age_days + male + preterm_recorded'
    mixed_warnings=[];mixed=None
    with warnings.catch_warnings(record=True) as caught:
     warnings.simplefilter('always')
     for optimizer in ['lbfgs','powell']:
      try:
       mixed=smf.mixedlm(mixed_formula,repeat,groups=repeat.subject_id).fit(reml=True,method=optimizer,maxiter=1000,disp=False)
       if mixed.converged:break
      except Exception as e:mixed_warnings.append(optimizer+': '+str(e))
     mixed_warnings.extend(str(w.message) for w in caught)
    longitudinal={'n_subjects':77,'n_observations':162,'mixed_warnings':mixed_warnings,'formula':mixed_formula,'mixed_converged':bool(mixed.converged) if mixed else False}
    if mixed is not None:
     for term in ['within_FT4','between_FT4']:
      b=float(mixed.params[term]);lo,hi=mixed.conf_int().loc[term]
      longitudinal[term]={'beta':b,'ci_low':float(lo),'ci_high':float(hi),'p':float(mixed.pvalues[term]),'percent_IQR':float(100*np.expm1(b*iqr['FT4'])),'percent_ci_low':float(100*np.expm1(lo*iqr['FT4'])),'percent_ci_high':float(100*np.expm1(hi*iqr['FT4']))}
    changes=[]
    for subject,g in repeat.groupby('subject_id'):
     g=g.sort_values('charttime').iloc[:2];r0,r1=g.iloc[0],g.iloc[1]
     changes.append({'dCys':r1.log_CysC-r0.log_CysC,'dFT4':r1.log_FT4-r0.log_FT4,'dCr':r1.log_Cr-r0.log_Cr,'days':(r1.charttime-r0.charttime).total_seconds()/86400})
    change=pd.DataFrame(changes);cm=smf.ols('dCys ~ dFT4 + dCr + days',change).fit(cov_type='HC3',use_t=True)
    clo,chi=cm.conf_int().loc['dFT4'];longitudinal['paired_change']={'n':len(change),'beta':float(cm.params.dFT4),'ci_low':float(clo),'ci_high':float(chi),'p':float(cm.pvalues.dFT4)}

    # Exploratory nonlinearity: three-knot restricted cubic spline for log FT4.
    xkn=primary.log_FT4.quantile([.1,.5,.9]).to_numpy()
    formula_nl=f"log_CysC ~ cr(log_FT4, knots=({xkn[1]!r},), lower_bound={xkn[0]!r}, upper_bound={xkn[2]!r}, constraints='center') + {adjust}"
    nl=smf.ols(formula_nl,primary).fit(cov_type='HC3',use_t=True)
    # Nested restriction: linear contrast of three fitted values has zero curvature.
    grid=pd.DataFrame([primary.iloc[0].to_dict() for _ in range(3)])
    grid['log_FT4']=xkn
    X=np.asarray(build_design_matrices([nl.model.data.design_info],grid)[0])
    fraction=(xkn[1]-xkn[0])/(xkn[2]-xkn[0]);contrast=X[1]-(1-fraction)*X[0]-fraction*X[2]
    nltest=nl.t_test(contrast)

    qual['diagnostics']={'cook_gt4_over_n_n':int(influential.sum()),'max_cook':float(cook.max()),'max_leverage':float(influence.hat_matrix_diag.max()),'VIF':vif,'nonlinearity_p_exploratory':float(nltest.pvalue)}
    qual['primary_eligible_raw_records_verified']=len(selected)
    result={'executed_at':datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),'source_snapshot_extracted_at':E.get('executed_at'),
     'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__},
     'quality':qual,'spearman':correlations,'regressions':results,'robust_regression':robust_result,'longitudinal':longitudinal,
     'age_spline_knots_days':knots.tolist(),'log_FT4_spline_knots':xkn.tolist(),'preliminary':True}
    (OUT/'analysis_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
    modeltext=[]
    for name,model in models.items():modeltext.append(name+'\n'+model.summary().as_text())
    if mixed:modeltext.append('mixed repeated observations\n'+mixed.summary().as_text())
    (OUT/'model_summaries.txt').write_text('\n\n'.join(modeltext),encoding='utf-8')

    # A static research figure using the frozen sample and fitted models.


def prepare_covariates():
    from pathlib import Path
    import json,re,hashlib
    from datetime import datetime
    import numpy as np
    import pandas as pd
    BASE=WORK;OUT=BASE/'covariate_analysis'
    # Reuse the audited frozen-cohort preparation, stopping before any results output.
    scope=initial_analysis(preparation_only=True)
    d=scope['primary'].copy();old=d.copy()
    X=json.loads((OUT/'extended_source_snapshot.json').read_text(encoding='utf-8'))
    tables={k:pd.DataFrame(v['rows']) for k,v in X['queries'].items()}
    dx=tables['diagnoses'];dr=tables['drugs'];vi=tables['vitals'];la=tables['labs'];we=tables['weight'];sy=tables['symptoms']
    for t,c in [(dr,'startdate'),(dr,'enddate'),(vi,'charttime'),(la,'charttime'),(we,'charttime'),(we,'dob'),(sy,'recordtime')]:t[c]=pd.to_datetime(t[c])
    dr['text']=dr[['drug_name','drug_name_cn','drug_name_generic']].fillna('').agg(' '.join,axis=1)
    d['anchor']=d[['charttime','renal_time']].min(axis=1)
    qa={'frozen_n':len(d),'source_extracted_at':X['extracted_at'],'readonly':X['queries']['connection']['rows'],
     'source_sha256':{f:hashlib.sha256((OUT/f).read_bytes()).hexdigest() for f in ['metadata_audit.json','extended_source_snapshot.json']},
     'definitions':{'anchor':'Earlier of thyroid and renal panel charttimes; covariates use no later charttime.',
     'weight_kg':'Original last positive kg measurement within 7 days before thyroid panel; not birthweight.',
     'strict_weight_kg':'Last positive kg measurement within 7 days before earlier panel; not birthweight.',
     'lbw_recorded':'P07.0/P07.1 discharge diagnostic record; no record is not normal birthweight.',
     'sga_recorded':'P05.0/P05.1 diagnostic record; no record is not AGA.',
     'vitals':'Median in prior 24h, within admission. HR 30-350 bpm, RR 1-150/min, SpO2 (0,100]. Physiological-range filtering is exploratory and independent of hormone association.',
     'labs':'Latest exact positive value in prior 24h, consistent units; censored results are not exact values.',
     'drugs':'Prescription start <=anchor and end >=anchor-48h or missing; ordered medication, not verified administration; injection formulation is not an administration route.'}}
    vi['valid']=False
    for item,lo,hi in [('1001',25,45),('1002',30,350),('1003',30,350),('1004',1,150),('1006',0.01,100),('1015',1,150),('1016',1,250)]:
     vi.loc[vi.itemid.eq(item),'valid']=vi.loc[vi.itemid.eq(item),'valuenum'].between(lo,hi)
    qa['invalid_vitals_by_item']=vi.loc[~vi.valid].groupby('itemid').size().to_dict()
    qa['vital_units']=vi.groupby('itemid').valueuom.apply(lambda x:list(x.dropna().unique())).to_dict()
    qa['lab_units']=la.groupby('itemid').valueuom.apply(lambda x:list(x.dropna().unique())).to_dict()
    qa['drug_form_counts']=dr.drug_form.value_counts().to_dict()
    qa['drug_missing_enddate']=int(dr.enddate.isna().sum());qa['drug_end_before_start']=int((dr.enddate<dr.startdate).sum())
    la['exact']=la.valuenum.gt(0)&~la.value.fillna('').str.contains(r'[<>≤≥]|LESS|GREATER',case=False,regex=True)
    qa['nonexact_labs']=la.loc[~la.exact].groupby('itemid').size().to_dict()
    qa['discordant_lab_duplicate_groups']=int((la[la.exact].groupby(['hadm_id','itemid','charttime']).valuenum.nunique()>1).sum())
    med_patterns={'vaso_strict':r'dopamin|多巴胺|epinephrine|肾上腺素|dobutamin|多巴酚丁胺|milrinone|米力农',
     'dopamine_strict':r'dopamin|多巴胺','steroid_strict':r'hydrocortis|氢化可的松|dexameth|地塞米松|methylpred|甲泼|prednis|泼尼'}
    new=[]
    for r in d.itertuples():
     row={'subject_id':r.subject_id};g=dx[dx.hadm_id.eq(r.hadm_id)];codes=g.icd10_code.fillna('').str.strip()
     row['preterm_live']=int(codes.str.match(r'P07\.[23]').any())
     for key,pat in [('lbw_recorded',r'P07\.[01]'),('sga_recorded',r'P05\.[01]'),('extreme_prematurity_recorded',r'P07\.2')]:row[key]=int(codes.str.match(pat).any())
     g=dr[dr.hadm_id.eq(r.hadm_id)]
     assert (g.startdate<=r.anchor).all()
     for key,pat in med_patterns.items():row[key]=int(g.text.str.contains(pat,case=False,regex=True).any())
     recent_known=g[(g.startdate>=r.anchor-pd.Timedelta(hours=48))|g.enddate.notna()]
     uncertain=g[g.enddate.isna()&(g.startdate<r.anchor-pd.Timedelta(hours=48))]
     for key,pat in med_patterns.items():
      row[key.replace('_strict','_bounded')]=int(recent_known.text.str.contains(pat,case=False,regex=True).any())
     row['uncertain_drug_duration']=int(len(uncertain)>0)
     g=vi[vi.hadm_id.eq(r.hadm_id)&vi.valid]
     for key,ids in [('hr',['1002','1003']),('rr',['1004']),('spo2',['1006']),('temperature',['1001']),('sbp',['1016']),('dbp',['1015'])]:
      z=g[g.itemid.isin(ids)].copy()
      # Same-time HR takes precedence over pulse; equal item/time duplicates use median.
      if key=='hr':z=z.sort_values('itemid').drop_duplicates('charttime',keep='last')
      else:z=z.groupby('charttime',as_index=False).valuenum.median()
      row[key]=float(z.valuenum.median()) if len(z) else np.nan
     g=la[la.hadm_id.eq(r.hadm_id)&la.exact].sort_values(['charttime','itemid'])
     for key,ids in [('albumin_strict',[5024,6465]),('lactate',[5227]),('crp',[5626,5821])]:
      z=g[g.itemid.isin(ids)];row[key]=float(z.valuenum.iloc[-1]) if len(z) else np.nan
     g=we[we.hadm_id.eq(r.hadm_id)&we.valuenum.gt(0)&we.valueuom.eq('kg')].sort_values('charttime')
     z=g[g.charttime>=r.anchor-pd.Timedelta(days=7)]
     row['strict_weight_kg']=float(z.valuenum.iloc[-1]) if len(z) else np.nan
     row['any_prior_weight_kg']=float(g.valuenum.iloc[-1]) if len(g) else np.nan
     z=g[(g.charttime>=g.dob)&(g.charttime<g.dob+pd.Timedelta(days=1))]
     row['firstday_weight_kg']=float(z.valuenum.iloc[0]) if len(z) else np.nan
     g=sy[sy.hadm_id.eq(r.hadm_id)&(sy.recordtime<=r.anchor)&(sy.recordtime>=r.anchor-pd.Timedelta(hours=24))]
     row['symptom_documented_24h']=int(len(g)>0)
     row['oxygen_positive_documentation']=int((g.symptom_attribute.eq('+')&g.symptom_name.fillna('').str.replace('\xa0',' ').str.contains('nasal catheter to snuff oxygen',case=False)).any())
     new.append(row)
    added=pd.DataFrame(new).set_index('subject_id')
    for c in added:d[c]=d.subject_id.map(added[c])
    assert len(d)==d.subject_id.nunique()==760
    assert (d.preterm_live==d.preterm_recorded).all()
    qa['drug_flag_disagreement_vs_original']={v:int((d[v+'_strict']!=d[v+'_recent']).sum()) for v in ['vaso','steroid']}
    qa['available_counts']={c:int(d[c].notna().sum()) for c in ['albumin','weight_kg','albumin_strict','strict_weight_kg','any_prior_weight_kg','firstday_weight_kg','lactate','crp','hr','rr','spo2','temperature','sbp','dbp']}
    qa['recorded_flag_counts']={c:int(d[c].sum()) for c in ['preterm_recorded','lbw_recorded','sga_recorded','extreme_prematurity_recorded','vaso_strict','dopamine_strict','steroid_strict','symptom_documented_24h','oxygen_positive_documentation']}
    qa['bounded_drug_flags']={c:int(d[c].sum()) for c in ['vaso_bounded','dopamine_bounded','steroid_bounded','uncertain_drug_duration']}
    qa['strict_values_changed_vs_original']={c:int((d[c].notna()&d[c+'_strict'].notna()&~np.isclose(d[c],d[c+'_strict'],equal_nan=True)).sum()) for c in ['albumin']}
    qa['severity_complete_n']=int(d[['hr','rr','spo2','lactate']].notna().all(axis=1).sum())
    qa['severity_albumin_complete_n']=int(d[['hr','rr','spo2','lactate','albumin_strict']].notna().all(axis=1).sum())
    qa['repeat_audit']={'all_eligible_pairs':len(scope['allpairs']),'all_eligible_subjects':int(scope['allpairs'].subject_id.nunique()),'subjects_with_any_repeat':int((scope['allpairs'].groupby('subject_id').size()>1).sum()),'original_48h_repeated_subjects':77}
    qa['unavailable']=['No numerical gestational-age or birthweight source column/item found in pic_icu metadata.','Selected-cohort emr_symptoms contains only +/- symptoms; no numerical gestational-age or birthweight record.','No structured ventilator status, FiO2 or validated complete neonatal severity score. Sparse oxygen symptom documentation cannot establish absence of support.','Routine first-day weight is not verified birthweight; whole-stay maturity diagnoses have no onset timestamp.']
    d.to_json(OUT/'enriched_frozen_cohort_local.json',orient='records',date_format='iso',force_ascii=False,indent=2)
    (OUT/'covariate_quality.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')


def extended_analysis():
    from pathlib import Path
    from datetime import datetime
    import json,hashlib,platform
    import numpy as np,pandas as pd
    import scipy,statsmodels
    import statsmodels.formula.api as smf
    from statsmodels.stats.outliers_influence import variance_inflation_factor
    BASE=WORK;OUT=BASE/'covariate_analysis'
    d=pd.read_json(OUT/'enriched_frozen_cohort_local.json')
    old=json.loads((BASE/'initial_analysis/analysis_results.json').read_text(encoding='utf-8'))
    assert len(d)==d.subject_id.nunique()==760
    assert np.isfinite(d[['log_FT4','log_CysC','log_Cr','age_days']]).all().all()
    d['log_lactate']=np.log(d.lactate)
    knots=old['age_spline_knots_days'];iqr=np.log(d.FT4.quantile(.75)/d.FT4.quantile(.25))
    age=f"cr(age_days, knots=({knots[1]!r},), lower_bound={knots[0]!r}, upper_bound={knots[2]!r}, constraints='center')"
    B=age+' + male + preterm_recorded';C=B+' + log_Cr + lag_hours';D=C+' + vaso_strict'
    results=[];mods={}
    def fit(frame,label,terms):
     f='log_CysC ~ log_FT4'+(' + '+terms if terms else '')
     m=smf.ols(f,frame,missing='drop').fit(cov_type='HC3',use_t=True)
     b=float(m.params.log_FT4);lo,hi=m.conf_int().loc['log_FT4'];x=m.model.exog
     assert np.linalg.matrix_rank(x)==x.shape[1],label+' rank deficient'
     row={'label':label,'n':int(m.nobs),'parameters':x.shape[1],'beta':b,'beta_low':float(lo),'beta_high':float(hi),'p':float(m.pvalues.log_FT4),
     'percent':float(100*np.expm1(b*iqr)),'low':float(100*np.expm1(lo*iqr)),'high':float(100*np.expm1(hi*iqr)),
     'r2':float(m.rsquared),'formula':f,'max_vif':float(max(variance_inflation_factor(x,i) for i,name in enumerate(m.model.exog_names) if name!='Intercept'))}
     row['max_leverage']=float(m.get_influence().hat_matrix_diag.max())
     row['hc3_inference_usable']=row['max_leverage']<1-1e-8
     if not row['hc3_inference_usable']:row['inference_limitation']='Singleton vasoactive-prescription indicator yields leverage approximately1; do not interpret HC3 interval or P value. See separately labelled no-vaso subset analysis.'
     mods[label]=m;results.append(row);return m
    for name,terms in [('A_crude',''),('B_maturity',B),('C_creatinine',C),('D_temporal_drugs',D),('E_maturity_proxies',D+' + lbw_recorded + sga_recorded')]:fit(d,name,terms)
    original=fit(d,'original_primary_reproduced',C+' + vaso_recent')
    orig=next(r for r in old['regressions'] if r['label']=='fully_adjusted' and r['hormone']=='FT4')
    assert np.isclose(original.params.log_FT4,orig['coefficient'],rtol=1e-7)
    extensions={'lactate':['log_lactate'],'albumin':['albumin_strict'],'weight':['strict_weight_kg'],
     'lactate_albumin':['log_lactate','albumin_strict'],'severity':['hr','rr','spo2','log_lactate'],
     'severity_albumin':['hr','rr','spo2','log_lactate','albumin_strict']}
    skipped=[]
    for label,cols in extensions.items():
     sub=d.dropna(subset=cols)
     if len(sub)<100:skipped.append({'label':label,'n':len(sub),'reason':'Predeclared minimum100 complete subjects'});continue
     fit(sub,label+'_same_subset_base',D)
     fit(sub,label+'_extended',D+' + '+' + '.join(cols))
    fit(d,'bounded_drug_definition',C+' + vaso_bounded')
    fit(d[~d.uncertain_drug_duration.astype(bool)],'exclude_uncertain_drug_duration',C+' + vaso_bounded')
    # Integrity repair after the first extension run: singleton indicator makes HC3 undefined.
    repair={'written_at':datetime.now().astimezone().isoformat(),'status':'Post-result numerical diagnostic repair, not a pre-result specification.',
     'finding':'Severity complete subset has only1 vasoactive-prescription positive subject; the corresponding leverage is approximately1, making HC3 inference unusable.',
     'action':'Keep original model output flagged unusable. Restrict severity and severity+albumin comparisons to patients without broad prior/possibly recent vasoactive prescriptions; remove the constant medication covariate in BOTH paired models. No original primary subject removed. Report changed target subset.'}
    if not (OUT/'severity_numerical_repair.json').exists():(OUT/'severity_numerical_repair.json').write_text(json.dumps(repair,ensure_ascii=False,indent=2),encoding='utf-8')
    for name in ['severity','severity_albumin']:
     cols=extensions[name];sub=d.dropna(subset=cols);sub=sub[sub.vaso_strict.eq(0)]
     fit(sub,name+'_no_vaso_base',C)
     fit(sub,name+'_no_vaso_extended',C+' + '+' + '.join(cols))
    # Independent NumPy coefficient and HC3 covariance calculation.
    m=mods['D_temporal_drugs'];x=m.model.exog;y=m.model.endog
    inv=np.linalg.inv(x.T@x);beta=inv@x.T@y;res=y-x@beta;h=np.sum(x*(x@inv),axis=1)
    cov=inv@(x.T@((res/(1-h))[:,None]**2*x))@inv
    np.testing.assert_allclose(beta,m.params,rtol=1e-8,atol=1e-10)
    np.testing.assert_allclose(cov,m.cov_params(),rtol=1e-8,atol=1e-10)
    reduced=smf.ols('log_CysC ~ '+D,d).fit()
    delta=m.rsquared-reduced.rsquared
    rng=np.random.default_rng(SEED);betas=[];failed=0
    for i in range(2000):
     idx=rng.integers(0,len(d),len(d));xb=x[idx];yb=y[idx]
     if np.linalg.matrix_rank(xb)<x.shape[1]:failed+=1;continue
     betas.append(float(np.linalg.lstsq(xb,yb,rcond=None)[0][m.model.exog_names.index('log_FT4')]))
    ci=np.quantile(betas,[.025,.975]);bootstrap={'requested':2000,'successful':len(betas),'rank_failed':failed,'beta_low':float(ci[0]),'beta_high':float(ci[1]),'percent_low':float(100*np.expm1(ci[0]*iqr)),'percent_high':float(100*np.expm1(ci[1]*iqr)),'knots_fixed':True}

    # Missingness descriptions use standardized differences, not hypothesis tests.
    binary=['male','preterm_recorded','lbw_recorded','sga_recorded','vaso_strict','steroid_strict']
    continuous=['age_days','FT4','CysC','Cr','TSH','log_FT4','log_CysC','log_Cr']
    groups={'weight_kg':d.weight_kg.notna(),'albumin_strict':d.albumin_strict.notna(),'lactate':d.lactate.notna(),
     'severity_complete':d[extensions['severity']].notna().all(axis=1)}
    comparisons=[]
    for key,mask in groups.items():
     a=d[mask];b=d[~mask]
     for col in binary+continuous:
      av=a[col].dropna();bv=b[col].dropna();p1=av.mean();p0=bv.mean()
      sd=np.sqrt((p1*(1-p1)+p0*(1-p0))/2) if col in binary else np.sqrt((av.var(ddof=1)+bv.var(ddof=1))/2)
      smd=float((p1-p0)/sd) if sd>0 else (0. if p1==p0 else None)
      if col in binary:sa={'n':int(av.sum()),'percent':float(p1*100)};sb={'n':int(bv.sum()),'percent':float(p0*100)}
      else:sa={str(q):float(av.quantile(q)) for q in [.25,.5,.75]};sb={str(q):float(bv.quantile(q)) for q in [.25,.5,.75]}
      comparisons.append({'availability_group':key,'variable':col,'available_n':len(a),'missing_n':len(b),'available':sa,'missing':sb,'smd_available_minus_missing':smd})
    missingcols=['weight_kg','albumin_strict','lactate','hr','rr','spo2']
    patterns=d[missingcols].notna().astype(int).groupby(missingcols,dropna=False).size().reset_index(name='n')
    patterns=patterns.sort_values('n',ascending=False)
    out={'executed_at':datetime.now().astimezone().isoformat(),'exploratory_after_initial_results':True,'n':len(d),'iqr_from':float(d.FT4.quantile(.25)),'iqr_to':float(d.FT4.quantile(.75)),
     'models':results,'skipped':skipped,'bootstrap_D':bootstrap,'missingness_comparisons':comparisons,'missingness_patterns':patterns.to_dict('records'),
     'diagnostics_D':{'added_R2':float(delta),'partial_R2':float(delta/(1-reduced.rsquared)),'numpy_HC3_verified':True,'cook_gt4_over_n':int((m.get_influence().cooks_distance[0]>4/len(d)).sum()),'patients_removed_from_primary':0},
     'source_sha256':{str(f.relative_to(BASE)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [OUT/'enriched_frozen_cohort_local.json',BASE/'initial_analysis/analysis_results.json']},
     'software':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__}}
    (OUT/'strengthened_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    (OUT/'all_model_summaries.txt').write_text('\n\n'.join(name+'\n'+mod.summary().as_text() for name,mod in mods.items()),encoding='utf-8')


def admission_weight_analysis():
    from pathlib import Path
    from datetime import datetime
    import json,hashlib
    import numpy as np,pandas as pd
    import statsmodels.formula.api as smf
    BASE=WORK;OUT=BASE/'covariate_analysis'

    d=pd.read_json(OUT/'enriched_frozen_cohort_local.json')
    for c in ['charttime','renal_time','anchor']:d[c]=pd.to_datetime(d[c])
    x=json.loads((OUT/'extended_source_snapshot.json').read_text(encoding='utf-8'))['queries']
    a=pd.DataFrame(json.loads((BASE/'screening/feasibility_raw.json').read_text(encoding='utf-8'))['queries']['cohort']['rows'])
    for c in ['admittime','dob']:a[c]=pd.to_datetime(a[c])
    d=d.merge(a[['subject_id','hadm_id','admittime','dob']],on=['subject_id','hadm_id'],validate='one_to_one')
    d['admission_age_days']=(d.admittime-d.dob).dt.total_seconds()/86400
    w=pd.DataFrame(x['weight']['rows']);w.charttime=pd.to_datetime(w.charttime)
    w=w[w.valuenum.gt(0)&w.valueuom.eq('kg')&~w.value.fillna('').str.contains(r'[<>≤≥]|LESS|GREATER',regex=True,case=False)]
    discordant=int((w.groupby(['hadm_id','charttime']).valuenum.nunique()>1).sum())
    w=w.groupby(['hadm_id','charttime'],as_index=False).valuenum.median().sort_values('charttime')
    values=[]
    for r in d.itertuples():
     g=w[(w.hadm_id==r.hadm_id)&(w.charttime>=r.admittime)&(w.charttime<=r.anchor)]
     row={'subject_id':r.subject_id}
     for h in [24,48]:
      z=g[g.charttime<=r.admittime+pd.Timedelta(hours=h)]
      row[f'admission_weight_{h}h']=float(z.valuenum.iloc[0]) if len(z) else np.nan
      if h==24:
       row['first_weight_delay_hours']=(z.charttime.iloc[0]-r.admittime).total_seconds()/3600 if len(z) else np.nan
       row['first_weight_age_days']=(z.charttime.iloc[0]-r.dob).total_seconds()/86400 if len(z) else np.nan
     values.append(row)
    d=d.merge(pd.DataFrame(values),on='subject_id',validate='one_to_one')
    sub=d.dropna(subset=['admission_weight_24h'])
    r=json.loads((OUT/'strengthened_results.json').read_text(encoding='utf-8'))
    base=next(x['formula'] for x in r['models'] if x['label']=='D_temporal_drugs');iqr=np.log(r['iqr_to']/r['iqr_from'])
    models=[]
    for name,extra in [('same_subset_D',''),('same_subset_D_admission_age',' + admission_age_days'),('admission_weight_and_age',' + admission_age_days + admission_weight_24h')]:
     m=smf.ols(base+extra,sub).fit(cov_type='HC3',use_t=True)
     assert np.linalg.matrix_rank(m.model.exog)==m.model.exog.shape[1]
     b=m.params.log_FT4;lo,hi=m.conf_int().loc['log_FT4']
     models.append({'label':name,'n':int(m.nobs),'formula':base+extra,'beta':float(b),'p':float(m.pvalues.log_FT4),'percent':float(100*np.expm1(b*iqr)),'low':float(100*np.expm1(lo*iqr)),'high':float(100*np.expm1(hi*iqr))})
    binary=['male','preterm_recorded','vaso_strict'];comparisons=[]
    for col in binary+['age_days','admission_age_days','FT4','CysC','Cr']:
     av=sub[col];bv=d.loc[d.admission_weight_24h.isna(),col];p1=av.mean();p0=bv.mean()
     sd=np.sqrt((p1*(1-p1)+p0*(1-p0))/2) if col in binary else np.sqrt((av.var()+bv.var())/2)
     def desc(v):return {'n':int(v.sum()),'percent':float(v.mean()*100)} if col in binary else {str(q):float(v.quantile(q)) for q in [.25,.5,.75]}
     comparisons.append({'variable':col,'available':desc(av),'missing':desc(bv),'smd_available_minus_missing':float((p1-p0)/sd)})
    out={'executed_at':datetime.now().astimezone().isoformat(),'primary_n':len(d),'available24h':len(sub),'available48h':int(d.admission_weight_48h.notna().sum()),'discordant_weight_time_groups':discordant,
     'weight24h_range_kg':[float(sub.admission_weight_24h.min()),float(sub.admission_weight_24h.max())],
     'admission_age_days_available':{str(q):float(sub.admission_age_days.quantile(q)) for q in [0,.25,.5,.75,1]},
     'weight_delay_hours':{str(q):float(sub.first_weight_delay_hours.quantile(q)) for q in [0,.25,.5,.75,1]},
     'weight_age_days':{str(q):float(sub.first_weight_age_days.quantile(q)) for q in [0,.25,.5,.75,1]},
     'firstday_weight_count_in_24h_group':int(sub.first_weight_age_days.lt(1).sum()),'admitted_age7d_or_more':int(sub.admission_age_days.ge(7).sum()),
     'models':models,'missingness':comparisons,'source_sha256':{f:hashlib.sha256((OUT/f).read_bytes()).hexdigest() for f in ['enriched_frozen_cohort_local.json','extended_source_snapshot.json']}}
    (OUT/'admission_weight_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    d.to_json(OUT/'admission_weight_cohort_local.json',orient='records',date_format='iso',indent=2)


def audit_cohort(require_tsh=True, write_outputs=True):
    from pathlib import Path
    from datetime import datetime
    from zoneinfo import ZoneInfo
    import json,re,hashlib
    import numpy as np,pandas as pd
    B=WORK;O=B/'primary_analysis'
    X=json.loads((O/'reaudit_source_local.json').read_text(encoding='utf-8'));Q=X['queries']
    A=pd.DataFrame(Q['cohort']['rows']);L=pd.DataFrame(Q['labs']['rows']);DX=pd.DataFrame(Q['diagnoses']['rows']);RX=pd.DataFrame(Q['thyroid_drugs_expanded']['rows']);I=pd.DataFrame(Q['icustays']['rows'])
    for df,cols in [(A,['admittime','dischtime','dob']),(L,['charttime']),(RX,['startdate','enddate']),(I,['intime','outtime'])]:
     for c in cols:df[c]=pd.to_datetime(df[c])
    assert not A.duplicated(['subject_id','hadm_id']).any()
    A=A.sort_values(['subject_id','admittime','hadm_id']).drop_duplicates('subject_id')
    valid=L.valuenum.gt(0)&~L.value.fillna('').str.contains(r'[<>≤≥]|LESS|GREATER',case=False,regex=True)
    mapping={6095:'TSH',6193:'TSH',6164:'FT3',6165:'FT4',6183:'FT4',5032:'Cr',5041:'Cr',6954:'Cr',6261:'CysC'}
    V=L[valid&L.itemid.isin(mapping)].copy();V['test']=V.itemid.map(mapping)
    V=V.groupby(['subject_id','hadm_id','charttime','test']).valuenum.median().unstack('test').reset_index()
    V=V.merge(A[['subject_id','hadm_id','admittime','dischtime','dob','gender']],on=['subject_id','hadm_id'],validate='many_to_one')
    V=V[(V.charttime>=V.admittime)&(V.charttime<=V.dischtime)]
    V['age_days']=(V.charttime-V.dob).dt.total_seconds()/86400
    T=V[V.age_days.between(3,28,inclusive='left')].dropna(subset=['TSH','FT4'] if require_tsh else ['FT4'])
    K=V[V.age_days.between(3,28,inclusive='left')].dropna(subset=['CysC','Cr']);KG={h:g for h,g in K.groupby('hadm_id')}
    DX['code']=DX.icd10_code.fillna('').str.strip()
    thy_dx=DX[DX.code.str.match(r'^E0[0-7](\.|$)')]
    renal_dx=DX[DX.code.str.match(r'^Q6[0-4](\.|$)')]
    text_dx=DX[DX.title_cn.fillna('').str.contains(r'甲状腺|甲状腺素')]
    rxtext=RX[['drug_name','drug_name_cn','drug_name_generic']].fillna('').agg(' '.join,axis=1)
    RX['replacement']=rxtext.str.contains(r'levothyrox|liothyron|thyroid|左甲状腺素|甲状腺片',case=False,regex=True)
    RX['antithyroid']=rxtext.str.contains(r'methimazole|thiamazole|carbimazole|propylthiouracil|甲巯咪唑|丙硫氧嘧啶|卡比马唑',case=False,regex=True)
    sx=pd.DataFrame(Q['rrt_ecmo_surgery']['rows'])
    if len(sx):sx['surgery_begin_time']=pd.to_datetime(sx.surgery_begin_time)
    records=[]
    for r in T.itertuples():
     g=KG.get(r.hadm_id)
     if g is None:continue
     g=g.assign(lag=(g.charttime-r.charttime).dt.total_seconds()/3600)
     g=g[g.lag.abs()<=6].assign(abs_lag=lambda x:x.lag.abs()).sort_values(['abs_lag','charttime'])
     if not len(g):continue
     v=g.iloc[0];anchor=min(v.charttime,r.charttime);end=max(v.charttime,r.charttime)
     drug=RX[(RX.hadm_id==r.hadm_id)&(RX.startdate<=end)]
     row={c:getattr(r,c) for c in ['subject_id','hadm_id','charttime','FT4','FT3','TSH','age_days','gender']}
     row.update({'renal_time':v.charttime,'CysC':v.CysC,'Cr':v.Cr,'lag_hours':v.lag,
     'thyroid_dx':r.hadm_id in set(thy_dx.hadm_id),'renal_malformation':r.hadm_id in set(renal_dx.hadm_id),
     'thyroid_rx_before_later_panel':bool((drug.replacement|drug.antithyroid).any()),
     'thyroid_rx_before_thyroid_panel':bool(((drug.replacement|drug.antithyroid)&(drug.startdate<=r.charttime)).any()),
     'thyroid_rx_between_panels':bool(((drug.replacement|drug.antithyroid)&(drug.startdate>anchor)).any()),
     'extracorporeal':bool(((sx.hadm_id==r.hadm_id)&(sx.surgery_begin_time<=end)).any()) if len(sx) else False})
     records.append(row)
    F=pd.DataFrame(records).sort_values(['subject_id','charttime']);all_f=F.copy()
    flow=[{'stage':'Neonatal hospitalizations with recorded NICU stay','n':len(Q['cohort']['rows']),'unit':'admissions'},
     {'stage':'First eligible neonatal NICU admission per subject','n':len(A),'unit':'subjects'},
     {'stage':'Exact positive TSH and FT4 at same timestamp, age72h to<28d','n':T.subject_id.nunique(),'unit':'subjects'},
     {'stage':'CysC and creatinine paired within6h; both age72h to<28d','n':F.subject_id.nunique(),'unit':'subjects'}]
    exclusions=[]
    for flag in ['thyroid_dx','renal_malformation','thyroid_rx_before_later_panel','extracorporeal']:
     before=set(F.subject_id);F=F[~F[flag]];removed=before-set(F.subject_id)
     exclusions.append({'reason':flag,'excluded_subjects':len(removed),'remaining_subjects':F.subject_id.nunique()})
     flow.append({'stage':'After '+flag,'n':F.subject_id.nunique(),'unit':'subjects'})
    if not write_outputs:
     return {'A':A,'L':L,'V':V,'T':T,'K':K,'F':F,'all_f':all_f,'flow':flow,'exclusions':exclusions,'DX':DX,'mapping':mapping,'Q':Q}
    C=F.drop_duplicates('subject_id');prior=pd.read_json(B/'screening/cystatin_primary_6h_patient_audit.json')
    assert set(C.subject_id)==set(prior.subject_id),'Cohort differs; investigate before model lock.'
    for col in ['FT4','FT3','TSH','Cr','CysC','lag_hours','age_days']:
     np.testing.assert_allclose(C.set_index('subject_id')[col].sort_index(),prior.set_index('subject_id')[col].sort_index(),rtol=1e-7,equal_nan=True)
    for col in ['charttime','renal_time']:
     assert (C.set_index('subject_id')[col].sort_index()==pd.to_datetime(prior.set_index('subject_id')[col]).sort_index()).all()
    candidate_h=set(all_f.hadm_id);code_rows=[]
    for label,g in [('thyroid',thy_dx),('urinary',renal_dx)]:
     g=g[g.hadm_id.isin(candidate_h)]
     for keys,z in g.groupby(['icd10_code_cn','code','title_cn']):code_rows.append({'domain':label,'local_code':keys[0],'icd10_code':keys[1],'title_cn':keys[2],'subjects':int(z.subject_id.nunique())})
    priorad=pd.DataFrame(Q['all_admissions_of_subjects']['rows']);priorad.admittime=pd.to_datetime(priorad.admittime)
    earliest=priorad.groupby('subject_id').admittime.min();first=A.set_index('subject_id').admittime
    notfirst=set(first[first>first.index.map(earliest)].index)
    in_icu=[]
    for r in C.itertuples():
     g=I[(I.hadm_id==r.hadm_id)&((I.first_careunit=='NICU')|(I.last_careunit=='NICU'))]
     in_icu.append({'subject_id':r.subject_id,'both_in_recorded_nicu_interval':bool(((g.intime<=min(r.charttime,r.renal_time))&(g.outtime>=max(r.charttime,r.renal_time))).any())})
    membership=pd.DataFrame(in_icu);C=C.merge(membership,on='subject_id',validate='one_to_one')
    out={'executed_at':datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),'fresh_readonly':Q['connection']['rows'],
     'final_n':len(C),'same_subjects_values_times_as_frozen760':True,'flow':flow,'sequential_exclusions':exclusions,'exclusion_codes_in_paired_candidates':code_rows,
     'code_text_hadm_disagreement':sorted(set(thy_dx.hadm_id)^set(text_dx.hadm_id)),
     'thyroid_med_names':RX[['drug_name','drug_name_cn']].drop_duplicates().to_dict('records'),'antithyroid_records':int(RX.antithyroid.sum()),
     'pairs_with_thyroid_med_start_between_panels':int(all_f.thyroid_rx_between_panels.sum()),
     'final_both_panels_in_recorded_nicu_interval':int(C.both_in_recorded_nicu_interval.sum()),
     'first_eligible_not_first_ever_in_source':{'all_neonatal_subjects':len(notfirst),'in_final_cohort':len(set(C.subject_id)&notfirst)},
     'dictionary_rrt_items':len(Q['rrt_ecmo_dictionary']['rows']),'rrt_surgery_records_neonatal_source':len(sx),
     'diagnosis_unmatched_dictionary_records':int(DX.icd10_code.isna().sum()),
     'limitations':['Whole-admission diagnoses lack onset timestamps; do not call them known pre-existing conditions.','Q60-Q64 excludes urinary congenital malformations including urachal anomalies; not all are renal parenchymal disease.','NICU first/last unit labels do not reconstruct all transfers.','No complete RRT registry or verified drug administration log.','No inferred normal thyroid classification from absence of diagnosis.'],
     'sha256_source':hashlib.sha256((O/'reaudit_source_local.json').read_bytes()).hexdigest()}
    (O/'cohort_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2,default=int),encoding='utf-8')
    C.to_json(O/'final_cohort_verified_local.json',orient='records',date_format='iso',force_ascii=False,indent=2)
    all_f.to_json(O/'paired_candidates_local.json',orient='records',date_format='iso',force_ascii=False,indent=2)

    return {'A':A,'L':L,'V':V,'T':T,'K':K,'F':F,'all_f':all_f,'flow':flow,'exclusions':exclusions,'DX':DX,'mapping':mapping,'Q':Q}


def primary_analysis():
    from pathlib import Path
    from datetime import datetime
    from zoneinfo import ZoneInfo
    import json,hashlib,platform,warnings
    import numpy as np,pandas as pd
    import scipy.stats as st
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests
    B=WORK;O=B/'primary_analysis'
    plan={'age_knots':AGE_KNOTS}
    d=pd.read_json(B/'covariate_analysis/admission_weight_cohort_local.json')
    verified=pd.read_json(O/'final_cohort_verified_local.json');d=d.merge(verified[['subject_id','both_in_recorded_nicu_interval']],on='subject_id',validate='one_to_one')
    assert len(d)==d.subject_id.nunique()==760
    for c in ['charttime','renal_time']:d[c]=pd.to_datetime(d[c])
    raw=json.loads((O/'reaudit_source_local.json').read_text(encoding='utf-8'));rx=pd.DataFrame(raw['queries']['drugs']['rows'])
    for c in ['startdate','enddate']:rx[c]=pd.to_datetime(rx[c])
    rx['text']=rx[['drug_name','drug_name_cn','drug_name_generic']].fillna('').agg(' '.join,axis=1);rg={h:g for h,g in rx.groupby('hadm_id')}
    def prepare(a):
     a=a.copy()
     for c in ['charttime','renal_time']:a[c]=pd.to_datetime(a[c])
     for c in ['FT4','FT3','TSH','Cr','CysC']:a['log_'+c]=np.log(a[c])
     a['male']=(a.gender=='M').astype(int);a.preterm_recorded=a.preterm_recorded.astype(int)
     flags=[]
     for r in a.itertuples():
      anchor=min(r.charttime,r.renal_time);g=rg.get(r.hadm_id,rx.iloc[:0]);g=g[g.startdate<=anchor]
      bounded=g[(g.startdate>=anchor-pd.Timedelta(hours=48))|(g.enddate>=anchor-pd.Timedelta(hours=48))]
      broad=g[(g.enddate>=anchor-pd.Timedelta(hours=48))|g.enddate.isna()]
      pat=r'dopamin|多巴胺|epinephrine|肾上腺素|dobutamin|多巴酚丁胺|milrinone|米力农'
      steroid=r'hydrocortis|氢化可的松|dexameth|地塞米松|methylpred|甲泼|prednis|泼尼'
      flags.append({'vaso_bounded':int(bounded.text.str.contains(pat,case=False).any()),'vaso_broad':int(broad.text.str.contains(pat,case=False).any()),'steroid_broad':int(broad.text.str.contains(steroid,case=False).any())})
     for c in flags[0]:a[c]=[f[c] for f in flags]
     return a
    d=prepare(d);d['log_lactate']=np.log(d.lactate)
    kn=plan['age_knots'];age=f"cr(age_days, knots=({kn[1]},), lower_bound={kn[0]}, upper_bound={kn[2]}, constraints='center')"
    maturity=age+' + male + preterm_recorded';noCr=maturity+' + lag_hours';withCr=noCr+' + log_Cr';adj=withCr+' + vaso_bounded'
    models=[];fits={};iqrs={h:np.log(d[h].quantile(.75)/d[h].quantile(.25)) for h in ['FT4','FT3','TSH']}
    def fit(a,label,terms=adj,h='FT4'):
     m=smf.ols('log_CysC ~ log_'+h+(' + '+terms if terms else ''),a).fit(cov_type='HC3',use_t=True);term='log_'+h;b=m.params[term];lo,hi=m.conf_int().loc[term]
     assert np.linalg.matrix_rank(m.model.exog)==m.model.exog.shape[1],label
     lev=float(m.get_influence().hat_matrix_diag.max());assert lev<1-1e-8,label
     row={'label':label,'hormone':h,'n':int(m.nobs),'beta':float(b),'beta_low':float(lo),'beta_high':float(hi),'p':float(m.pvalues[term]),'percent':float(100*np.expm1(b*iqrs[h])),'low':float(100*np.expm1(lo*iqrs[h])),'high':float(100*np.expm1(hi*iqrs[h])),'r2':float(m.rsquared),'formula':m.model.formula,'max_leverage':lev}
     models.append(row);fits[label]=m;return m
    for label,terms in [('crude',''),('maturity',maturity),('maturity_lag',noCr),('plus_creatinine',withCr),('final_primary',adj)]:fit(d,label,terms)
    for h in ['FT3','TSH']:fit(d,'secondary_'+h,adj,h)
    fit(d,'historical_primary',withCr+' + vaso_recent');fit(d,'broad_drug_definition',withCr+' + vaso_broad')
    fit(d[d.uncertain_drug_duration==0],'no_uncertain_drug_duration')
    fit(d[d.both_in_recorded_nicu_interval],'nicu_interval')
    fit(d[(d.vaso_broad==0)&(d.steroid_broad==0)],'no_broad_vaso_steroid',withCr)
    scope=build_pairs(6, write_outputs=False)
    ap=scope['F'];short=ap[ap.lag_hours.abs()<=1].sort_values(['subject_id','charttime']).drop_duplicates('subject_id')
    late=ap[(ap.age_days>=7)&(ap.age_days+ap.lag_hours/24>=7)].sort_values(['subject_id','charttime']).drop_duplicates('subject_id')
    fit(prepare(short),'pairing_1h');fit(prepare(late),'age7d');fit(prepare(pd.read_json(B/'screening/cystatin_primary_24h_patient_audit.json')),'pairing_24h')
    extensions={'albumin':['albumin_strict'],'lactate':['log_lactate'],'albumin_lactate':['albumin_strict','log_lactate'],'recent_weight':['strict_weight_kg']}
    for name,cols in extensions.items():
     sub=d.dropna(subset=cols);fit(sub,name+'_base');fit(sub,name+'_extended',adj+' + '+' + '.join(cols))
    sub=d.dropna(subset=['admission_weight_24h']);fit(sub,'admission_weight_base',adj+' + admission_age_days');fit(sub,'admission_weight_extended',adj+' + admission_age_days + admission_weight_24h')
    for name,cols in [('severity',['hr','rr','spo2','log_lactate']),('severity_albumin',['hr','rr','spo2','log_lactate','albumin_strict'])]:
     sub=d.dropna(subset=cols);sub=sub[sub.vaso_broad==0]
     fit(sub,name+'_base',withCr);fit(sub,name+'_extended',withCr+' + '+' + '.join(cols))
    ck=d.log_Cr.quantile([.1,.5,.9]).to_list();crterm=f"cr(log_Cr,knots=({ck[1]},),lower_bound={ck[0]},upper_bound={ck[2]},constraints='center')"
    fit(d,'nonlinear_creatinine',noCr+' + vaso_bounded + '+crterm)
    rlm=smf.rlm('log_CysC ~ log_FT4 + '+adj,d,M=sm.robust.norms.HuberT()).fit(cov='H3');b=rlm.params.log_FT4;lo,hi=rlm.conf_int().loc['log_FT4']
    robust={'label':'Huber_RLM','n':760,'percent':float(100*np.expm1(b*iqrs['FT4'])),'low':float(100*np.expm1(lo*iqrs['FT4'])),'high':float(100*np.expm1(hi*iqrs['FT4'])),'p':float(rlm.pvalues.log_FT4)}
    # Algebraic direction-reversal audit. Same covariates in reduced/full models.
    z=smf.ols('log_Cr ~ log_FT4 + '+noCr,d).fit();m0=fits['maturity_lag'];m1=fits['plus_creatinine']
    decomp={'reduced_beta':float(m0.params.log_FT4),'full_beta':float(m1.params.log_FT4),'Cr_outcome_coefficient':float(m1.params.log_Cr),'FT4_Cr_conditional_coefficient':float(z.params.log_FT4),'omitted_Cr_component':float(m1.params.log_Cr*z.params.log_FT4)}
    assert np.isclose(decomp['reduced_beta'],decomp['full_beta']+decomp['omitted_Cr_component'])
    yr=smf.ols('log_CysC ~ '+adj,d).fit().resid;xr=smf.ols('log_FT4 ~ '+adj,d).fit().resid
    fwl=float(xr@yr/(xr@xr));assert np.isclose(fwl,fits['final_primary'].params.log_FT4)
    decomp['FWL_slope']=fwl;decomp['FWL_verified']=True
    decomp['spearman']={a+'__'+b:dict(zip(['rho','p'],map(float,st.spearmanr(d[a],d[b])))) for a,b in [('FT4','CysC'),('FT4','Cr'),('Cr','CysC'),('age_days','Cr')]}
    pd.DataFrame({'x_logFT4_residual':xr,'y_logCysC_residual':yr}).to_json(O/'residuals_local.json',orient='records')
    # Bootstrap uncertainty with all rows retained and fixed spline knots.
    m=fits['final_primary'];xx=m.model.exog;yy=m.model.endog;j=m.model.exog_names.index('log_FT4');rng=np.random.default_rng(20261015);boot=[]
    for _ in range(2000):
     ix=rng.integers(0,760,760);assert np.linalg.matrix_rank(xx[ix])==xx.shape[1];boot.append(np.linalg.lstsq(xx[ix],yy[ix],rcond=None)[0][j])
    bl,bh=np.quantile(boot,[.025,.975]);bootstrap={'resamples':2000,'low':float(100*np.expm1(bl*iqrs['FT4'])),'high':float(100*np.expm1(bh*iqrs['FT4']))}
    from statsmodels.stats.outliers_influence import variance_inflation_factor
    inv=np.linalg.inv(xx.T@xx);bb=inv@xx.T@yy;hh=np.sum(xx*(xx@inv),axis=1);res=yy-xx@bb;cov=inv@(xx.T@((res/(1-hh))[:,None]**2*xx))@inv
    np.testing.assert_allclose(cov,m.cov_params(),rtol=1e-7,atol=1e-10)
    diag={'HC3_independently_verified':True,'max_vif':float(max(variance_inflation_factor(xx,i) for i in range(1,xx.shape[1]))),'cook_gt4_over_n':int((m.get_influence().cooks_distance[0]>4/760).sum()),'removed_primary':0,'added_R2':float(m.rsquared-smf.ols('log_CysC ~ '+adj,d).fit().rsquared)}
    for row,q in zip([r for r in models if r['label'] in ['final_primary','secondary_FT3','secondary_TSH']],multipletests([r['p'] for r in models if r['label'] in ['final_primary','secondary_FT3','secondary_TSH']],method='fdr_bh')[1]):row['BH_q_three_hormones']=float(q)
    baseline=[]
    cont=[('age_days','索引甲功检测日龄','天'),('admission_age_days','入院日龄','天'),('FT4','FT4','pmol/L'),('FT3','FT3精确值','pmol/L'),('TSH','TSH','mIU/L'),('CysC','胱抑素C','mg/L'),('Cr','肌酐','μmol/L'),('albumin_strict','检测前24小时白蛋白','g/L'),('lactate','检测前24小时乳酸','mmol/L'),('strict_weight_kg','检测前7天最近体重','kg'),('admission_weight_24h','入院24小时首次体重','kg')]
    for c,label,u in cont:
     v=d[c].dropna();qs=v.quantile([.25,.5,.75]).tolist();baseline.append({'variable':c,'label':label,'unit':u,'available':len(v),'missing':760-len(v),'median':qs[1],'q1':qs[0],'q3':qs[2],'display':f'{qs[1]:.2f} [{qs[0]:.2f}, {qs[2]:.2f}]'})
    for c,label in [('male','男性'),('preterm_recorded','有早产诊断记录'),('lbw_recorded','有低出生体重诊断记录'),('sga_recorded','有小于胎龄诊断记录'),('vaso_bounded','有时间可判定的近期血管活性药处方'),('steroid_bounded','有时间可判定的近期激素处方'),('both_in_recorded_nicu_interval','两组检测均处记录的NICU相关ICU时段')]:
     n=int(d[c].sum());baseline.append({'variable':c,'label':label,'unit':'','available':760,'missing':0,'count':n,'percent':100*n/760,'display':f'{n} ({100*n/760:.1f}%)'})
    older=json.loads((B/'initial_analysis/analysis_results.json').read_text(encoding='utf-8'))
    out={'executed_at':datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),'lock_version':'1.0','models':models,'robust':robust,'direction_reversal':decomp,'bootstrap_primary':bootstrap,'diagnostics':diag,'baseline':baseline,'longitudinal_historical':older['longitudinal'],'qa':{'cohort760':True,'primary_vaso_bounded':int(d.vaso_bounded.sum()),'primary_in_nicu_interval':int(d.both_in_recorded_nicu_interval.sum())},'software':{'python':platform.python_version(),'numpy':np.__version__},'analysis_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (O/'primary_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
    (O/'final_model_summaries.txt').write_text('\n\n'.join(k+'\n'+v.summary().as_text() for k,v in fits.items()),encoding='utf-8')
    d.to_json(O/'final_analysis_dataset_local.json',orient='records',date_format='iso',force_ascii=False,indent=2)

    return {'fits':fits,'models':models,'d':d,'iqrs':iqrs,'rlm':rlm}


def supplementary_analysis():
    from pathlib import Path
    import json, warnings, hashlib, platform
    from datetime import datetime
    from zoneinfo import ZoneInfo
    import numpy as np
    import pandas as pd
    import scipy, scipy.stats as st
    import statsmodels, statsmodels.formula.api as smf
    from patsy import build_design_matrices
    B=WORK; O=WORK/'supplementary_analysis'
    for f in ['evidence','figures','documents','reports','qa','code']: (O/f).mkdir(exist_ok=True)
    def read(p): return json.loads((B/p).read_text(encoding='utf-8'))
    def save(name,x): (O/'evidence'/name).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda z: z.item() if hasattr(z,'item') else str(z)),encoding='utf-8')
    d=pd.read_json(B/'primary_analysis/final_analysis_dataset_local.json')
    frozen=read('primary_analysis/primary_results.json'); row=next(x for x in frozen['models'] if x['label']=='final_primary')
    primary=smf.ols(row['formula'],d).fit(cov_type='HC3',use_t=True)
    assert len(d)==d.subject_id.nunique()==760
    assert np.isclose(primary.params.log_FT4,row['beta'],rtol=1e-9)
    q25,q50,q75=d.FT4.quantile([.25,.5,.75]); contrast=np.log(q75/q25)
    rep=pd.read_json(B/'screening/cystatin_repeat_6h_patient_audit.json')
    counts=rep.groupby('subject_id').size();rep=rep[rep.subject_id.isin(counts[counts>=2].index)].copy()
    assert len(rep)==162 and rep.subject_id.nunique()==77
    rep.charttime=pd.to_datetime(rep.charttime)
    for v in ['FT4','Cr','CysC']:rep['log_'+v]=np.log(rep[v])
    rep['male']=(rep.gender=='M').astype(int);rep.preterm_recorded=rep.preterm_recorded.astype(int)
    for v in ['FT4','Cr']:
     rep['between_'+v]=rep.groupby('subject_id')['log_'+v].transform('mean');rep['within_'+v]=rep['log_'+v]-rep['between_'+v]
    formula='log_CysC ~ within_FT4 + between_FT4 + within_Cr + between_Cr + age_days + male + preterm_recorded'
    mixed=[];summaries=[]
    for method in ['lbfgs','powell','bfgs']:
     out={'optimizer':method}
     with warnings.catch_warnings(record=True) as caught:
      warnings.simplefilter('always')
      try:
       m=smf.mixedlm(formula,rep,groups=rep.subject_id).fit(reml=True,method=method,maxiter=1000,disp=False)
       h,sing=m.model.hessian(m.params_object); cp=np.asarray(m.cov_params())
       out.update(converged=bool(m.converged),random_intercept_variance=float(m.cov_re.iloc[0,0]),residual_variance=float(m.scale),log_likelihood=float(m.llf),hessian_singular=bool(sing),negative_hessian_eigenvalues=np.linalg.eigvalsh(-h).tolist(),covariance_eigenvalues=np.linalg.eigvalsh(cp).tolist(),covariance_condition_number=float(np.linalg.cond(cp)),terms={})
       for t in ['within_FT4','between_FT4']:
        lo,hi=m.conf_int().loc[t];out['terms'][t]={'beta':float(m.params[t]),'se':float(m.bse[t]),'p':float(m.pvalues[t]),'percent':float(100*np.expm1(m.params[t]*contrast)),'low':float(100*np.expm1(lo*contrast)),'high':float(100*np.expm1(hi*contrast))}
       summaries.append(method+'\n'+m.summary().as_text())
      except Exception as e:out['error']=str(e)
      out['warnings']=[str(w.message) for w in caught]
     mixed.append(out)
    changes=[]
    for _,g in rep.groupby('subject_id'):
     a,b=g.sort_values('charttime').iloc[:2].itertuples();changes.append(dict(dCys=b.log_CysC-a.log_CysC,dFT4=b.log_FT4-a.log_FT4,dCr=b.log_Cr-a.log_Cr,days=(b.charttime-a.charttime).total_seconds()/86400))
    change=pd.DataFrame(changes);cm=smf.ols('dCys ~ dFT4 + dCr + days',change).fit(cov_type='HC3',use_t=True)
    lo,hi=cm.conf_int().loc['dFT4'];changeout={'n':77,'formula':cm.model.formula,'beta':float(cm.params.dFT4),'se':float(cm.bse.dFT4),'p':float(cm.pvalues.dFT4),'percent':float(100*np.expm1(cm.params.dFT4*contrast)),'low':float(100*np.expm1(lo*contrast)),'high':float(100*np.expm1(hi*contrast)),'rank':int(np.linalg.matrix_rank(cm.model.exog)),'n_parameters':cm.model.exog.shape[1],'max_leverage':float(cm.get_influence().hat_matrix_diag.max()),'covariance_min_eigenvalue':float(np.linalg.eigvalsh(cm.cov_params()).min())}
    summaries.append('First two changes\n'+cm.summary().as_text())
    save('longitudinal_diagnostics.json',{'formula':formula,'n_subjects':77,'n_observations':162,'optimizers':mixed,'paired_change':changeout})
    # Three-knot natural spline: two exposure df, including one nonlinear df.
    kn=d.log_FT4.quantile([.1,.5,.9]).to_numpy()
    term=f"cr(log_FT4, knots=({kn[1]!r},), lower_bound={kn[0]!r}, upper_bound={kn[2]!r}, constraints='center')"
    nlform=row['formula'].replace('log_FT4',term,1);nl=smf.ols(nlform,d).fit(cov_type='HC3',use_t=True)
    def design(values):
     g=pd.DataFrame([d.iloc[0].to_dict() for _ in values]);g['log_FT4']=values
     return np.asarray(build_design_matrices([nl.model.data.design_info],g)[0])
    x=design(kn);fraction=(kn[1]-kn[0])/(kn[2]-kn[0]);curvature=x[1]-(1-fraction)*x[0]-fraction*x[2]
    test=nl.t_test(curvature);delta=design(np.log([q25,q75]));iqtest=nl.t_test(delta[1]-delta[0]);ci=np.asarray(iqtest.conf_int()).ravel()
    limits=d.FT4.quantile([.05,.95]).to_numpy();grid=np.geomspace(*limits,201);dx=design(np.log(grid))-design([np.log(q50)])[0]
    estimate=dx@nl.params;se=np.sqrt(np.maximum(0,np.einsum('ij,jk,ik->i',dx,nl.cov_params(),dx)));crit=st.t.ppf(.975,nl.df_resid)
    pd.DataFrame({'FT4':grid,'percent':100*np.expm1(estimate),'low':100*np.expm1(estimate-crit*se),'high':100*np.expm1(estimate+crit*se),'linear_percent':100*np.expm1(primary.params.log_FT4*np.log(grid/q50))}).to_csv(O/'evidence/FT4_spline_curve.csv',index=False)
    hist,edges=np.histogram(d.FT4,bins=np.linspace(d.FT4.min(),d.FT4.max(),21));pd.DataFrame({'lower':edges[:-1],'upper':edges[1:],'count':hist}).to_csv(O/'evidence/FT4_histogram.csv',index=False)
    nlout={'formula':nlform,'n':760,'lnFT4_knots':kn.tolist(),'FT4_knots':np.exp(kn).tolist(),'nonlinearity_p':float(np.asarray(test.pvalue).item()),'contrast_Q25_Q75':[float(q25),float(q75)],'IQR_percent':float(100*np.expm1(np.asarray(iqtest.effect).item())),'IQR_low':float(100*np.expm1(ci[0])),'IQR_high':float(100*np.expm1(ci[1])),'display_5_95_percentiles':limits.tolist(),'below_display':int((d.FT4<limits[0]).sum()),'above_display':int((d.FT4>limits[1]).sum()),'below_lower_knot':int((d.log_FT4<kn[0]).sum()),'above_upper_knot':int((d.log_FT4>kn[2]).sum()),'rank':int(np.linalg.matrix_rank(nl.model.exog)),'parameters':nl.model.exog.shape[1],'max_leverage':float(nl.get_influence().hat_matrix_diag.max()),'covariance_min_eigenvalue':float(np.linalg.eigvalsh(nl.cov_params()).min()),'historical_legacy_nonlinearity_p':read('initial_analysis/analysis_results.json')['quality']['diagnostics']['nonlinearity_p_exploratory'],'interpretation':'Failure to reject the single nonlinear component does not prove linearity; tail uncertainty remains.'}
    save('FT4_functional_form.json',nlout);summaries.append('Final-adjustment FT4 spline\n'+nl.summary().as_text())
    (O/'evidence/review_model_summaries.txt').write_text('\n\n'.join(summaries),encoding='utf-8')
    save('analysis_receipt.json',{'executed_at':datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),'analysis_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'frozen_main_beta_reproduced':float(primary.params.log_FT4),'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__}})


def thyroid_selection_analysis():
    from pathlib import Path
    import json, hashlib
    import numpy as np,pandas as pd
    import statsmodels.formula.api as smf
    B=WORK;O=WORK/'supplementary_analysis'
    sc=audit_cohort(write_outputs=False)
    main=sc['F'].drop_duplicates('subject_id').copy();assert len(main)==760
    A,L,V,T,K=sc['A'],sc['L'],sc['V'],sc['T'],sc['K']
    raw=L.merge(A[['subject_id','hadm_id','dob','admittime','dischtime']],on=['subject_id','hadm_id'],validate='many_to_one')
    raw=raw[(raw.charttime>=raw.admittime)&(raw.charttime<=raw.dischtime)].copy()
    raw['test']=raw.itemid.map(sc['mapping']);raw['age']=(raw.charttime-raw.dob).dt.total_seconds()/86400
    remaining=set(A.subject_id);fail=[]
    def stage(label,eligible):
     nonlocal remaining
     gone=remaining-set(eligible);fail.append({'reason':label,'excluded':len(gone),'remaining':len(remaining)-len(gone)});remaining-=gone
    both=set(raw.loc[raw.test=='TSH','subject_id']) & set(raw.loc[raw.test=='FT4','subject_id'])
    stage('One or both thyroid tests absent throughout eligible hospitalization',both)
    rt=raw[raw.test.isin(['TSH','FT4'])].groupby(['subject_id','hadm_id','charttime']).test.nunique();same=rt[rt==2].reset_index()
    stage('Both tests recorded but never at the same timestamp',same.subject_id)
    validsame=V.dropna(subset=['TSH','FT4']);stage('Same-time thyroid records but no exact positive pair',validsame.subject_id)
    stage('Exact positive thyroid pair exists only outside age 72 h to <28 d',T.subject_id)
    assert len(remaining)==1649 and sum(x['excluded'] for x in fail)==1007
    # Hierarchical renal failures among eligible thyroid-panel patients.
    renal_fail=[];rr=remaining.copy()
    for label,ids in [('No exact positive same-time CysC/creatinine panel during hospitalization',V.dropna(subset=['CysC','Cr']).subject_id),('Renal panel exists but none at age 72 h to <28 d',K.subject_id),('Both age-eligible panels exist but none within 6 h',sc['all_f'].subject_id)]:
     excluded=rr-set(ids);rr-=excluded;renal_fail.append({'reason':label,'excluded':len(excluded),'remaining':len(rr)})
    assert len(rr)==767 and sum(x['excluded'] for x in renal_fail)==882
    ns=audit_cohort(require_tsh=False, write_outputs=False)
    expanded=ns['F'].drop_duplicates('subject_id').copy();extra=set(expanded.subject_id)-set(main.subject_id)
    ns['flow'][2]['stage']='Exact positive FT4, age72h to<28d; TSH not required'
    shared=main[['subject_id','charttime','renal_time']].merge(expanded[['subject_id','charttime','renal_time']],on='subject_id',suffixes=('_main','_expanded'))
    changed=((shared.charttime_main!=shared.charttime_expanded)|(shared.renal_time_main!=shared.renal_time_expanded)).sum()
    out={'source_sha256':hashlib.sha256((B/'primary_analysis/reaudit_source_local.json').read_bytes()).hexdigest(),'source':'Audited extraction dated 2026-10-15; no new database transaction or write','selection_fraction':760/2656,'thyroid_failures_hierarchical':fail,'renal_failures_hierarchical':renal_fail,'without_TSH_n':len(expanded),'additional_patients':len(extra),'changed_index_pairs_shared_patients':int(changed),'main_n':760,'expanded_flow':ns['flow'],'expanded_exclusions':ns['exclusions']}
    if extra or changed:
     dx=sc['DX'];preterm=set(dx.loc[dx.code.str.match(r'^P07\.[23]'),'hadm_id']);expanded['preterm_recorded']=expanded.hadm_id.isin(preterm).astype(int);expanded['male']=(expanded.gender=='M').astype(int)
     rx=pd.DataFrame(sc['Q']['drugs']['rows']);rx.startdate=pd.to_datetime(rx.startdate);rx.enddate=pd.to_datetime(rx.enddate);rx['text']=rx[['drug_name','drug_name_cn','drug_name_generic']].fillna('').agg(' '.join,axis=1)
     flags=[]
     for r in expanded.itertuples():
      anchor=min(r.charttime,r.renal_time);g=rx[(rx.hadm_id==r.hadm_id)&(rx.startdate<=anchor)];g=g[(g.startdate>=anchor-pd.Timedelta(hours=48))|(g.enddate>=anchor-pd.Timedelta(hours=48))]
      flags.append(int(g.text.str.contains(r'dopamin|多巴胺|epinephrine|肾上腺素|dobutamin|多巴酚丁胺|milrinone|米力农',case=False).any()))
     expanded['vaso_bounded']=flags
     for c in ['FT4','Cr','CysC']:expanded['log_'+c]=np.log(expanded[c])
     frozen=json.loads((B/'primary_analysis/primary_results.json').read_text(encoding='utf-8'));r=next(x for x in frozen['models'] if x['label']=='final_primary')
     m=smf.ols(r['formula'],expanded).fit(cov_type='HC3',use_t=True);lo,hi=m.conf_int().loc['log_FT4'];delta=np.log(18.485/12.0975)
     out['sensitivity']={'formula':r['formula'],'n':int(m.nobs),'beta':float(m.params.log_FT4),'p':float(m.pvalues.log_FT4),'percent':float(100*np.expm1(m.params.log_FT4*delta)),'low':float(100*np.expm1(lo*delta)),'high':float(100*np.expm1(hi*delta)),'rank':int(np.linalg.matrix_rank(m.model.exog)),'max_leverage':float(m.get_influence().hat_matrix_diag.max())}
    else:out['sensitivity']='Not refit: removing TSH yielded identical patients and paired timestamps.'
    (O/'evidence/TSH_selection_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2,default=int),encoding='utf-8')


def export_coefficients(scope):
    from pathlib import Path
    import json,hashlib,io,contextlib,math
    import numpy as np,pandas as pd
    import statsmodels.formula.api as smf
    from scipy.stats import t
    B=WORK;O=WORK/'verification'
    original=B/'primary_analysis/primary_results.json'
    before=hashlib.sha256(original.read_bytes()).hexdigest()
    frozen=json.loads(original.read_text(encoding='utf-8'));rows=[];checks=[]
    def add(label,m):
     ci=m.conf_int()
     for k in m.params.index:
      rows.append(dict(model=label,term=k,n_observations=int(m.nobs),df_resid=float(m.df_resid),beta=float(m.params[k]),standard_error=float(m.bse[k]),ci_low=float(ci.loc[k,0]),ci_high=float(ci.loc[k,1]),p=float(m.pvalues[k]),formula=m.model.formula,inference='HC3/t' if label!='Huber_RLM' else 'H3/normal'))
    for row in frozen['models']:
     m=scope['fits'][row['label']];term='log_'+row['hormone'];r=next(x for x in scope['models'] if x['label']==row['label'])
     for k in ['n','beta','beta_low','beta_high','p','percent','low','high','r2','max_leverage']:np.testing.assert_allclose(r[k],row[k],rtol=1e-9,atol=1e-11)
     add(row['label'],m);checks.append(row['label'])
    add('Huber_RLM',scope['rlm'])
    nl=json.loads((B/'supplementary_analysis/evidence/FT4_functional_form.json').read_text(encoding='utf-8'))
    m=smf.ols(nl['formula'],scope['d']).fit(cov_type='HC3',use_t=True);add('exploratory_FT4_spline',m)
    pd.DataFrame(rows).to_csv(O/'evidence/all_model_parameters_full_precision.csv',index=False)
    ts=json.loads((B/'supplementary_analysis/evidence/TSH_selection_audit.json').read_text(encoding='utf-8'))['sensitivity']
    contrast=scope['iqrs']['FT4'];lo=math.log1p(ts['low']/100)/contrast;hi=math.log1p(ts['high']/100)/contrast;df=ts['n']-ts['rank'];se=(hi-lo)/(2*t.ppf(.975,df))
    np.testing.assert_allclose(2*t.sf(abs(ts['beta']/se),df),ts['p'],rtol=1e-10)
    pd.DataFrame([dict(model='no_TSH_requirement',term='log_FT4',n=ts['n'],df_resid=df,beta=ts['beta'],standard_error=se,ci_low=lo,ci_high=hi,p=ts['p'],formula=ts['formula'],source='SE recovered algebraically from the stored full-precision HC3/t CI and rank; P independently agrees')]).to_csv(O/'evidence/no_TSH_exposure_precision.csv',index=False)
    assert before==hashlib.sha256(original.read_bytes()).hexdigest()
    (O/'evidence/parameter_export_receipt.json').write_text(json.dumps(dict(frozen_models_matched=checks,coefficient_rows=len(rows),frozen_results_sha256=before,frozen_results_unchanged=True,new_analysis_definitions=0,bootstrap='Existing exact 2000-resample reproduction receipt retained; not repeated in this display-precision audit.'),indent=2),encoding='utf-8')


def selection_comparison():
    from pathlib import Path
    import json,hashlib
    import numpy as np,pandas as pd
    B=WORK;O=WORK/'selection_analysis';E=O/'evidence'
    E.mkdir(parents=True,exist_ok=True)
    source=B/'primary_analysis/reaudit_source_local.json';x=json.loads(source.read_text(encoding='utf-8'));q=x['queries']
    a=pd.DataFrame(q['cohort']['rows']);dx=pd.DataFrame(q['diagnoses']['rows']);final=pd.read_json(B/'primary_analysis/final_cohort_verified_local.json');pairs=pd.read_json(B/'primary_analysis/paired_candidates_local.json')
    for col in ['admittime','dischtime','dob']:a[col]=pd.to_datetime(a[col])
    a=a.sort_values(['subject_id','admittime','hadm_id']).drop_duplicates('subject_id').copy()
    assert len(a)==a.subject_id.nunique()==2656
    a['included']=a.subject_id.isin(final.subject_id);a['paired']=a.subject_id.isin(pairs.subject_id)
    assert (a.included.sum(),(~a.included).sum(),(~a.paired).sum(),(a.paired&~a.included).sum())==(760,1896,1889,7)
    dx=dx.merge(a[['subject_id','hadm_id']],on=['subject_id','hadm_id'],validate='many_to_one')
    dx['code']=dx.icd10_code.fillna('').str.strip()
    rules={'prematurity':r'^P07\.[23](\.|$)','sepsis':r'^(P36|A40|A41)(\.|$)','RDS':r'^P22\.0(\.|$)','low_birthweight':r'^P07\.[01](\.|$)'}
    for k,pat in rules.items():a[k]=a.hadm_id.isin(dx.loc[dx.code.str.match(pat),'hadm_id']).astype(int)
    a['admission_age_days']=(a.admittime-a.dob).dt.total_seconds()/86400;a['male']=a.gender.map({'M':1,'F':0})
    a['any_diagnosis']=a.hadm_id.isin(dx.hadm_id).astype(int)
    principal=dx[(dx.diag_category=='Admitting Diagnosis')&(dx.seq_num==1)].copy()
    principal=principal.drop_duplicates(['hadm_id','icd10_code_cn'])
    multi=principal.groupby('hadm_id').size();ambiguous=set(multi[multi>1].index)
    a['admitting_ambiguous']=a.hadm_id.isin(ambiguous).astype(int)
    principal=principal[~principal.hadm_id.isin(ambiguous)]
    a['admitting_primary_record']=a.hadm_id.isin(principal.hadm_id).astype(int)
    for k,pat in rules.items():a['admitting_'+k]=a.hadm_id.isin(principal.loc[principal.code.str.match(pat),'hadm_id']).astype(int)
    a['admitting_other']=a.admitting_primary_record-a[[f'admitting_{k}' for k in rules]].sum(axis=1)
    assert a.admitting_other.ge(0).all()
    labels={'admission_age_days':'Age at admission, days','male':'Male sex','any_diagnosis':'Any recorded diagnosis in selected hospitalization','prematurity':'Recorded prematurity (P07.2/P07.3)','sepsis':'Recorded sepsis (P36/A40/A41)','RDS':'Recorded neonatal RDS (P22.0)','low_birthweight':'Recorded low birthweight (P07.0/P07.1)','admitting_primary_record':'Unambiguous principal admitting diagnosis available','admitting_ambiguous':'Multiple principal admitting diagnosis codes','admitting_prematurity':'Principal admitting diagnosis: prematurity','admitting_sepsis':'Principal admitting diagnosis: sepsis','admitting_RDS':'Principal admitting diagnosis: RDS','admitting_low_birthweight':'Principal admitting diagnosis: low birthweight','admitting_other':'Principal admitting diagnosis: other'}
    labels={k:v.replace('Unambiguous principal admitting diagnosis available','Single admitting-diagnosis code with seq_num=1').replace('Multiple principal admitting diagnosis codes','Multiple admitting-diagnosis codes with seq_num=1').replace('Principal admitting diagnosis:','Single seq_num=1 admitting code:') for k,v in labels.items()}
    def summ(v,continuous):
     v=v.dropna();out={'available':len(v)}
     if continuous:out.update(median=float(v.median()),q25=float(v.quantile(.25)),q75=float(v.quantile(.75)),mean=float(v.mean()),sd=float(v.std(ddof=1)))
     else:out.update(count=int(v.sum()),proportion=float(v.mean()))
     return out
    rows=[]
    for cmp,mask in [('all_not_included',~a.included),('without_qualifying_pair',~a.paired)]:
     for k,label in labels.items():
      v=a.loc[a.included,k];w=a.loc[mask,k];cont=k=='admission_age_days';s1=summ(v,cont);s0=summ(w,cont)
      sd=np.sqrt((v.var(ddof=1)+w.var(ddof=1))/2) if cont else np.sqrt((v.mean()*(1-v.mean())+w.mean()*(1-w.mean()))/2)
      smd=float((v.mean()-w.mean())/sd) if sd>0 else None
      row={'comparison':cmp,'variable':k,'label':label,'included_n':len(v),'comparator_n':len(w),'included_missing':int(v.isna().sum()),'comparator_missing':int(w.isna().sum()),'smd':smd}
      row.update({'included_'+kk:vv for kk,vv in s1.items()});row.update({'comparator_'+kk:vv for kk,vv in s0.items()})
      fmt=lambda ss:f'{ss["median"]:.2f} [{ss["q25"]:.2f}, {ss["q75"]:.2f}]' if cont else f'{ss["count"]} ({ss["proportion"]*100:.1f}%)'
      row.update(included_display=fmt(s1),comparator_display=fmt(s0));rows.append(row)
    pd.DataFrame(rows).to_csv(E/'selection_comparison_full_precision.csv',index=False)
    # Aggregate code audit, no identifiers or individual rows exported.
    codes=dx[dx.code.str.match('|'.join('('+p+')' for p in rules.values()))].groupby(['code','icd10_code_cn','title_cn','title'],dropna=False).agg(diagnosis_records=('row_id','size'),distinct_hospitalizations=('hadm_id','nunique')).reset_index()
    codes.to_csv(E/'selection_diagnosis_code_audit.csv',index=False,encoding='utf-8-sig')
    result={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'analysis_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'denominator':2656,'included':760,'all_not_included':1896,'without_qualifying_pair':1889,'paired_clinical_exclusions':7,'unmatched_diagnosis_dictionary_records':int(dx.icd10_code.isna().sum()),'rules':rules,'principal_admitting_definition':"diag_category == 'Admitting Diagnosis' and seq_num == 1; duplicate identical local codes collapsed; hospitalizations with multiple distinct principal codes labelled ambiguous and not assigned a disease category",'principal_ambiguous_hospitalizations':len(ambiguous),'results':rows,'database_writes':0,'new_hypothesis_tests':0,'limitations':['Recorded diagnoses lack onset timestamps; whole-admission prevalence is not baseline disease prevalence.','Principal-admitting rows use entire group denominator and must be read with availability; missing/ambiguous is not other diagnosis.','Actual gestational age and verified birthweight unavailable; codes are imperfect proxies.','Observed differences describe selection but do not identify its mechanism or correct selection bias.','This analysis was added after review of the existing results.']}
    (E/'selection_comparison.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description='Reproduce the PIC neonatal FT4-cystatin C analyses.')
    parser.add_argument('--extract', action='store_true', help='Read the authorized PIC database before analysis.')
    args = parser.parse_args()
    initialize()
    if args.extract:
        extract_sources()
    build_pairs(6)
    build_pairs(24)
    print('1/5 Paired cohorts prepared.', flush=True)
    initial_analysis()
    if args.extract:
        extract_covariates()
    prepare_covariates()
    extended_analysis()
    admission_weight_analysis()
    print('2/5 Clinical covariates prepared.', flush=True)
    audit_cohort()
    fitted = primary_analysis()
    print('3/5 Primary and sensitivity analyses completed.', flush=True)
    supplementary_analysis()
    thyroid_selection_analysis()
    export_coefficients(fitted)
    selection_comparison()
    print('4/5 Supplementary analyses completed.', flush=True)
    subprocess.run([sys.executable, str(HERE/'cohort_flowchart.py'),
                    '--config', str(HERE/'cohort_config.json'),
                    '--output-dir', str(WORK/'figures/cohort')], check=True)
    print('5/5 Cohort flowchart generated. Run figures.R for the remaining figures.', flush=True)

if __name__ == '__main__':
    main()
