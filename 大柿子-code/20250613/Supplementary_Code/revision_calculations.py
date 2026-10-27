"""Reproduce the revision-only calculations using Python 3 (standard library)."""
from pathlib import Path
import csv,math,json
BASE=Path(__file__).resolve().parent
OUT=BASE/'calculated';OUT.mkdir(exist_ok=True)
def load(name,delimiter=','):
    with (BASE/'inputs'/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter=delimiter))
def save(name,rows):
    with (OUT/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
priors=[1e-6,5e-6,1e-5,5e-5];sens=[]
for r in load('finnGen_coloc_input.csv'):
    h=float(r['pph4'])
    for prior in priors:
        ratio=prior/1e-5;r['PPH4_p12_'+str(prior)]=ratio*h/(1-h+ratio*h)
    r['p12_at_PPH4_070']=1e-5*.7*(1-h)/(.3*h) if h>0 else ''
    sens.append(r)
save('coloc_prior_sensitivity.csv',sens)
vr=[]
for endpoint in ['L40','20002_1477']:
    for r in load('VARS2_mQTL_'+endpoint+'.smr','\t'):
        b,se=float(r['b_SMR']),float(r['se_SMR'])
        r['endpoint']=endpoint;r['beta_lower_95']=b-1.96*se;r['beta_upper_95']=b+1.96*se
        r['F_QTL']=(float(r['b_eQTL'])/float(r['se_eQTL']))**2
        r['absolute_ratio_error']=abs(b-float(r['b_GWAS'])/float(r['b_eQTL']))
        vr.append(r)
save('VARS2_coefficients_and_checks.csv',vr)
primary=[r for r in sens if r['baseline_primary']=='True']
summary={'baseline_candidate_annotation_rows':len(primary),'p12_values':priors,'PPH4_above_070_counts':[sum(r['PPH4_p12_'+str(p)]>.7 for r in primary) for p in priors],'VARS2_rows':len(vr),'maximum_effect_ratio_error':max(r['absolute_ratio_error'] for r in vr)}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
