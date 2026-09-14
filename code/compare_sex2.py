import pandas as pd, numpy as np, json
from scipy import stats
R='results'
m=pd.read_parquet(f'{R}/channel_axons_male.parquet'); f=pd.read_parquet(f'{R}/channel_axons_female.parquet')
m['dimcat']=m.dimorphism.fillna('isomorphic').replace({'potentially male-specific':'male-specific','potentially sexually dimorphic':'sexually dimorphic'})
fd=f.dim.astype(str).str.lower(); f['dimcat']=np.where(fd.str.contains('female'),'female-specific',np.where(fd.str.contains('dimorphic'),'sexually dimorphic','isomorphic'))
print('female dim raw:',f.dim.value_counts(dropna=False).to_dict())
# artifacts: male axons < 0.1 um
art=m[m.diam_um<0.1]; print('male channel axons < 0.1 um:',len(art), art.superclass.value_counts().to_dict())
out={}
out['female_measured']=int(len(f)); out['female_counts']=f.groupby('superclass').size().to_dict()
for lab,df in [('male',m),('female',f)]:
    dn=df[df.dir=='down']; up=df[df.dir=='up']
    out[f'{lab}_share_count']=round(len(dn)/len(df),3); out[f'{lab}_share_sumd']=round(dn.diam_um.sum()/df.diam_um.sum(),3); out[f'{lab}_share_area']=round(dn.area_um2.sum()/df.area_um2.sum(),3)
    betas=np.linspace(0.25,2.5,226); ratio=np.array([(dn.diam_um**b).sum()/(up.diam_um**b).sum() for b in betas]); out[f'{lab}_parity_beta']=round(float(betas[np.argmin(np.abs(ratio-1))]),2)
    out[f'{lab}_ratio_beta1']=round(float(ratio[np.argmin(np.abs(betas-1))]),3)
    iso=df[df.dimcat=='isomorphic']; out[f'{lab}_iso_median']={sc:round(iso[iso.superclass==sc].diam_um.median(),3) for sc in ['descending_neuron','ascending_neuron','sensory_ascending']}
    out[f'{lab}_gini']=round(1-2*np.trapz(np.sort(df.area_um2.values).cumsum()/df.area_um2.sum(),dx=1/len(df)),3)
    c=df.sort_values('area_um2',ascending=False); cs=c.area_um2.cumsum()/c.area_um2.sum(); out[f'{lab}_n_half_area']=int((cs<0.5).sum()+1); out[f'{lab}_frac_half_area']=round(((cs<0.5).sum()+1)/len(df),3)
# relative caliber (normalized to isomorphic DN/AN median of same sex & class)
def rel(df):
    med=df[df.dimcat=='isomorphic'].groupby('superclass').diam_um.median(); return df.diam_um/df.superclass.map(med)
m['rel']=rel(m); f['rel']=rel(f)
for lab,df,spec in [('male',m,'male-specific'),('female',f,'female-specific')]:
    for sc in ['descending_neuron','ascending_neuron']:
        a=df[(df.superclass==sc)&(df.dimcat==spec)]; b=df[(df.superclass==sc)&(df.dimcat=='isomorphic')]
        if len(a)>=3:
            mw=stats.mannwhitneyu(a.rel,b.rel,alternative='two-sided'); out[f'{lab}_{spec}_{sc}']={'n':int(len(a)),'rel_median':round(a.rel.median(),2),'abs_median':round(a.diam_um.median(),3),'iso_median':round(b.diam_um.median(),3),'MWU_p':float(f'{mw.pvalue:.2g}')}
        d=df[(df.superclass==sc)&(df.dimcat=='sexually dimorphic')]
        if len(d)>=3:
            mw=stats.mannwhitneyu(d.rel,b.rel,alternative='two-sided'); out[f'{lab}_dimorphic_{sc}']={'n':int(len(d)),'rel_median':round(d.rel.median(),2),'MWU_p':float(f'{mw.pvalue:.2g}')}
# type-matched relative calibers (exclude male artifacts)
mt=m[(m.diam_um>=0.1)&m.type.notna()].groupby('type').agg(male=('rel','median'),male_abs=('diam_um','median'),sc=('superclass','first'),dim=('dimcat','first'))
ft=f[f.mtype.notna()].groupby('mtype').agg(female=('rel','median'),female_abs=('diam_um','median'))
j=mt.join(ft,how='inner'); out['n_matched_types']=int(len(j))
out['matched_abs_spearman']=round(stats.spearmanr(j.male_abs,j.female_abs)[0],3); out['matched_abs_median_log2_f_over_m']=round(float(np.median(np.log2(j.female_abs/j.male_abs))),3)
boot=[np.median(np.log2(j.female_abs.sample(frac=1,replace=True,random_state=i)/j.male_abs.sample(frac=1,replace=True,random_state=i))) for i in range(200)]
out['matched_abs_median_log2_ci']=[round(np.percentile(boot,2.5),3),round(np.percentile(boot,97.5),3)]
out['matched_rel_spearman']=round(stats.spearmanr(j.male,j.female)[0],3); out['matched_rel_median_abs_log2']=round(float(np.median(np.abs(np.log2(j.female/j.male)))),3)
# within-sex L-R for same types (male) for comparison
lr=m[(m.diam_um>=0.1)&m.type.notna()&m.somaSide.isin(['L','R'])].groupby(['type','somaSide']).diam_um.median().unstack().dropna(); out['male_LR_median_abs_log2']=round(float(np.median(np.abs(np.log2(lr.L/lr.R)))),3)
j['log2_rel_m_over_f']=np.log2(j.male/j.female); j=j.sort_values('log2_rel_m_over_f',ascending=False)
print('top male-enlarged (relative):\n',j.head(10).round(2).to_string()); print('top female-enlarged (relative):\n',j.tail(8).round(2).to_string())
print('dimorphic types relative:\n',j[j.dim=='sexually dimorphic'].round(2).head(20).to_string())
j.to_parquet(f'{R}/type_matched_relative.parquet'); f.to_parquet(f'{R}/channel_axons_female.parquet'); m.to_parquet(f'{R}/channel_axons_male.parquet')
print(json.dumps(out,indent=1)); json.dump(out,open(f'{R}/compare_sex2.json','w'),indent=1)
# Table S1 female sheet
f[['root_888','superclass','cell_type','malecns_cell_type','side','dim','n_planes','area_um2','diam_um']].rename(columns={'area_um2':'neck_area_um2','diam_um':'neck_diameter_um','dim':'dimorphism'}).to_csv('submission/TableS1_female_neck_axons.csv',index=False)
