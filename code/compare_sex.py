"""Male (MaleCNS, EM segmentation) vs female (BANC, mesh sections) neck calibers."""
import os
DATA_DIR=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
import pandas as pd, numpy as np, json
from scipy import stats
R='results'; B=DATA_DIR+'/banc'
male=pd.read_parquet(f'{R}/channel_axons_male.parquet')
b=pd.read_parquet(f'{B}/banc_neck_sections.parquet')
b=b[b.plane_y>0]
fb=b.groupby('root_888').agg(n_planes=('plane_y',lambda v: v.notna().sum()),area_um2=('area_um2','median'),diam_um=('diam_um','median'),seed_dist=('seed_dist_nm','median'),super_class=('super_class','first'),cell_type=('cell_type','first'),malecns_cell_type=('malecns_cell_type','first'),side=('side','first'),dim=('sexually_dimorphic','first')).reset_index()
fb=fb[fb.area_um2.notna()]
print('female measured',len(fb),fb.super_class.value_counts().to_dict())
print('seed distance to chosen polygon (nm) median',fb.seed_dist.median(),' 90%',fb.seed_dist.quantile(.9))
fb=fb[fb.seed_dist<3000]  # polygon must be within 3 um of the annotated neck point
print('after seed filter',len(fb))
sc_map={'descending':'descending_neuron','ascending':'ascending_neuron','sensory_ascending':'sensory_ascending','sensory_descending':'sensory_descending'}
fb['superclass']=fb.super_class.map(sc_map); fb['dir']=np.where(fb.superclass.isin(['descending_neuron','sensory_descending']),'down','up')
out={}
out['female_counts']=fb.groupby('superclass').size().to_dict(); out['female_area']=fb.groupby('superclass').area_um2.sum().round(1).to_dict()
out['female_median_diam']=fb.groupby('superclass').diam_um.median().round(3).to_dict(); out['female_sum_diam']=fb.groupby('superclass').diam_um.sum().round(1).to_dict()
for beta in [1.0,2.0]:
    dn=(fb[fb.dir=='down'].diam_um**beta).sum(); up=(fb[fb.dir=='up'].diam_um**beta).sum(); out[f'female_ratio_down_up_beta{beta}']=round(dn/up,3)
# type-matched comparison: female type -> male type via malecns_cell_type or cell_type
mt=male[male.type.notna()].groupby('type').diam_um.median()
fb['mtype']=fb.malecns_cell_type.where(fb.malecns_cell_type.notna(),fb.cell_type)
ft=fb[fb.mtype.notna()].groupby('mtype').diam_um.median()
j=pd.concat([mt.rename('male'),ft.rename('female')],axis=1).dropna()
out['n_matched_types']=len(j); out['spearman_male_female_types']=[round(x,3) for x in stats.spearmanr(j.male,j.female)]
out['pearson_log']=round(stats.pearsonr(np.log(j.male),np.log(j.female))[0],3)
out['median_log2_ratio_f_over_m']=round(float(np.median(np.log2(j.female/j.male))),3); out['median_abs_log2_ratio']=round(float(np.median(np.abs(np.log2(j.female/j.male)))),3)
# L-R within female for comparison
lr=fb[fb.side.isin(['left','right'])&fb.mtype.notna()].groupby(['mtype','side']).diam_um.median().unstack().dropna()
out['female_LR_spearman']=[round(x,3) for x in stats.spearmanr(lr.left,lr.right)]; out['female_LR_median_abs_log2']=round(float(np.median(np.abs(np.log2(lr.left/lr.right)))),3)
# sex-specific types
print('female dimorphism categories:',fb.dim.value_counts(dropna=False).to_dict())
fs=fb[fb.dim.astype(str).str.contains('female',case=False)]; print('female-specific/dimorphic neck neurons:\n',fs.groupby('superclass').agg(n=('root_888','size'),d_med=('diam_um','median'),d_mean=('diam_um','mean')).round(3).to_string())
print(fs.sort_values('diam_um',ascending=False).head(12)[['cell_type','superclass','diam_um','dim']].to_string())
# types with largest sex differences (beyond L-R)
j['log2_f_m']=np.log2(j.female/j.male); print('largest female>male:\n',j.sort_values('log2_f_m',ascending=False).head(10).round(3).to_string()); print('largest male>female:\n',j.sort_values('log2_f_m').head(10).round(3).to_string())
print(json.dumps(out,indent=1)); json.dump(out,open(f'{R}/compare_sex.json','w'),indent=1)
fb.to_parquet(f'{R}/channel_axons_female.parquet'); j.to_parquet(f'{R}/type_matched_male_female.parquet')
