import os
DATA_DIR=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
import pandas as pd, numpy as np
D=DATA_DIR; R='results'
a0=pd.read_parquet(f'{D}/neck_slabs/areas_mip0_B.parquet'); a0=a0[a0.body!=0]
g0=a0.groupby('body').agg(n=('z','nunique'),area=('area_vox','median')); g0=g0[g0.n>=32]; g0['area_um2_mip0']=g0.area*0.008**2
s=pd.read_parquet(f'{R}/caliber_per_slab.parquet'); s1=s[(s.slab=='B')&(s.n_planes>=32)].set_index('body')
j=g0.join(s1[['area_um2_raw','diam_um']],how='inner'); j['ratio']=j.area_um2_mip0/j.area_um2_raw
ch=pd.read_parquet(f'{R}/channel_axons_male.parquet'); j=j[j.index.isin(ch.body)]
j['d0']=2*np.sqrt(j.area_um2_mip0/np.pi); j['bin']=pd.cut(j.d0,[0,0.2,0.3,0.5,1,3,20])
print('n',len(j)); print(j.groupby('bin').agg(n=('ratio','size'),ratio_med=('ratio','median'),ratio_q25=('ratio',lambda v:v.quantile(.25)),ratio_q75=('ratio',lambda v:v.quantile(.75))).round(3).to_string())
from scipy import stats; print('Spearman',stats.spearmanr(j.area_um2_mip0,j.area_um2_raw)); print('overall median ratio',round(j.ratio.median(),3))
j.drop(columns=['bin']).to_parquet(f'{R}/mip_validation.parquet')
