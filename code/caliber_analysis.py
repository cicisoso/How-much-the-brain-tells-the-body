"""Merge per-plane cross-sectional areas with annotations, correct for axon tilt, and produce per-axon caliber table."""
import os
DATA_DIR=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
import numpy as np, pandas as pd, glob, os
D=DATA_DIR; OUT='../results'; os.makedirs(OUT,exist_ok=True)
VOX=16e-3  # um per mip-1 voxel
ann=pd.read_feather(f'{D}/body-annotations-male-cns-v1.0-minconf-0.5.feather')
neck=pd.read_csv(f'{D}/neck_neurons.csv')
areas=pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f'{D}/neck_slabs/areas_*.parquet'))])
areas=areas[areas.body!=0]
# per body, per slab: median plane area; require presence in >= 32 of 64 planes
g=areas.groupby(['slab','body']).agg(n_planes=('z','nunique'),area_med=('area_vox','median'),area_min=('area_vox','min'),area_max=('area_vox','max')).reset_index()
g['area_um2_raw']=g.area_med*VOX*VOX
# tilt correction from skeleton direction near each slab centre (mip0 z = 2*mip1 z)
def tilt(b, zc_mip0):
    p=f'{D}/swc/{b}.swc'
    if not os.path.exists(p): return np.nan
    d=pd.read_csv(p,sep=r'\s+',comment='#',header=None,names=['id','t','x','y','z','r','p'])
    s=d[(d.z>zc_mip0-400)&(d.z<zc_mip0+400)&(d.x>45000)&(d.x<54600)&(d.y>48800)&(d.y<55300)]
    if len(s)<3: return np.nan
    X=s[['x','y','z']].values.astype(float); X-=X.mean(0)
    u,sv,vt=np.linalg.svd(X,full_matrices=False); v=vt[0]
    return abs(v[2])/np.linalg.norm(v)  # cos(theta) between axon direction and z
slab_z={'A':26240+32,'B':26944+32,'C':27648+32}
g['cos_theta']=[tilt(b,2*slab_z[s]) if b in set(neck.bodyId) else np.nan for s,b in zip(g.slab,g.body)]
g['area_um2']=g.area_um2_raw*g.cos_theta.fillna(1.0)
g['diam_um']=2*np.sqrt(g.area_um2/np.pi)
g.to_parquet(f'{OUT}/caliber_per_slab.parquet')
# per-body summary across slabs
s=g[g.n_planes>=32].groupby('body').agg(n_slabs=('slab','nunique'),area_um2=('area_um2','median'),diam_um=('diam_um','median'),diam_cv=('diam_um',lambda v: v.std()/v.mean() if len(v)>1 else np.nan)).reset_index()
s=s.merge(ann[['bodyId','superclass','class','type','instance','somaSide','rootSide','status','statusLabel','dimorphism','birthtime','fruDsx','flywireType','mancType']],left_on='body',right_on='bodyId',how='left')
s['neck_listed']=s.body.isin(neck.bodyId)
s.to_parquet(f'{OUT}/caliber_per_body.parquet')
print('bodies with >=32 planes in >=1 slab:',len(s)); print(s.groupby(s.neck_listed).size())
print(s[s.neck_listed].groupby('superclass').diam_um.describe().to_string())
print('unlisted bodies status:'); print(s[~s.neck_listed].status.value_counts(dropna=False).head(10).to_string())
print('GF:', s[s.type=='DNp01'][['instance','diam_um','area_um2','n_slabs','diam_cv']].to_string())
