import os
DATA_DIR=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
import numpy as np, pandas as pd, time, os
from cloudvolume import CloudVolume
OUT=DATA_DIR+'/neck_slabs'
vol=CloudVolume('gs://flyem-male-cns/v1.0/segmentation', mip=0, use_https=True, progress=False, fill_missing=True, cache=DATA_DIR+'/cv_cache')
x0,x1=44800,54800; y0,y1=48600,55400; z0=53952   # mip0 coords, chunk aligned (multiples of 64)
t=time.time(); rows=[]
for zs in range(z0,z0+64,4):
    cut=np.asarray(vol[x0:x1, y0:y1, zs:zs+4])[...,0]
    for k in range(cut.shape[2]):
        ids,cnt=np.unique(cut[:,:,k],return_counts=True); rows.append(pd.DataFrame({'body':ids,'area_vox':cnt,'z':zs+k}))
    if zs==z0+32: np.savez_compressed(f'{OUT}/plane_mip0_z{zs}.npz',plane=cut[:,:,0].astype(np.uint64),x0=x0,y0=y0,z=zs)
    print(zs,round(time.time()-t),'s',flush=True)
df=pd.concat(rows); df['slab']='B0'; df.to_parquet(f'{OUT}/areas_mip0_B.parquet'); print('DONE',len(df),round(time.time()-t),'s')
