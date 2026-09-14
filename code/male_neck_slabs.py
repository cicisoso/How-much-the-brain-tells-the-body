import os
DATA_DIR=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
import numpy as np, pandas as pd, time, os
from cloudvolume import CloudVolume
OUT=DATA_DIR+'/neck_slabs'; os.makedirs(OUT,exist_ok=True)
vol=CloudVolume('gs://flyem-male-cns/v1.0/segmentation', mip=1, use_https=True, progress=False, fill_missing=True, cache=DATA_DIR+'/cv_cache')
# mip1 voxel = 16 nm. bbox from skeleton spread (mip0 x 45000-54600, y 48800-55300) with margin
x0,x1=22400,27400; y0,y1=24300,27700
slabs={'A':26240,'B':26944,'C':27648}   # mip1 z starts (chunk aligned: multiples of 64); mip0 z ≈ 52480, 53888, 55296
for name,z0 in slabs.items():
    t=time.time(); rows=[]
    for zs in range(z0,z0+64,8):
        cut=np.asarray(vol[x0:x1, y0:y1, zs:zs+8])[...,0]   # (x,y,z)
        for k in range(cut.shape[2]):
            plane=cut[:,:,k]
            ids,cnt=np.unique(plane,return_counts=True)
            rows.append(pd.DataFrame({'body':ids,'area_vox':cnt,'z':zs+k}))
        if zs==z0+32:
            np.savez_compressed(f'{OUT}/plane_{name}_z{zs}.npz',plane=cut[:,:,0].astype(np.uint64),x0=x0,y0=y0,z=zs)
    df=pd.concat(rows); df['slab']=name
    df.to_parquet(f'{OUT}/areas_{name}.parquet'); print(name,'done',len(df),'rows',round(time.time()-t),'s',flush=True)
print('ALL DONE')
