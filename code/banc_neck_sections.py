import pandas as pd, numpy as np, time, trimesh, os
from cloudvolume import CloudVolume
from multiprocessing import Pool
m=pd.read_feather('banc_888_meta.feather'); nk=pd.read_parquet('neck_connective_y92500.parquet')
m=m[m.root_888.notna()].copy(); m['r888']=m.root_888.astype('int64')
keep=m[m.super_class.isin(['descending','ascending','sensory_ascending','sensory_descending'])]
nk2=nk[nk.pt_root_id.isin(keep.r888)].copy()
nk2=nk2.merge(keep[['r888','super_class','cell_type','malecns_cell_type','side','sexually_dimorphic','cell_class','volume_nm3','l2_cable_length_um']],left_on='pt_root_id',right_on='r888',how='left')
print('neurons to section',len(nk2),nk2.super_class.value_counts().to_dict(),flush=True)
SC=np.array([4,4,45.0]); PLANES=[92000,92500,93000]  # y in 4-nm voxels; +-2 um
def work(args):
    rid,pt=args
    try:
        cv=CloudVolume('gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_meshes',use_https=True,progress=False)
        mesh=cv.mesh.get(int(rid),fuse=False)[int(rid)]
        tm=trimesh.Trimesh(vertices=mesh.vertices,faces=mesh.faces,process=False)
        seed=np.array(pt,dtype=float)*SC; out=[]
        for yv in PLANES:
            y=yv*4.0
            sec=tm.section(plane_origin=[0,y,0],plane_normal=[0,1,0])
            if sec is None: out.append((rid,yv,np.nan,np.nan,0)); continue
            p2,T=sec.to_2D(); polys=p2.polygons_full
            if len(polys)==0: out.append((rid,yv,np.nan,np.nan,0)); continue
            # choose polygon whose centroid (back in 3D) is closest to seed (in x,z)
            best=None
            for pg in polys:
                c2=np.array(pg.centroid.coords[0]); c3=trimesh.transform_points(np.array([[c2[0],c2[1],0.0]]),T)[0]
                dist=np.hypot(c3[0]-seed[0],c3[2]-seed[2])
                if best is None or dist<best[0]: best=(dist,pg.area)
            out.append((rid,yv,best[1],best[0],len(polys)))
        return out
    except Exception as e:
        return [(rid,-1,np.nan,np.nan,-1)]
if __name__=='__main__':
    t=time.time(); args=list(zip(nk2.pt_root_id.astype('int64'),nk2.pt_position))
    with Pool(10) as p: res=p.map(work,args,chunksize=4)
    rows=[r for rr in res for r in rr]
    df=pd.DataFrame(rows,columns=['root_888','plane_y','area_nm2','seed_dist_nm','n_polys'])
    df['area_um2']=df.area_nm2/1e6; df['diam_um']=2*np.sqrt(df.area_um2/np.pi)
    df=df.merge(nk2[['pt_root_id','super_class','cell_type','malecns_cell_type','side','sexually_dimorphic','cell_class','volume_nm3','l2_cable_length_um','tag']],left_on='root_888',right_on='pt_root_id',how='left')
    df.to_parquet('banc_neck_sections.parquet'); print('done',len(df),'rows',round(time.time()-t),'s',flush=True)
    print(df.groupby('super_class').diam_um.describe().to_string())
