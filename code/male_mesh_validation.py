import pandas as pd, numpy as np, trimesh, time
from cloudvolume import CloudVolume
from multiprocessing import Pool
s=pd.read_parquet('results/caliber_per_body.parquet'); lst=s[s.neck_listed]
rng=np.random.default_rng(0)
# stratified sample across diameters: 40 from each quartile + all > 3 um
q=pd.qcut(lst.diam_um,4,labels=False); sel=pd.concat([lst[q==i].sample(40,random_state=i) for i in range(4)]+[lst[lst.diam_um>3]])
Z=(26976*2)*8.0
def work(b):
    try:
        vol=CloudVolume('gs://flyem-male-cns/v1.0/segmentation', mip=1, use_https=True, progress=False)
        mesh=vol.mesh.get(int(b))[int(b)]; tm=trimesh.Trimesh(vertices=mesh.vertices,faces=mesh.faces,process=False)
        out=[]
        for dz in (-320.0,0.0,320.0):
            sec=tm.section(plane_origin=[0,0,Z+dz],plane_normal=[0,0,1])
            if sec is None: out.append(np.nan); continue
            p2,T=sec.to_2D(); polys=p2.polygons_full
            # polygon inside neck bbox (x 45000-54600, y 48800-55300 mip0 -> nm *8)
            best=np.nan
            for pg in polys:
                c2=np.array(pg.centroid.coords[0]); c3=trimesh.transform_points(np.array([[c2[0],c2[1],0.0]]),T)[0]
                if 45000*8<c3[0]<54600*8 and 48800*8<c3[1]<55300*8: best=pg.area/1e6 if np.isnan(best) else max(best,pg.area/1e6)
            out.append(best)
        return (b,np.nanmedian(out),len(mesh.vertices))
    except Exception as e: return (b,np.nan,-1)
if __name__=='__main__':
    t=time.time()
    with Pool(8) as p: res=p.map(work,sel.body.astype('int64').tolist(),chunksize=2)
    df=pd.DataFrame(res,columns=['body','area_mesh_um2','n_verts']).merge(sel[['body','area_um2','diam_um','superclass']],on='body')
    df.to_parquet('results/male_mesh_validation.parquet'); print('done',len(df),round(time.time()-t),'s')
    d=df.dropna(); print('n',len(d),' median ratio mesh/voxel',round((d.area_mesh_um2/d.area_um2).median(),3),' IQR',(d.area_mesh_um2/d.area_um2).quantile([.25,.75]).round(3).tolist())
    from scipy import stats; print('Spearman',stats.spearmanr(d.area_mesh_um2,d.area_um2))
    d['bin']=pd.cut(d.diam_um,[0,0.3,0.5,1,3,20]); print(d.groupby('bin').apply(lambda x: pd.Series({'n':len(x),'ratio_med':(x.area_mesh_um2/x.area_um2).median()})).to_string())
