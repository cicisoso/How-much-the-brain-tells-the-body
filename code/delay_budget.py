import os
DATA_DIR=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
import numpy as np, pandas as pd, os
from collections import defaultdict, deque
D=DATA_DIR
neck=pd.read_csv(f'{D}/neck_neurons.csv')
ZN=54000.0; VOX=0.008  # um per voxel
def analyze(b):
    d=pd.read_csv(f'{D}/swc/{b}.swc',sep=r'\s+',comment='#',header=None,names=['id','t','x','y','z','r','p'])
    idx={i:k for k,i in enumerate(d.id.values)}
    xyz=d[['x','y','z']].values*VOX; par=d.p.values; ids=d.id.values
    n=len(d); children=defaultdict(list)
    for k in range(n):
        if par[k]!=-1 and par[k] in idx: children[idx[par[k]]].append(k)
    # node nearest neck plane along z with x,y inside neck bbox
    inbox=(d.x.values>45000)&(d.x.values<54600)&(d.y.values>48800)&(d.y.values<55300)
    cand=np.where(inbox)[0]
    if len(cand)==0: return None
    k0=cand[np.argmin(np.abs(xyz[cand,2]-ZN*VOX))]
    if abs(xyz[k0,2]-ZN*VOX)>3: return None
    # undirected adjacency
    adj=defaultdict(list)
    for k in range(n):
        if par[k]!=-1 and par[k] in idx:
            j=idx[par[k]]; adj[k].append(j); adj[j].append(k)
    # walk from k0 in the two directions along the unbranched path until a 'real' branch point (side branch total cable >= 3 um)
    def subtree_len(start,forbid):
        seen={forbid,start}; q=deque([start]); L=0.0
        while q:
            u=q.popleft()
            for v in adj[u]:
                if v not in seen:
                    seen.add(v); L+=np.linalg.norm(xyz[u]-xyz[v]); q.append(v)
        return L
    def walk(start,prev):
        L=0.0; u=start; p=prev
        while True:
            nb=[v for v in adj[u] if v!=p]
            if len(nb)==0: return L,u,'end'
            if len(nb)==1:
                L+=np.linalg.norm(xyz[u]-xyz[nb[0]]); p,u=u,nb[0]; continue
            # branch: keep going along the longest branch if the others are tiny
            lens=[(subtree_len(v,u),v) for v in nb]; lens.sort(reverse=True)
            if lens[1][0]<3.0:
                v=lens[0][1]; L+=np.linalg.norm(xyz[u]-xyz[v]); p,u=u,v; continue
            return L,u,'branch'
    nbrs=adj[k0]
    if len(nbrs)<2: return None
    # direction: neighbor with smaller z is toward brain (z decreases toward brain? brain at low z)
    nb_sorted=sorted(nbrs,key=lambda v: xyz[v,2])
    Lb,kb,tb=walk(nb_sorted[0],k0); Lv,kv,tv=walk(nb_sorted[-1],k0)
    rmed=d.r.values[k0]
    return dict(bodyId=b,len_to_brain_um=Lb,len_to_vnc_um=Lv,end_brain=tb,end_vnc=tv,z_brain_end=xyz[kb,2],z_vnc_end=xyz[kv,2],r_skel_neck=rmed,n_nodes=n)
out=[]
for i,b in enumerate(neck.bodyId):
    try:
        r=analyze(b)
        if r: out.append(r)
    except Exception as e:
        print('err',b,e)
    if i%500==0: print(i,flush=True)
R=pd.DataFrame(out); R.to_parquet(f'{D}/neck_delay_budget.parquet'); print(R.describe().to_string()); print('n',len(R))
