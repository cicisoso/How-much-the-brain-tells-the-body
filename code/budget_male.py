import os
DATA_DIR=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
import pandas as pd, numpy as np, json
from scipy import stats
D=DATA_DIR; R='results'
s=pd.read_parquet(f'{R}/caliber_per_body.parquet'); ann=pd.read_feather(f'{D}/body-annotations-male-cns-v1.0-minconf-0.5.feather')
lst=s[s.neck_listed].merge(ann[['bodyId','subclass','somaNeuromere']],on='bodyId',how='left')
lst['dir']=np.where(lst.superclass.isin(['descending_neuron','sensory_descending','efferent_descending']),'down','up')
chan=lst[lst.superclass.isin(['descending_neuron','ascending_neuron','sensory_ascending','sensory_descending'])].copy()
# information models: rate = r1 * d^beta (Hz), bits/spike b; report structural shares too
def budget(df,beta=1.0,r1=10.0,b=2.0):
    return (r1*df.diam_um**beta*b).sum()
out={}
out['counts']=chan.groupby('superclass').size().to_dict()
out['area_um2']=chan.groupby('superclass').area_um2.sum().round(1).to_dict()
out['sum_diam_um']=chan.groupby('superclass').diam_um.sum().round(1).to_dict()
out['median_diam']=chan.groupby('superclass').diam_um.median().round(3).to_dict()
out['mean_diam']=chan.groupby('superclass').diam_um.mean().round(3).to_dict()
for beta in [0.5,1.0,1.5,2.0]:
    dn=budget(chan[chan.dir=='down'],beta); up=budget(chan[chan.dir=='up'],beta)
    out[f'ratio_down_up_beta{beta}']=round(dn/up,3)
    out[f'bits_s_down_beta{beta}']=round(dn); out[f'bits_s_up_beta{beta}']=round(up)
# capacity-bound model (MacKay-McCulloch): H = r log2(e/(r dt)), dt=1ms, r=10*d
dt=1e-3; r=10*chan.diam_um; H=r*np.log2(np.e/(r*dt)); chan['H_bits_s']=H
out['capacity_down']=round(H[chan.dir=='down'].sum()); out['capacity_up']=round(H[chan.dir=='up'].sum())
# per-axon share: top-k axons carrying X% of area
c=chan.sort_values('area_um2',ascending=False); cs=c.area_um2.cumsum()/c.area_um2.sum()
out['n_axons_for_half_area']=int((cs<0.5).sum()+1); out['top10_area_share']=round(c.area_um2.head(10).sum()/c.area_um2.sum(),3)
out['gini_area']=round(1-2*np.trapz(np.sort(chan.area_um2.values).cumsum()/chan.area_um2.sum(),dx=1/len(chan)),3)
# distribution shape: lognormal fit per class
for sc in ['descending_neuron','ascending_neuron','sensory_ascending']:
    d=chan[chan.superclass==sc].diam_um; ld=np.log(d)
    out[f'lognorm_{sc}']={'mu':round(ld.mean(),3),'sigma':round(ld.std(),3),'skew_log':round(stats.skew(ld),3),'d_mode_um':round(np.exp(ld.mean()-ld.std()**2),3)}
# DN by motor domain subclass
sub=chan[chan.superclass=='descending_neuron'].groupby('subclass').agg(n=('body','size'),d_med=('diam_um','median'),area=('area_um2','sum'),sumd=('diam_um','sum')).round(2)
sub['area_share']=(sub.area/sub.area.sum()).round(3); sub['count_share']=(sub.n/sub.n.sum()).round(3)
print('DN by target subclass:\n',sub.sort_values('area',ascending=False).to_string())
# AN by soma neuromere
an=chan[chan.superclass=='ascending_neuron'].groupby('somaNeuromere').agg(n=('body','size'),d_med=('diam_um','median'),area=('area_um2','sum')).round(2); an['area_share']=(an.area/an.area.sum()).round(3)
print('AN by soma neuromere:\n',an.sort_values('area',ascending=False).to_string())
sa=chan[chan.superclass=='sensory_ascending'].groupby('subclass').agg(n=('body','size'),d_med=('diam_um','median'),area=('area_um2','sum')).round(2); sa['area_share']=(sa.area/sa.area.sum()).round(3)
print('SA by organ:\n',sa.sort_values('area',ascending=False).to_string())
# dimorphism
dm=chan.copy(); dm['dim']=dm.dimorphism.fillna('isomorphic').replace({'potentially male-specific':'male-specific','potentially sexually dimorphic':'sexually dimorphic'})
print('dimorphism:\n',dm.groupby(['superclass','dim']).agg(n=('body','size'),d_med=('diam_um','median'),d_mean=('diam_um','mean')).round(3).to_string())
ms=dm[(dm.dim=='male-specific')&(dm.superclass=='descending_neuron')][['type','instance','diam_um','area_um2']]; print(ms.to_string())
mw=stats.mannwhitneyu(dm[(dm.dim!='isomorphic')&(dm.superclass=='descending_neuron')].diam_um, dm[(dm.dim=='isomorphic')&(dm.superclass=='descending_neuron')].diam_um); print('MWU DN dimorphic vs iso',mw)
print('top 15 axons:\n',chan.sort_values('diam_um',ascending=False).head(15)[['type','instance','superclass','subclass','diam_um','area_um2','dimorphism']].to_string())
print(json.dumps(out,indent=1)); json.dump(out,open(f'{R}/budget_male.json','w'),indent=1)
chan.to_parquet(f'{R}/channel_axons_male.parquet')
