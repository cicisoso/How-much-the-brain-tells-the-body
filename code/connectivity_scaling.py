"""Caliber versus synaptic connectivity, cell-type population size and neurotransmitter for every brain-body axon (MaleCNS).
Writes results/channel_axons_male_connectivity.parquet and results/connectivity_scaling.json."""
import os
DATA_DIR=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
import pandas as pd, numpy as np, json
from scipy import stats
import statsmodels.api as sm
D=DATA_DIR; R='results'
ch=pd.read_parquet(f'{R}/channel_axons_male.parquet'); ch['body']=ch.body.astype('int64')
pre=pd.read_parquet(f'{D}/neck_syn_pre.parquet',columns=['body_pre','body_post','primary_post'])
post=pd.read_parquet(f'{D}/neck_syn_post.parquet',columns=['body_pre','body_post','primary_post'])
VNC_PREFIX=('LegNp','IntTct','LTct','WTct','HTct','NTct','ANm','mVAC','Ov','VNC-unspecified','DMetaN','MesoAN','MesoLN','ProAN','ProLN','MetaAN','MetaLN','CV-unspecified','AbN')
def region(r): return 'VNC' if str(r).startswith(VNC_PREFIX) else 'brain'
rmap={c:region(c) for c in pre.primary_post.cat.categories}
pre['reg']=pre.primary_post.map(rmap); post['reg']=post.primary_post.map(rmap)
g=pre.groupby('body_pre'); out=pd.DataFrame({'n_out':g.size(),'n_out_partners':g.body_post.nunique()})
outr=pre.groupby(['body_pre','reg']).size().unstack(fill_value=0); out['n_out_vnc']=outr.get('VNC',0); out['n_out_brain']=outr.get('brain',0)
pc=pre.groupby(['body_pre','body_post']).size(); out['n_out_partners5']=(pc>=5).groupby(level=0).sum()
g=post.groupby('body_post'); inp=pd.DataFrame({'n_in':g.size(),'n_in_partners':g.body_pre.nunique()})
inr=post.groupby(['body_post','reg']).size().unstack(fill_value=0); inp['n_in_vnc']=inr.get('VNC',0); inp['n_in_brain']=inr.get('brain',0)
pc2=post.groupby(['body_post','body_pre']).size(); inp['n_in_partners5']=(pc2>=5).groupby(level=0).sum()
m=ch.merge(out,left_on='body',right_index=True,how='left').merge(inp,left_on='body',right_index=True,how='left')
for c in list(out.columns)+list(inp.columns): m[c]=m[c].fillna(0).astype(int)
nt=pd.read_feather(f'{D}/body-neurotransmitters-male-cns-v1.0.feather')[['body','predicted_nt','predicted_nt_confidence','consensus_nt']]
m=m.merge(nt,on='body',how='left')
ann=pd.read_feather(f'{D}/body-annotations-male-cns-v1.0-minconf-0.5.feather')
pop=ann[ann.superclass.isin(['descending_neuron','ascending_neuron','sensory_ascending'])].groupby('type').size(); m['type_n']=m.type.map(pop)
m.to_parquet(f'{R}/channel_axons_male_connectivity.parquet')
rng=np.random.default_rng(0)
def fit(x,y,n=1000):
    ok=(x>0)&(y>0); lx=np.log10(x[ok].values); ly=np.log10(y[ok].values)
    s,i=np.polyfit(lx,ly,1); bs=[]
    for k in range(n):
        idx=rng.integers(0,len(lx),len(lx)); bs.append(np.polyfit(lx[idx],ly[idx],1)[0])
    rho,p=stats.spearmanr(x[ok],y[ok])
    return dict(n=int(ok.sum()),slope=round(float(s),2),slope_ci=[round(float(v),2) for v in np.percentile(bs,[2.5,97.5])],intercept=round(float(i),3),rho=round(float(rho),3),p=float(f'{p:.2g}'))
out={}
sets={'descending_neuron':['n_out_vnc','n_in_brain','n_out_partners5','n_in_partners5','n_out','n_in'],'ascending_neuron':['n_out_brain','n_in_vnc','n_out_partners5','n_in_partners5','n_out','n_in'],'sensory_ascending':['n_out_brain','n_out','n_in','n_out_partners5']}
for sc,cols in sets.items():
    d=m[(m.superclass==sc)&(m.diam_um>=0.1)]
    for c in cols: out[f'{sc}:{c}']=fit(d.diam_um,d[c].astype(float))
    q=pd.qcut(d.diam_um,4,labels=['Q1','Q2','Q3','Q4'])
    out[f'{sc}:quartiles']=d.assign(q=q).groupby('q',observed=True).agg(n=('body','size'),d_med=('diam_um','median'),out_med=(cols[0],'median'),partners_med=('n_out_partners5','median'),type_n_med=('type_n','median')).round(3).to_dict('index')
# caliber explains VNC output beyond target class (DNs)
d=m[(m.superclass=='descending_neuron')&(m.diam_um>=0.1)&m.subclass.notna()].copy()
X=pd.get_dummies(d.subclass,drop_first=True).astype(float); X['logd']=np.log10(d.diam_um); X=sm.add_constant(X)
y=np.log10(d.n_out_vnc.clip(lower=1)); r1=sm.OLS(y,X).fit(); r0=sm.OLS(y,X.drop(columns='logd')).fit()
out['dn_out_vnc_model']={'n':int(len(d)),'coef_logd':round(float(r1.params['logd']),2),'coef_ci':[round(float(v),2) for v in r1.conf_int().loc['logd']],'p_logd':float(f'{r1.pvalues["logd"]:.2g}'),'R2_with':round(float(r1.rsquared),3),'R2_subclass_only':round(float(r0.rsquared),3)}
# synapses per unit cross-section area by quartile (DN VNC outputs; AN brain outputs)
for sc,c in [('descending_neuron','n_out_vnc'),('ascending_neuron','n_out_brain')]:
    d=m[(m.superclass==sc)&(m.diam_um>=0.1)]; q=pd.qcut(d.diam_um,4,labels=['Q1','Q2','Q3','Q4'])
    out[f'{sc}:out_per_um2_by_quartile']=d.assign(q=q,v=d[c]/d.area_um2).groupby('q',observed=True).v.median().round(0).to_dict()
# population size of cell types
for sc in ['descending_neuron','ascending_neuron']:
    d=m[(m.superclass==sc)&(m.diam_um>=0.1)&m.type.notna()]
    t=d.groupby('type').agg(n=('body','size'),d_med=('diam_um','median'),area=('area_um2','sum'),sumd=('diam_um','sum'))
    t['pop']=pd.cut(t.n,[0,2,4,8,1000],labels=['1-2','3-4','5-8','>8'])
    tab=t.groupby('pop',observed=True).agg(types=('n','size'),axons=('n','sum'),d_med=('d_med','median'),d_q25=('d_med',lambda v:v.quantile(.25)),d_q75=('d_med',lambda v:v.quantile(.75)),area=('area','sum'),sumd=('sumd','sum')).round(3)
    rho,p=stats.spearmanr(t.n,t.d_med)
    out[f'{sc}:population']={'n_types':int(len(t)),'rho_n_vs_median_d':round(float(rho),3),'p':float(f'{p:.2g}'),'table':tab.to_dict('index'),'kruskal_p':float(f'{stats.kruskal(*[t[t["pop"]==k].d_med for k in tab.index]).pvalue:.2g}')}
    t.to_csv(f'{R}/type_population_{sc}.csv')
# neurotransmitter
for sc in ['descending_neuron','ascending_neuron','sensory_ascending']:
    d=m[(m.superclass==sc)&(m.diam_um>=0.1)]; g=d[d.consensus_nt.isin(['acetylcholine','gaba','glutamate'])]
    tab=g.groupby('consensus_nt').agg(n=('body','size'),d_med=('diam_um','median'),d_q25=('diam_um',lambda v:v.quantile(.25)),d_q75=('diam_um',lambda v:v.quantile(.75)),area=('area_um2','sum')).round(3)
    res={'table':tab.to_dict('index')}
    ks=[k for k in ['acetylcholine','gaba','glutamate'] if (g.consensus_nt==k).sum()>=10]
    if len(ks)>=2: res['kruskal_p']=float(f'{stats.kruskal(*[g[g.consensus_nt==k].diam_um for k in ks]).pvalue:.2g}')
    for a in ['gaba','glutamate']:
        if (g.consensus_nt==a).sum()>=10:
            mw=stats.mannwhitneyu(g[g.consensus_nt==a].diam_um,g[g.consensus_nt=='acetylcholine'].diam_um,alternative='two-sided'); res[f'{a}_vs_ach']={'U':float(mw.statistic),'p':float(f'{mw.pvalue:.2g}'),'n':[int((g.consensus_nt==a).sum()),int((g.consensus_nt=='acetylcholine').sum())]}
    out[f'{sc}:nt']=res
# direction totals by synapses
dn=m[m.dir=='down']; up=m[m.dir=='up']
out['synapse_totals']={'down_out_vnc':int(dn.n_out_vnc.sum()),'up_out_brain':int(up.n_out_brain.sum()),'down_in_brain':int(dn.n_in_brain.sum()),'up_in_vnc':int(up.n_in_vnc.sum()),'ratio_out_down_up':round(dn.n_out_vnc.sum()/up.n_out_brain.sum(),3),'down_partners5':int(dn.n_out_partners5.sum()),'up_partners5':int(up.n_out_partners5.sum())}
# SA per organ correlations
sa=m[m.superclass=='sensory_ascending']
out['sa_by_organ']={org:{'n':int(len(g)),'rho_d_out':round(float(stats.spearmanr(g.diam_um,g.n_out)[0]),2),'p':float(f'{stats.spearmanr(g.diam_um,g.n_out)[1]:.2g}'),'d_med':round(float(g.diam_um.median()),2)} for org,g in sa.groupby('subclass') if len(g)>=15}
json.dump(out,open(f'{R}/connectivity_scaling.json','w'),indent=1); print(json.dumps(out,indent=1)[:6000])
