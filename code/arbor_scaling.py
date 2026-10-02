"""Is the caliber-synapse scaling more than 'big neurons are big everywhere'?
For every brain-body axon of the MaleCNS, measure from its skeleton the cable length, terminal tips and branch points that it
builds in the brain and in the VNC, then (1) decompose the synapse-caliber exponent into a cable term and a density term,
(2) ask whether neck diameter still predicts synaptic output once arbor cable (or whole-neuron cable) is controlled, and
(3) compare the scaling of terminal tips with diameter with the exponents expected from branching rules
(Rall 3/2 impedance matching, area preservation 2, Murray's law 3).
Writes results/arbor_scaling.parquet and results/arbor_scaling.json."""
import os, json, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy import stats
import statsmodels.api as sm
DATA_DIR=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
D=DATA_DIR; R='results'
VOX=0.008                       # um per voxel (skeleton coordinates are in 8 nm voxels)
Z_BRAIN, Z_VNC = 51500, 56500   # connective proper (minimal lateral spread) lies between these z planes; brain anterior (low z)

def arbor(b):
    p=f'{D}/swc/{b}.swc'
    if not os.path.exists(p): return None
    d=pd.read_csv(p,sep=r'\s+',comment='#',header=None,names=['id','t','x','y','z','r','p'])
    ids=d.id.values; par=d.p.values; idx=pd.Series(np.arange(len(d)),index=ids)
    has=par!=-1; ch=np.where(has)[0]; pa=idx.reindex(par[has]).values
    ok=~np.isnan(pa); ch=ch[ok]; pa=pa[ok].astype(int)
    xyz=d[['x','y','z']].values.astype(float)
    seg=np.linalg.norm((xyz[ch]-xyz[pa])*VOX,axis=1); zm=(xyz[ch,2]+xyz[pa,2])/2
    deg=np.bincount(np.concatenate([ch,pa]),minlength=len(d))
    tip=deg==1; bp=deg>=3; z=xyz[:,2]; rad=d.r.values*VOX
    mr=lambda msk: float(rad[msk].mean()) if msk.any() else np.nan   # skeleton radii are floored at 32 voxels (0.256 um); used only for rank correlations
    out=dict(body=int(b),cable_total_um=float(seg.sum()),skel_diam_brain_um=2*mr(z<Z_BRAIN),skel_diam_vnc_um=2*mr(z>Z_VNC),
             cable_brain_um=float(seg[zm<Z_BRAIN].sum()),cable_vnc_um=float(seg[zm>Z_VNC].sum()),cable_neck_um=float(seg[(zm>=Z_BRAIN)&(zm<=Z_VNC)].sum()),
             tips_brain=int((tip&(z<Z_BRAIN)).sum()),tips_vnc=int((tip&(z>Z_VNC)).sum()),
             bp_brain=int((bp&(z<Z_BRAIN)).sum()),bp_vnc=int((bp&(z>Z_VNC)).sum()),n_nodes=int(len(d)))
    return out

rng=np.random.default_rng(0)
def boot_ols(y,X,names,n=1000):
    """OLS of y on X (with constant); returns coefficients with percentile bootstrap 95% CIs (resampling axons)."""
    Xc=sm.add_constant(X); r=sm.OLS(y,Xc).fit(); co=np.asarray(r.params); pv=np.asarray(r.pvalues); bs=[]
    for k in range(n):
        i=rng.integers(0,len(y),len(y)); bs.append(np.asarray(sm.OLS(y[i],Xc[i]).fit().params))
    bs=np.array(bs); res={'n':int(len(y)),'R2':round(float(r.rsquared),3)}
    for j,nm in enumerate(['const']+names):
        res[nm]={'coef':round(float(co[j]),3),'ci':[round(float(v),3) for v in np.percentile(bs[:,j],[2.5,97.5])],'p':float(f'{pv[j]:.2g}')}
    return res

def slope(x,y,n=1000):
    lx=np.log10(x); ly=np.log10(y); return boot_ols(ly,lx[:,None],['slope'],n)

if __name__=='__main__':
    m=pd.read_parquet(f'{R}/channel_axons_male_connectivity.parquet')
    if os.path.exists(f'{R}/arbor_scaling.parquet') and os.environ.get('REBUILD')!='1':
        A=pd.read_parquet(f'{R}/arbor_scaling.parquet')
    else:
        with Pool(16) as pool: rows=pool.map(arbor,m.body.astype('int64').tolist(),chunksize=20)
        A=pd.DataFrame([r for r in rows if r]); A.to_parquet(f'{R}/arbor_scaling.parquet')
    x=m.merge(A,on='body',how='inner'); x=x[x.diam_um>=0.1]
    out={'n_skeletons':int(len(A)),'boundaries_z_vox':[Z_BRAIN,Z_VNC]}
    # direction-specific target region: DN outputs in VNC (inputs in brain), AN/SA outputs in brain (inputs in VNC)
    cfg={'descending_neuron':dict(out='n_out_vnc',Lout='cable_vnc_um',Tout='tips_vnc',inp='n_in_brain',Lin='cable_brain_um',Tin='tips_brain'),
         'ascending_neuron': dict(out='n_out_brain',Lout='cable_brain_um',Tout='tips_brain',inp='n_in_vnc',Lin='cable_vnc_um',Tin='tips_vnc'),
         'sensory_ascending':dict(out='n_out_brain',Lout='cable_brain_um',Tout='tips_brain',inp='n_in_vnc',Lin='cable_vnc_um',Tin='tips_vnc')}
    for sc,c in cfg.items():
        d=x[(x.superclass==sc)&(x[c['out']]>0)&(x[c['Lout']]>0)&(x[c['Tout']]>0)].copy()
        dd=d.diam_um.values; S=d[c['out']].values.astype(float); L=d[c['Lout']].values; T=d[c['Tout']].values.astype(float); Ltot=d.cable_total_um.values
        r={'n':int(len(d))}
        r['S_vs_d']=slope(dd,S); r['L_vs_d']=slope(dd,L); r['T_vs_d']=slope(dd,T); r['Ltot_vs_d']=slope(dd,Ltot)
        r['density_vs_d']=slope(dd,S/L)            # synapses per um of target-region cable
        r['S_vs_L']=slope(L,S); r['S_vs_T']=slope(T,S); r['L_per_tip_vs_d']=slope(dd,L/T)
        ly=np.log10(S)
        r['model_d_L']=boot_ols(ly,np.column_stack([np.log10(dd),np.log10(L)]),['logd','logL'])
        r['model_L_only']=boot_ols(ly,np.log10(L)[:,None],['logL'],n=200)
        r['model_d_Ltot']=boot_ols(ly,np.column_stack([np.log10(dd),np.log10(Ltot)]),['logd','logLtot'])
        r['model_Ltot_only']=boot_ols(ly,np.log10(Ltot)[:,None],['logLtot'],n=200)
        r['model_d_T']=boot_ols(ly,np.column_stack([np.log10(dd),np.log10(T)]),['logd','logT'])
        # partial Spearman correlation of diameter and output given target-region cable (rank residuals)
        rk=lambda v: stats.rankdata(v)
        def resid(a,b): bb=np.polyfit(b,a,1); return a-np.polyval(bb,b)
        pr=stats.pearsonr(resid(rk(dd),rk(L)),resid(rk(S),rk(L)))
        r['partial_spearman_d_S_given_L']={'r':round(float(pr[0]),3),'p':float(f'{pr[1]:.2g}')}
        # quartiles of diameter: median cable, tips, density
        q=pd.qcut(d.diam_um,4,labels=['Q1','Q2','Q3','Q4'])
        r['quartiles']=d.assign(q=q,dens=S/L).groupby('q',observed=True).agg(n=('body','size'),d_med=('diam_um','median'),L_med=(c['Lout'],'median'),T_med=(c['Tout'],'median'),S_med=(c['out'],'median'),dens_med=('dens','median')).round(3).to_dict('index')
        # inputs on the other side (dendritic field) for DN/AN
        if sc!='sensory_ascending':
            e=x[(x.superclass==sc)&(x[c['inp']]>0)&(x[c['Lin']]>0)]
            r['Sin_vs_d']=slope(e.diam_um.values,e[c['inp']].values.astype(float)); r['Lin_vs_d']=slope(e.diam_um.values,e[c['Lin']].values)
            r['density_in_vs_d']=slope(e.diam_um.values,(e[c['inp']]/e[c['Lin']]).values)
            r['model_in_d_L']=boot_ols(np.log10(e[c['inp']].values.astype(float)),np.column_stack([np.log10(e.diam_um.values),np.log10(e[c['Lin']].values)]),['logd','logL'])
        out[sc]=r
    # pooled DN+AN central axons: tips vs diameter (branching-rule test) with both arbors summed
    ca=x[x.superclass.isin(['descending_neuron','ascending_neuron'])].copy(); ca['T_target']=np.where(ca.superclass=='descending_neuron',ca.tips_vnc,ca.tips_brain)
    ca=ca[ca.T_target>0]; out['central_T_target_vs_d']=slope(ca.diam_um.values,ca.T_target.values.astype(float))
    # density = presynaptic sites (T-bars) per um x postsynaptic contacts per T-bar (polyadicity)
    VNC_PREFIX=('LegNp','IntTct','LTct','WTct','HTct','NTct','ANm','mVAC','Ov','VNC-unspecified','DMetaN','MesoAN','MesoLN','ProAN','ProLN','MetaAN','MetaLN','CV-unspecified','AbN')
    pre=pd.read_parquet(f'{D}/neck_syn_pre.parquet',columns=['body_pre','x_pre','y_pre','z_pre','primary_post'])
    rmap={c:('VNC' if str(c).startswith(VNC_PREFIX) else 'brain') for c in pre.primary_post.cat.categories}; pre['reg']=pre.primary_post.map(rmap).astype(str)
    tb=pre.drop_duplicates(['body_pre','x_pre','y_pre','z_pre','reg']).groupby(['body_pre','reg']).size().unstack(fill_value=0)
    x=x.merge(tb.rename(columns={'VNC':'tbar_vnc','brain':'tbar_brain'}),left_on='body',right_index=True,how='left'); x[['tbar_vnc','tbar_brain']]=x[['tbar_vnc','tbar_brain']].fillna(0)
    for sc,S_,L_,T_ in [('descending_neuron','n_out_vnc','cable_vnc_um','tbar_vnc'),('ascending_neuron','n_out_brain','cable_brain_um','tbar_brain'),('sensory_ascending','n_out_brain','cable_brain_um','tbar_brain')]:
        d=x[(x.superclass==sc)&(x[S_]>0)&(x[L_]>0)&(x[T_]>0)]
        A_=('skel_diam_vnc_um' if T_=='tbar_vnc' else 'skel_diam_brain_um')
        out[sc]['arbor_calibre']={'rho_neck_d_vs_arbor_skel_d':round(float(stats.spearmanr(d.diam_um,d[A_])[0]),3),'rho_tbar_density_vs_arbor_skel_d':round(float(stats.spearmanr(d[A_],d[T_]/d[L_])[0]),3),'median_arbor_skel_d_um':round(float(d[A_].median()),3)}
        out[sc]['tbar']={'n':int(len(d)),'tbar_vs_d':slope(d.diam_um.values,d[T_].values.astype(float)),'tbar_per_um_vs_d':slope(d.diam_um.values,(d[T_]/d[L_]).values),
                         'contacts_per_tbar_vs_d':slope(d.diam_um.values,(d[S_]/d[T_]).values),'median_contacts_per_tbar':round(float((d[S_]/d[T_]).median()),2),'median_tbar_per_um':round(float((d[T_]/d[L_]).median()),3)}
        # exact decomposition on one sample: log S = log L + log(T/L) + log(S/T), so the OLS slopes on log d add up
        dd_=d.diam_um.values
        dec={'S':slope(dd_,d[S_].values.astype(float)),'L':slope(dd_,d[L_].values),'T_per_L':slope(dd_,(d[T_]/d[L_]).values),'S_per_T':slope(dd_,(d[S_]/d[T_]).values)}
        dec['sum_of_parts']=round(sum(dec[k]['slope']['coef'] for k in ['L','T_per_L','S_per_T']),3)
        q=pd.qcut(d.diam_um,4,labels=['Q1','Q2','Q3','Q4'])
        dec['quartiles']=d.assign(q=q,L=d[L_],T=d[T_],tpl=d[T_]/d[L_],spt=d[S_]/d[T_],S=d[S_]).groupby('q',observed=True).agg(d_med=('diam_um','median'),L_med=('L','median'),T_med=('T','median'),tbar_per_um_med=('tpl','median'),contacts_per_tbar_med=('spt','median'),S_med=('S','median')).round(3).to_dict('index')
        out[sc]['decomposition']=dec
    # target class, cable and diameter together (DNs): does diameter survive both?
    d=x[(x.superclass=='descending_neuron')&(x.n_out_vnc>0)&(x.cable_vnc_um>0)&x.subclass.notna()].copy()
    X=pd.get_dummies(d.subclass,drop_first=True).astype(float); X['logL']=np.log10(d.cable_vnc_um); X['logd']=np.log10(d.diam_um); X=sm.add_constant(X); y=np.log10(d.n_out_vnc)
    r2=sm.OLS(y,X).fit(); r1=sm.OLS(y,X.drop(columns='logd')).fit()
    out['dn_model_subclass_L_d']={'n':int(len(d)),'coef_logd':round(float(r2.params['logd']),3),'ci_logd':[round(float(v),3) for v in r2.conf_int().loc['logd']],'p_logd':float(f'{r2.pvalues["logd"]:.2g}'),
                                  'coef_logL':round(float(r2.params['logL']),3),'R2_with_d':round(float(r2.rsquared),3),'R2_subclass_L':round(float(r1.rsquared),3)}
    # saved table keeps every channel axon (the d >= 0.1 um filter applies to the fits only)
    full=m.merge(A,on='body',how='left').merge(tb.rename(columns={'VNC':'tbar_vnc','brain':'tbar_brain'}),left_on='body',right_index=True,how='left'); full[['tbar_vnc','tbar_brain']]=full[['tbar_vnc','tbar_brain']].fillna(0)
    full.to_parquet(f'{R}/channel_axons_male_arbor.parquet')
    json.dump(out,open(f'{R}/arbor_scaling.json','w'),indent=1)
    def show(sc):
        r=out[sc]; f=lambda k: f"{r[k]['slope']['coef']:.2f} [{r[k]['slope']['ci'][0]:.2f}, {r[k]['slope']['ci'][1]:.2f}]"
        print(f"\n== {sc} n={r['n']}")
        for k in ['S_vs_d','L_vs_d','density_vs_d','T_vs_d','L_per_tip_vs_d','Ltot_vs_d','S_vs_L','S_vs_T']: print(f'  {k:16s} {f(k)}')
        for k in ['model_d_L','model_d_Ltot','model_d_T']:
            mm=r[k]; ks=[kk for kk in mm if kk not in ('n','R2','const')]
            print(f'  {k:16s} R2={mm["R2"]}  '+'  '.join(f'{kk}={mm[kk]["coef"]:.2f} [{mm[kk]["ci"][0]:.2f},{mm[kk]["ci"][1]:.2f}] p={mm[kk]["p"]}' for kk in ks))
        print(f'  L only R2={r["model_L_only"]["R2"]}; Ltot only R2={r["model_Ltot_only"]["R2"]}; partial rho(d,S|L)={r["partial_spearman_d_S_given_L"]}')
        if 'Sin_vs_d' in r: print(f"  inputs: Sin~d {f('Sin_vs_d')}  Lin~d {f('Lin_vs_d')}  dens_in~d {f('density_in_vs_d')}  model_in d={r['model_in_d_L']['logd']['coef']} {r['model_in_d_L']['logd']['ci']} L={r['model_in_d_L']['logL']['coef']}")
        print('  quartiles',json.dumps(r['quartiles']))
        t=r['tbar']; g=lambda k: f"{t[k]['slope']['coef']:.2f} [{t[k]['slope']['ci'][0]:.2f}, {t[k]['slope']['ci'][1]:.2f}]"
        print(f"  T-bars n={t['n']}: tbar~d {g('tbar_vs_d')}  tbar/um~d {g('tbar_per_um_vs_d')}  contacts/tbar~d {g('contacts_per_tbar_vs_d')}  median contacts/tbar {t['median_contacts_per_tbar']}  median tbar/um {t['median_tbar_per_um']}")
    for sc in cfg: show(sc)
    print('\ncentral tips vs d', out['central_T_target_vs_d']['slope'])
    print('DN subclass+L+d model', out['dn_model_subclass_L_d'])
    for sc in cfg: print(sc,'arbor calibre',out[sc]['arbor_calibre'])
