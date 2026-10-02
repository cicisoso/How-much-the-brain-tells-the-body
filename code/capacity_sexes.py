"""Descending/ascending capacity ratio sum(d^beta) as a function of the rate-calibre exponent beta, in the male (MaleCNS)
and the female (BANC) connectomes, with the parity exponent and bootstrap 95% confidence intervals (2000 resamples of axons
within each direction). Writes results/capacity_sexes.json."""
import json, numpy as np, pandas as pd
from scipy.optimize import brentq
R='results'; rng=np.random.default_rng(0)
m=pd.read_parquet(f'{R}/channel_axons_male.parquet'); f=pd.read_parquet(f'{R}/channel_axons_female.parquet')
betas=np.round(np.linspace(0.25,2.5,46),3)
def ratio(dn,up,b): return (dn**b).sum()/(up**b).sum()
def parity(dn,up): return brentq(lambda b: ratio(dn,up,b)-1,0.2,4.0)
out={'betas':betas.tolist()}
for sex,df in [('male',m),('female',f)]:
    dn=df[df.dir=='down'].diam_um.values; up=df[df.dir=='up'].diam_um.values
    bs_par=[]; bs_r1=[]
    for k in range(2000):
        a=rng.choice(dn,len(dn)); u=rng.choice(up,len(up)); bs_par.append(parity(a,u)); bs_r1.append(ratio(a,u,1.0))
    out[sex]={'n_down':int(len(dn)),'n_up':int(len(up)),'ratio_curve':[round(float(ratio(dn,up,b)),4) for b in betas],
              'ratio_beta0.5':round(float(ratio(dn,up,0.5)),3),'ratio_beta1':round(float(ratio(dn,up,1.0)),3),'ratio_beta1_ci':[round(float(v),4) for v in np.percentile(bs_r1,[2.5,97.5])],
              'ratio_beta2':round(float(ratio(dn,up,2.0)),3),'parity_beta':round(float(parity(dn,up)),3),'parity_beta_ci':[round(float(v),4) for v in np.percentile(bs_par,[2.5,97.5])]}
json.dump(out,open(f'{R}/capacity_sexes.json','w'),indent=1)
for s in ['male','female']: print(s,{k:v for k,v in out[s].items() if k!='ratio_curve'})
