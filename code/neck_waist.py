import os
DATA_DIR=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
import pandas as pd, numpy as np
D=DATA_DIR
neck=pd.read_csv(f'{D}/neck_neurons.csv')
rows=[]
for b,sc in zip(neck.bodyId,neck.superclass):
    d=pd.read_csv(f'{D}/swc/{b}.swc',sep=r'\s+',comment='#',header=None,names=['id','t','x','y','z','r','p'])
    s=d[(d.z>40000)&(d.z<70000)].copy(); s['b']=b; s['sc']=sc; rows.append(s)
P=pd.concat(rows); P['zb']=(P.z//500)*500
g=P.groupby('zb').agg(n=('b','nunique'),sx=('x','std'),sy=('y','std'),xmin=('x',lambda v:v.quantile(.01)),xmax=('x',lambda v:v.quantile(.99)),ymin=('y',lambda v:v.quantile(.01)),ymax=('y',lambda v:v.quantile(.99)))
g['spread']=np.sqrt(g.sx**2+g.sy**2)
pd.set_option('display.width',200); print(g.to_string())
