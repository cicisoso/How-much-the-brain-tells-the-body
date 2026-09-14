"""Per-neuron synapse counts by ROI for neck neurons from the syn-partners table (streamed)."""
import os
DATA_DIR=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
import pyarrow as pa, pyarrow.feather as pf, pyarrow.ipc as ipc, pandas as pd, numpy as np, time
D=DATA_DIR
neck=pd.read_csv(f'{D}/neck_neurons.csv'); ids=set(neck.bodyId.astype('int64'))
t=time.time(); rf=ipc.open_file(pa.memory_map(f'{D}/syn-partners-male-cns-v1.0-minconf-0.5.feather','r'))
print('batches',rf.num_record_batches, rf.schema, flush=True)
pre_rows=[]; post_rows=[]
for i in range(rf.num_record_batches):
    b=rf.get_batch(i).to_pandas()
    cols=b.columns.tolist()
    if i==0: print(cols,flush=True)
    m1=b[b.body_pre.isin(ids)] if 'body_pre' in cols else b[b.pre_body.isin(ids)]
    m2=b[b.body_post.isin(ids)] if 'body_post' in cols else b[b.post_body.isin(ids)]
    pre_rows.append(m1); post_rows.append(m2)
    if i%20==0: print(i,round(time.time()-t),'s',flush=True)
pre=pd.concat(pre_rows); post=pd.concat(post_rows)
pre.to_parquet(f'{D}/neck_syn_pre.parquet'); post.to_parquet(f'{D}/neck_syn_post.parquet'); print('done',len(pre),len(post),round(time.time()-t),'s')
