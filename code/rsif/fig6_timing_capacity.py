"""Figure 6. What the cable buys in time, and how the capacity balance depends on the rate-calibre exponent."""
from figlib_rsif import *
import matplotlib.gridspec as gridspec
dl=pd.read_parquet(f'{R}/channel_axons_male_delay.parquet'); K=json.load(open(f'{R}/capacity_sexes.json')); m=pd.read_parquet(f'{R}/channel_axons_male_connectivity.parquet')
fig=plt.figure(figsize=(W,64*MM))
gs=gridspec.GridSpec(1,3,wspace=0.5,left=0.075,right=0.985,top=0.8,bottom=0.2)
axA,axB,axC=[fig.add_subplot(gs[0,i]) for i in range(3)]
d=dl[(dl.end_brain=='branch')&(dl.end_vnc=='branch')]; bins=np.linspace(0,2.5,40)
for sc in ['descending_neuron','ascending_neuron']:
    x=d[d.superclass==sc].latency_ms; axA.hist(x,bins=bins,histtype='step',color=CLASS_COLOR[sc],lw=1.0,label={'descending_neuron':'DN','ascending_neuron':'AN'}[sc]); axA.axvline(x.median(),color=CLASS_COLOR[sc],lw=0.6,ls=':')
gf=dl[dl.type=='DNp01'].latency_ms.max(); axA.plot([gf,gf],[10,205],color=C['grey'],lw=0.5); axA.text(gf,207,'GF',ha='center',va='bottom',color=C['DN'])
axA.set_xlabel('arbor-to-arbor conduction time (ms)'); axA.set_ylabel('axons'); axA.set_ylim(0,230); axA.legend(handlelength=1.2,loc='upper right',labelspacing=0.2)
b=np.array(K['betas']); axB.axvspan(0.8,1.2,color=C['light'],lw=0)
for sex,ls in [('male','-'),('female','--')]:
    y=np.array(K[sex]['ratio_curve']); k=y<=2.0; axB.plot(b[k],y[k],color=C[sex],lw=1.1,ls=ls,label=sex); axB.axhline(1,color=C['grey'],lw=0.5,ls=':')
for sex,y in [('male',1.0),('female',0.93)]:
    lo,hi=K[sex]['parity_beta_ci']; axB.plot([lo,hi],[y,y],color=C[sex],lw=2.0,solid_capstyle='butt')
axB.set_xlabel('rate–calibre exponent β'); axB.set_ylabel('descending / ascending\ncapacity ratio'); axB.set_ylim(0.4,2.0); axB.set_xlim(0.25,2.5)
axB.text(1.0,1.9,'empirical\nrange',ha='center',va='top',color=C['grey']); axB.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,handlelength=1.6,borderaxespad=0.1)
order=m.sort_values('diam_um',ascending=False); ks=np.concatenate([[0],np.unique(np.round(np.logspace(0,3,40)).astype(int))]); rat=[]; sh=[]
for k in ks:
    rest=order.iloc[k:]; dd=rest[rest.dir=='down']; uu=rest[rest.dir=='up']; rat.append(dd.diam_um.sum()/uu.diam_um.sum()); sh.append(dd.area_um2.sum()/rest.area_um2.sum())
axC.plot(ks+1,rat,color='k',lw=1.1,label='capacity ratio, β = 1'); axC.plot(ks+1,sh,color=C['DN'],lw=1.1,label='descending area share'); axC.axhline(0.5,color=C['grey'],lw=0.5,ls=':')
axC.set_xscale('log'); axC.set_xticks([1,11,101,1001]); axC.set_xticklabels(['0','10','100','1000']); axC.set_xlabel('thickest axons removed'); axC.set_ylabel('value'); axC.set_ylim(0.2,1.0)
axC.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=1,handlelength=1.2,labelspacing=0.2,borderaxespad=0.1)
panel_labels(fig,[axA,axB,axC],'abc',rows=['abc'])
export(fig,'Figure6',axes=[axA,axB,axC],ids=list('abc'))
