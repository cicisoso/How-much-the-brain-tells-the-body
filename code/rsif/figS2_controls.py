"""Figure S2. Robustness of the direction balance and arbor controls for the calibre-output scaling."""
from figlib_rsif import *
import matplotlib.gridspec as gridspec
m=pd.read_parquet(f'{R}/channel_axons_male_arbor.parquet'); K=json.load(open(f'{R}/capacity_sexes.json')); A=json.load(open(f'{R}/arbor_scaling.json'))
dn=m[m.dir=='down']; up=m[m.dir=='up']
fig=plt.figure(figsize=(W,112*MM)); gs=gridspec.GridSpec(2,2,wspace=0.42,hspace=0.6,left=0.09,right=0.93,top=0.9,bottom=0.1)
axA,axB,axC,axD=[fig.add_subplot(gs[i,j]) for i in range(2) for j in range(2)]
# (a) bootstrap distribution of the capacity ratio at beta = 1 (male); 2000 resamples within each direction
ch=pd.read_parquet(f'{R}/channel_axons_male.parquet'); di=ch[ch.dir=='down'].diam_um.values; ui=ch[ch.dir=='up'].diam_um.values
rng=np.random.default_rng(0); bs=np.array([rng.choice(di,len(di)).sum()/rng.choice(ui,len(ui)).sum() for i in range(2000)])
axA.hist(bs,bins=40,color=C['DN_light'],lw=0); axA.axvline(np.median(bs),color='k',lw=0.8); lo,hi=K['male']['ratio_beta1_ci']
axA.set_xlabel('descending / ascending capacity ratio (β = 1)'); axA.set_ylabel('bootstrap samples'); axA.set_xlim(0.7,0.9)
axA.text(0.97,0.96,f'2000 resamples\n95% CI {lo:.2f}–{hi:.2f}',transform=axA.transAxes,va='top',ha='right')
# (b) output synapses per square micrometre of neck cross-section by diameter quartile
pos=np.arange(4); rng=np.random.default_rng(1)
for sc,col,color,off in [('descending_neuron','n_out_vnc',C['DN'],-0.19),('ascending_neuron','n_out_brain',C['AN'],0.19)]:
    d=m[(m.superclass==sc)&(m.diam_um>=0.1)].copy(); d['q']=pd.qcut(d.diam_um,4,labels=False); d['v']=d[col]/d.area_um2; d=d[d.v>0]
    boxes(axB,[d[d.q==i].v.values for i in range(4)],pos+off,color,s=1.5,palpha=0.3,rng=rng)
    axB.plot([],[],color=color,lw=3,alpha=0.5,label={'descending_neuron':'DN (outputs in VNC)','ascending_neuron':'AN (outputs in brain)'}[sc])
axB.set_yscale('log'); plain_log(axB,'y'); axB.set_ylim(200,3e5); axB.set_xticks(pos); axB.set_xticklabels(['thinnest','2nd','3rd','thickest']); axB.set_xlabel('diameter quartile'); axB.set_ylabel('outputs per µm² of neck cross-section')
axB.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,handlelength=1.0,borderaxespad=0.1)
# (c) partial regression: diameter still predicts output at fixed cable in the target region
for sc,S_,L_,color,lab in [('descending_neuron','n_out_vnc','cable_vnc_um',C['DN'],'DN'),('ascending_neuron','n_out_brain','cable_brain_um',C['AN'],'AN')]:
    d=m[(m.superclass==sc)&(m.diam_um>=0.1)&(m[S_]>0)&(m[L_]>0)&(m['tips_vnc' if lab=='DN' else 'tips_brain']>0)]
    ly=np.log10(d[S_].values.astype(float)); lx=np.log10(d.diam_um.values); lL=np.log10(d[L_].values)
    ry=ly-np.polyval(np.polyfit(lL,ly,1),lL); rx=lx-np.polyval(np.polyfit(lL,lx,1),lL)
    axC.scatter(rx,ry,s=2.5,color=color,alpha=0.25,lw=0,rasterized=True); b=np.polyfit(rx,ry,1); xx=np.array([rx.min(),rx.max()]); axC.plot(xx,np.polyval(b,xx),color=color,lw=1.2)
    mm=A[sc]['model_d_L']['logd']; axC.plot([],[],color=color,lw=1.2,label=f'{lab}: partial exponent {mm["coef"]:.2f} {ci(*mm["ci"])}')
axC.axhline(0,color=C['grey'],lw=0.4,ls=':'); axC.axvline(0,color=C['grey'],lw=0.4,ls=':')
axC.set_xlabel('log diameter given cable (residual)'); axC.set_ylabel('log output given cable (residual)')
axC.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),handlelength=1.2,labelspacing=0.2,borderaxespad=0.1)
# (d) sensory afferents: output follows the arbor, not the calibre
d=m[(m.superclass=='sensory_ascending')&(m.diam_um>=0.1)&(m.n_out_brain>0)&(m.cable_brain_um>0)&(m.tips_brain>0)]
sc_=axD.scatter(d.cable_brain_um,d.n_out_brain,s=5,c=np.log10(d.diam_um),cmap='viridis',vmin=np.log10(0.15),vmax=np.log10(1.5),lw=0,rasterized=True)
fit_line(axD,d.cable_brain_um.values,d.n_out_brain.values.astype(float),C['SA'])
s=A['sensory_ascending']['S_vs_L']['slope']; md=A['sensory_ascending']['model_d_L']['logd']
axD.set_xscale('log'); axD.set_yscale('log'); plain_log(axD); axD.set_xlabel('SA cable in the brain (µm)'); axD.set_ylabel('SA output synapses in the brain')
axD.text(0.97,0.03,f'n = {len(d)}\nslope on cable {s["coef"]:.2f} {ci(*s["ci"])}\ndiameter at fixed cable {md["coef"]:.2f} {ci(*md["ci"])}'.replace('-','−'),transform=axD.transAxes,va='bottom',ha='right')
cb=fig.colorbar(sc_,ax=axD,fraction=0.05,pad=0.02); cb.set_label('neck diameter (µm)'); cb.set_ticks(np.log10([0.2,0.5,1])); cb.set_ticklabels(['0.2','0.5','1']); cb.outline.set_linewidth(0.4)
panel_labels(fig,[axA,axB,axC,axD],'abcd',rows=['ab','cd'])
export(fig,'FigureS2',axes=[axA,axB,axC,axD],ids=list('abcd'),exemptions=[{'panels':['d'],'checks':['panel-width','column','horizontal-gutter'],'reason':'panel d carries a colour bar'}])
