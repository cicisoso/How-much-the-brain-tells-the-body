from figlib import *
# publication settings (also set in figlib): Arial/Helvetica sans-serif, editable text, alignment gate before export
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans']
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})
try:
    from audit_panel_alignment import require_matplotlib_panel_alignment  # optional QA gate (nature-figure toolkit)
except ImportError:
    pass
import matplotlib.gridspec as gridspec, json
R='results'; m=pd.read_parquet(f'{R}/channel_axons_male_connectivity.parquet'); RB=json.load(open(f'{R}/robustness.json'))
dn=m[m.dir=='down']; up=m[m.dir=='up']
fig=plt.figure(figsize=(180*MM,66*MM)); gs=gridspec.GridSpec(1,3,wspace=0.5,left=0.07,right=0.985,top=0.8,bottom=0.25)
axA,axB,axC=[fig.add_subplot(gs[0,i]) for i in range(3)]
# a: remove the k thickest axons
order=m.sort_values('diam_um',ascending=False)
ks=np.unique(np.round(np.logspace(0,np.log10(1000),40)).astype(int)); ks=np.concatenate([[0],ks])
rat=[];sh=[]
for k in ks:
    rest=order.iloc[k:]; d=rest[rest.dir=='down']; u=rest[rest.dir=='up']
    rat.append(d.diam_um.sum()/u.diam_um.sum()); sh.append(d.area_um2.sum()/rest.area_um2.sum())
axA.plot(ks+1,rat,color='k',lw=1.2,label='capacity ratio, β = 1'); axA.plot(ks+1,sh,color=C['DN'],lw=1.2,label='descending area share'); axA.axhline(0.5,color=C['grey'],lw=0.5,ls='--')
axA.set_xscale('log'); plain_log(axA,'x'); axA.set_xticks([1,11,101,1001]); axA.set_xticklabels(['0','10','100','1,000']); axA.set_xlabel('thickest axons removed'); axA.set_ylabel('value'); axA.set_ylim(0.2,1.0); axA.legend(loc='lower left',fontsize=8,handlelength=1.2,labelspacing=0.2,borderaxespad=0.1)
# b: bootstrap of the ratio at beta=1
rng=np.random.default_rng(0); bs=[]
di=dn.diam_um.values; ui=up.diam_um.values
for i in range(2000): bs.append(rng.choice(di,len(di)).sum()/rng.choice(ui,len(ui)).sum())
bs=np.array(bs); axB.hist(bs,bins=40,color=C['DN_light'],lw=0); lo,hi=np.percentile(bs,[2.5,97.5])
axB.axvline(np.median(bs),color='k',lw=0.8)
axB.set_xlabel('descending / ascending capacity ratio (β = 1)'); axB.set_ylabel('bootstrap samples'); axB.set_xlim(0.7,0.96)
axB.text(0.97,0.96,f'2,000 resamples\n95% CI {lo:.2f}–{hi:.2f}',transform=axB.transAxes,va='top',ha='right',fontsize=8)
# c: output synapses per unit cross-section area by diameter quartile
pos=np.arange(4); rng=np.random.default_rng(1)
for k,(sc,col,color,off) in enumerate([('descending_neuron','n_out_vnc',C['DN'],-0.19),('ascending_neuron','n_out_brain',C['AN'],0.19)]):
    d=m[(m.superclass==sc)].copy(); d['q']=pd.qcut(d.diam_um,4,labels=False); d['v']=d[col]/d.area_um2
    data=[d[d.q==i].v.values for i in range(4)]
    bp=axC.boxplot(data,positions=pos+off,widths=0.32,whis=(5,95),showfliers=False,patch_artist=True,medianprops=dict(color='k',lw=0.8),whiskerprops=dict(lw=0.5),capprops=dict(lw=0.5),boxprops=dict(lw=0.5))
    for patch in bp['boxes']: patch.set_facecolor(color); patch.set_alpha(0.35)
    for i,dd in enumerate(data): axC.scatter(pos[i]+off+rng.uniform(-0.1,0.1,len(dd)),dd,s=1.5,color=color,alpha=0.3,lw=0,rasterized=True,zorder=3)
    axC.plot([],[],color=color,lw=3,alpha=0.5,label={'descending_neuron':'DN (outputs in VNC)','ascending_neuron':'AN (outputs in brain)'}[sc])
axC.set_yscale('log'); plain_log(axC,'y'); axC.set_ylim(200,3e5); axC.set_xticks(pos); axC.set_xticklabels(['thinnest','2nd','3rd','thickest']); axC.set_xlabel('diameter quartile'); axC.set_ylabel('outputs per µm² of cable'); axC.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=1,fontsize=8,handlelength=1.0,labelspacing=0.2,borderaxespad=0.1)
panel_labels_tight(fig,[axA,axB,axC],'abc')
base='figures/FigureS2'; require_matplotlib_panel_alignment(fig,json_out=base+'.alignment.json',overlay_svg=base+'.alignment.svg',tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=False,strict=True)
fig.savefig(base+'.pdf'); fig.savefig(base+'.svg'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300); print('saved',base)
