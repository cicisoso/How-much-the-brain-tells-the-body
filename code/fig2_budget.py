from figlib import *
# publication settings (also applied in figlib): Arial sans-serif, editable text, alignment gate before export
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans','Liberation Sans']
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})

import matplotlib.gridspec as gridspec, json
R='results'; ch=pd.read_parquet(f'{R}/channel_axons_male.parquet'); B=json.load(open(f'{R}/budget_male.json'))
fig=plt.figure(figsize=(174*MM,60*MM))
gs=gridspec.GridSpec(1,4,width_ratios=[1.35,1,1.1,1],wspace=0.6,left=0.06,right=0.985,top=0.84,bottom=0.24)
axA,axB,axC,axD=[fig.add_subplot(gs[0,i]) for i in range(4)]
# A: caliber distributions on log axis
bins=np.logspace(np.log10(0.08),np.log10(12),45)
for sc in ['ascending_neuron','sensory_ascending','descending_neuron']:
    d=ch[ch.superclass==sc].diam_um; short={'ascending_neuron':'AN','sensory_ascending':'SA','descending_neuron':'DN'}[sc]; axA.hist(d,bins=bins,histtype='step',lw=1.0,color=CLASS_COLOR[sc],label=f'{short} (n = {len(d):,})')
    axA.axvline(d.median(),color=CLASS_COLOR[sc],lw=0.6,ls=':')
axA.set_xscale('log'); plain_log(axA,'x'); axA.set_xlabel('axon diameter at neck (µm)'); axA.set_ylabel('axons'); axA.set_ylim(0,205); axA.legend(loc='upper right',fontsize=5.5,handlelength=1.0,borderaxespad=0.2,labelspacing=0.3)
top=ch.sort_values('diam_um',ascending=False).drop_duplicates('type').head(2)
for _,r in top.iterrows(): axA.annotate(r.type.replace('DNp01','DNp01 (GF)'),xy=(r.diam_um,2),xytext=(r.diam_um,55 if r.type=='DNp01' else 30),fontsize=5,ha='center',rotation=90,va='bottom',arrowprops=dict(arrowstyle='-',lw=0.4,color=C['grey']))
# B: Lorenz curve of cross-sectional area
a=np.sort(ch.area_um2.values)[::-1]; cs=np.cumsum(a)/a.sum(); x=np.arange(1,len(a)+1)/len(a)
axB.plot(x,cs,color=C['DN'],lw=1.2); axB.plot([0,1],[0,1],color=C['grey'],lw=0.5,ls='--')
n50=B['n_axons_for_half_area']; axB.plot([n50/len(a)],[0.5],'o',ms=3,color=C['AN']); axB.annotate(f'{n50} axons\n= 50% of area',xy=(n50/len(a),0.5),xytext=(0.42,0.18),fontsize=5.5,arrowprops=dict(arrowstyle='-',lw=0.4,color=C['grey']))
axB.set_xlabel('fraction of axons (thickest first)'); axB.set_ylabel('cumulative area fraction'); axB.set_xlim(0,1); axB.set_ylim(0,1); axB.text(0.6,0.12,f'Gini = {B["gini_area"]:.2f}',fontsize=6)
# C: down vs up under three measures
meas=['axon\ncount','Σ diameter\n(r ~ d)','Σ area\n(r ~ d²)']
dn=ch[ch.dir=='down']; up=ch[ch.dir=='up']
vals_dn=[len(dn),dn.diam_um.sum(),dn.area_um2.sum()]; vals_up=[len(up),up.diam_um.sum(),up.area_um2.sum()]
frac_dn=np.array(vals_dn)/(np.array(vals_dn)+np.array(vals_up))
xx=np.arange(3); axC.bar(xx,frac_dn,color=C['DN'],width=0.6,label='brain → body'); axC.bar(xx,1-frac_dn,bottom=frac_dn,color=C['AN'],width=0.6,label='body → brain')
for i,f in enumerate(frac_dn): axC.text(i,f/2,f'{100*f:.0f}%',ha='center',va='center',color='white',fontsize=6); axC.text(i,f+(1-f)/2,f'{100*(1-f):.0f}%',ha='center',va='center',color='white',fontsize=6)
axC.set_xticks(xx); axC.set_xticklabels(meas,fontsize=5.5); axC.set_ylabel('share of channel'); axC.set_ylim(0,1); axC.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=1,fontsize=5.5,handlelength=1,columnspacing=0.8,borderaxespad=0.1)
# D: ratio vs beta
betas=np.linspace(0.25,2.5,50); ratio=[(dn.diam_um**b).sum()/(up.diam_um**b).sum() for b in betas]
axD.axvspan(0.8,1.2,color=C['light'],lw=0); axD.plot(betas,ratio,color='k',lw=1.0); axD.axhline(1,color=C['grey'],lw=0.5,ls='--')
axD.set_xlabel('exponent β (rate $\\propto$ d$^{\\beta}$)',fontsize=7.5); axD.set_ylabel('descending / ascending\ncapacity ratio'); axD.set_ylim(0.4,2.0)
axD.text(1.0,2.03,'empirical range',ha='center',va='bottom',fontsize=5.5,color=C['grey'],clip_on=False); axD.text(1.75,0.95,'parity',fontsize=5.5,color=C['grey'],va='top')
panel_labels_tight(fig,[axA,axB,axC,axD],'ABCD')
base='figures/Figure2'; require_matplotlib_panel_alignment(fig,json_out=base+'.alignment.json',overlay_svg=base+'.alignment.svg',tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=False,strict=True,axes=[axA,axB,axC,axD],panel_ids=['A','B','C','D'],exemptions=[{'panels':['A','C'],'checks':['panel-width'],'reason':'distribution panel and stacked-bar panel intentionally wider'}])
fig.savefig(base+'.pdf'); fig.savefig(base+'.svg'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300)
