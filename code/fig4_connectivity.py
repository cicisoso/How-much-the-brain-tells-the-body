from figlib import *
# publication settings (also set in figlib): Arial/Helvetica sans-serif, editable text, alignment gate before export
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans']
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})
try:
    from audit_panel_alignment import require_matplotlib_panel_alignment  # optional QA gate (nature-figure toolkit)
except ImportError:
    pass
import matplotlib.gridspec as gridspec, json
from scipy import stats
R='results'; m=pd.read_parquet(f'{R}/channel_axons_male_connectivity.parquet'); J=json.load(open(f'{R}/connectivity_scaling.json'))
m=m[m.diam_um>=0.1]
fig=plt.figure(figsize=(180*MM,112*MM))
gs=gridspec.GridSpec(2,3,wspace=0.5,hspace=0.55,left=0.075,right=0.985,top=0.93,bottom=0.1)
ax=[fig.add_subplot(gs[i,j]) for i in range(2) for j in range(3)]
xlim=[0.09,12]
def scat(a,d,col,color,label,key):
    a.scatter(d.diam_um,d[col].clip(lower=1),s=3,color=C['grey'],alpha=0.45,lw=0,rasterized=True)
    s=fit_line(a,d.diam_um.values,d[col].values.astype(float),color,xlim=[d.diam_um.min(),d.diam_um.max()])
    j=J[key]; a.text(0.97,0.03,f'{label}, n = {j["n"]:,}\nslope {j["slope"]:.2f} [{j["slope_ci"][0]:.2f}, {j["slope_ci"][1]:.2f}]\nρ = {j["rho"]:.2f}',transform=a.transAxes,va='bottom',ha='right',fontsize=8)
    a.set_xscale('log'); a.set_yscale('log'); plain_log(a); a.set_xlim(xlim); a.set_ylim(8,3e5); a.set_xlabel('neck diameter (µm)')
dn=m[m.superclass=='descending_neuron']; an=m[m.superclass=='ascending_neuron']; sa=m[m.superclass=='sensory_ascending']
scat(ax[0],dn,'n_out_vnc',C['DN'],'DN','descending_neuron:n_out_vnc'); ax[0].set_ylabel('output synapses in VNC')
scat(ax[1],an,'n_out_brain',C['AN'],'AN','ascending_neuron:n_out_brain'); ax[1].set_ylabel('output synapses in brain')
scat(ax[2],sa,'n_out_brain',C['SA'],'SA','sensory_ascending:n_out_brain'); ax[2].set_ylabel('output synapses in brain')
# d: strong partners
a=ax[3]
for d,color,lab,key in [(dn,C['DN'],'DN','descending_neuron:n_out_partners5'),(an,C['AN'],'AN','ascending_neuron:n_out_partners5')]:
    a.scatter(d.diam_um,d.n_out_partners5.clip(lower=1),s=3,color=color,alpha=0.25,lw=0,rasterized=True); fit_line(a,d.diam_um.values,d.n_out_partners5.values.astype(float),color,xlim=[d.diam_um.min(),d.diam_um.max()])
    j=J[key]; a.plot([],[],color=color,lw=1.2,label=f'{lab}: slope {j["slope"]:.2f} [{j["slope_ci"][0]:.2f}, {j["slope_ci"][1]:.2f}]')
a.set_xscale('log'); a.set_yscale('log'); plain_log(a); a.set_xlim(xlim); a.set_ylim(0.8,3e3); a.set_xlabel('neck diameter (µm)'); a.set_ylabel('partners (≥5 synapses)'); a.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=8,handlelength=1.2,labelspacing=0.2,borderaxespad=0.1)
# e: population size of cell types vs median diameter
a=ax[4]; bins_=['1-2','3-4','5-8','>8']; pos=np.arange(4); rng=np.random.default_rng(0)
for k,(sc,color,off) in enumerate([('descending_neuron',C['DN'],-0.19),('ascending_neuron',C['AN'],0.19)]):
    t=pd.read_csv(f'{R}/type_population_{sc}.csv'); data=[t[t['pop']==b].d_med.values for b in bins_]
    bp=a.boxplot(data,positions=pos+off,widths=0.32,whis=(5,95),showfliers=False,patch_artist=True,medianprops=dict(color='k',lw=0.8),whiskerprops=dict(lw=0.5),capprops=dict(lw=0.5),boxprops=dict(lw=0.5))
    for patch in bp['boxes']: patch.set_facecolor(color); patch.set_alpha(0.35)
    for i,dd in enumerate(data): a.scatter(pos[i]+off+rng.uniform(-0.1,0.1,len(dd)),dd,s=2.5,color=color,alpha=0.6,lw=0,rasterized=True,zorder=3)
    a.plot([],[],color=color,lw=3,alpha=0.5,label={'descending_neuron':'DN types','ascending_neuron':'AN types'}[sc])
a.set_yscale('log'); plain_log(a,'y'); a.set_ylim(0.12,12); a.set_xticks(pos); a.set_xticklabels(bins_); a.set_xlabel('neurons per cell type'); a.set_ylabel('median diameter of type (µm)'); a.legend(loc='upper right',fontsize=8,handlelength=1.0,labelspacing=0.2,borderaxespad=0.1)
# f: neurotransmitter
a=ax[5]; nts=['acetylcholine','gaba','glutamate']; ntlab=['ACh','GABA','Glu']
for k,(sc,color,off) in enumerate([('descending_neuron',C['DN'],-0.19),('ascending_neuron',C['AN'],0.19)]):
    d=m[(m.superclass==sc)&m.consensus_nt.isin(nts)]; data=[d[d.consensus_nt==n].diam_um.values for n in nts]
    bp=a.boxplot(data,positions=np.arange(3)+off,widths=0.32,whis=(5,95),showfliers=False,patch_artist=True,medianprops=dict(color='k',lw=0.8),whiskerprops=dict(lw=0.5),capprops=dict(lw=0.5),boxprops=dict(lw=0.5))
    for patch in bp['boxes']: patch.set_facecolor(color); patch.set_alpha(0.35)
    for i,dd in enumerate(data): a.scatter(np.arange(3)[i]+off+rng.uniform(-0.1,0.1,len(dd)),dd,s=1.5,color=color,alpha=0.35,lw=0,rasterized=True,zorder=3)
a.set_yscale('log'); plain_log(a,'y'); a.set_ylim(0.12,12); a.set_xticks(np.arange(3)); a.set_xticklabels(ntlab); a.set_xlabel('predicted neurotransmitter'); a.set_ylabel('neck diameter (µm)')
panel_labels_tight(fig,ax,'abcdef')
base='figures/Figure4'; require_matplotlib_panel_alignment(fig,json_out=base+'.alignment.json',overlay_svg=base+'.alignment.svg',tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=False,strict=True)
fig.savefig(base+'.pdf'); fig.savefig(base+'.svg'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300); print('saved',base)
