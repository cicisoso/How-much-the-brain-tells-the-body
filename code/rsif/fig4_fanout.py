"""Figure 4. Fan-out, cell-type population size and predicted transmitter."""
from figlib_rsif import *
import matplotlib.gridspec as gridspec
m=pd.read_parquet(f'{R}/channel_axons_male_connectivity.parquet'); J=json.load(open(f'{R}/connectivity_scaling.json')); m=m[m.diam_um>=0.1]
dn=m[m.superclass=='descending_neuron']; an=m[m.superclass=='ascending_neuron']
fig=plt.figure(figsize=(W,62*MM))
gs=gridspec.GridSpec(1,3,wspace=0.5,left=0.075,right=0.985,top=0.8,bottom=0.17)
ax=[fig.add_subplot(gs[0,i]) for i in range(3)]; xlim=[0.09,12]; rng=np.random.default_rng(0)
a=ax[0]
for d,color,lab,key in [(dn,C['DN'],'DN','descending_neuron:n_out_partners5'),(an,C['AN'],'AN','ascending_neuron:n_out_partners5')]:
    a.scatter(d.diam_um,d.n_out_partners5.clip(lower=1),s=3,color=color,alpha=0.25,lw=0,rasterized=True); fit_line(a,d.diam_um.values,d.n_out_partners5.values.astype(float),color,xlim=[d.diam_um.min(),d.diam_um.max()])
    j=J[key]; a.plot([],[],color=color,lw=1.2,label=f'{lab}: slope {j["slope"]:.2f} {ci(*j["slope_ci"])}')
a.set_xscale('log'); a.set_yscale('log'); plain_log(a); a.set_xlim(xlim); a.set_ylim(0.8,3e3); a.set_xlabel('neck diameter (µm)'); a.set_ylabel('partners (≥ 5 synapses)')
a.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),handlelength=1.2,labelspacing=0.2,borderaxespad=0.1)
a=ax[1]; bins_=['1-2','3-4','5-8','>8']; pos=np.arange(4)
for sc,color,off in [('descending_neuron',C['DN'],-0.19),('ascending_neuron',C['AN'],0.19)]:
    t=pd.read_csv(f'{R}/type_population_{sc}.csv'); data=[t[t['pop']==b].d_med.values for b in bins_]
    boxes(a,data,pos+off,color,rng=rng); a.plot([],[],color=color,lw=3,alpha=0.5,label={'descending_neuron':'DN types','ascending_neuron':'AN types'}[sc])
a.set_yscale('log'); plain_log(a,'y'); a.set_ylim(0.12,12); a.set_xticks(pos); a.set_xticklabels(['1–2','3–4','5–8','> 8']); a.set_xlabel('neurons per cell type'); a.set_ylabel('median diameter of type (µm)')
a.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,handlelength=1.0,labelspacing=0.2,borderaxespad=0.1)
a=ax[2]; nts=['acetylcholine','gaba','glutamate']
for sc,color,off in [('descending_neuron',C['DN'],-0.19),('ascending_neuron',C['AN'],0.19)]:
    d=m[(m.superclass==sc)&m.consensus_nt.isin(nts)]; data=[d[d.consensus_nt==n].diam_um.values for n in nts]
    boxes(a,data,np.arange(3)+off,color,s=1.5,palpha=0.35,rng=rng); a.plot([],[],color=color,lw=3,alpha=0.5,label={'descending_neuron':'DNs','ascending_neuron':'ANs'}[sc])
a.set_yscale('log'); plain_log(a,'y'); a.set_ylim(0.12,12); a.set_xticks(np.arange(3)); a.set_xticklabels(['ACh','GABA','Glu']); a.set_xlabel('predicted transmitter'); a.set_ylabel('neck diameter (µm)')
a.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,handlelength=1.0,labelspacing=0.2,borderaxespad=0.1)
panel_labels(fig,ax,'abc',rows=['abc'])
export(fig,'Figure4',axes=ax,ids=list('abc'))
