import os
DATA_DIR=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
from figlib import *
import matplotlib.gridspec as gridspec
from scipy import stats
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans','Liberation Sans']
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})
from audit_panel_alignment import require_matplotlib_panel_alignment
D=DATA_DIR; R='results'
cal=pd.read_parquet(f'{R}/caliber_per_body.parquet'); g=pd.read_parquet(f'{R}/caliber_per_slab.parquet'); ch=pd.read_parquet(f'{R}/channel_axons_male.parquet')
dl=pd.read_parquet(f'{D}/neck_delay_budget.parquet'); mv=pd.read_parquet(f'{R}/male_mesh_validation.parquet').dropna()
fig=plt.figure(figsize=(174*MM,55*MM)); gs=gridspec.GridSpec(1,4,wspace=0.6,left=0.06,right=0.985,top=0.86,bottom=0.24)
axA,axB,axC,axD=[fig.add_subplot(gs[0,i]) for i in range(4)]
# A: consistency along the connective: slab A vs slab C diameters
gg=g[g.body.isin(ch.body)&(g.n_planes>=32)].pivot(index='body',columns='slab',values='diam_um').dropna()
axA.scatter(gg['A'],gg['C'],s=3,color=C['grey'],alpha=0.5,lw=0,rasterized=True); lim=[0.08,15]; axA.plot(lim,lim,'k--',lw=0.5); axA.set_xscale('log'); axA.set_yscale('log'); plain_log(axA); axA.set_xlim(lim); axA.set_ylim(lim)
axA.set_xlabel('diameter, anterior slab (µm)'); axA.set_ylabel('diameter, posterior slab (µm)'); axA.text(0.04,0.96,f'n = {len(gg)}\nρ = {stats.spearmanr(gg.A,gg.C)[0]:.2f}',transform=axA.transAxes,va='top',fontsize=6)
# B: tilt correction factor
ct=g[g.body.isin(ch.body)].cos_theta.dropna(); axB.hist(ct,bins=np.linspace(0.5,1,51),color=C['DN_light'],lw=0); axB.set_xlabel('cos θ (axon vs. slab normal)'); axB.set_ylabel('axon–slab pairs'); axB.text(0.04,0.96,f'median {ct.median():.3f}',transform=axB.transAxes,va='top',fontsize=6)
# C: mesh section vs voxel area
axC.scatter(mv.area_um2,mv.area_mesh_um2,s=5,color=C['grey'],alpha=0.7,lw=0,rasterized=True); lim2=[0.01,120]; axC.plot(lim2,lim2,'k--',lw=0.5); axC.set_xscale('log'); axC.set_yscale('log'); plain_log(axC); axC.set_xlim(lim2); axC.set_ylim(lim2)
axC.set_xlabel('voxel area (µm²)'); axC.set_ylabel('mesh-section area (µm²)'); axC.text(0.04,0.96,f'n = {len(mv)}\nmedian ratio {np.median(mv.area_mesh_um2/mv.area_um2):.2f}',transform=axC.transAxes,va='top',fontsize=6)
# D: skeleton radius vs EM diameter
m=ch.merge(dl[['bodyId','r_skel_neck']],on='bodyId'); axD.scatter(m.r_skel_neck*2*0.008,m.diam_um,s=3,color=C['grey'],alpha=0.4,lw=0,rasterized=True); axD.plot(lim,lim,'k--',lw=0.5); axD.set_xscale('log'); axD.set_yscale('log'); plain_log(axD); axD.set_xlim(lim); axD.set_ylim(lim)
axD.set_xlabel('skeleton-radius diameter (µm)'); axD.set_ylabel('EM diameter (µm)'); axD.text(0.04,0.96,f'n = {len(m)}\nρ = {stats.spearmanr(m.r_skel_neck,m.diam_um)[0]:.2f}',transform=axD.transAxes,va='top',fontsize=6)
panel_labels_tight(fig,[axA,axB,axC,axD],'ABCD')
base='figures/FigureS1'; require_matplotlib_panel_alignment(fig,json_out=base+'.alignment.json',overlay_svg=base+'.alignment.svg',tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=False,strict=True)
fig.savefig(base+'.pdf'); fig.savefig(base+'.svg'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300); print('saved',base)
