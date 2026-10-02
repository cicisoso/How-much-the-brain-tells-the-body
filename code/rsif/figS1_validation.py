"""Figure S1. Additional validation of the calibre measurements."""
from figlib_rsif import *
import matplotlib.gridspec as gridspec
from scipy import stats
g=pd.read_parquet(f'{R}/caliber_per_slab.parquet'); ch=pd.read_parquet(f'{R}/channel_axons_male.parquet')
dl=pd.read_parquet(f'{D}/neck_delay_budget.parquet'); mv=pd.read_parquet(f'{R}/male_mesh_validation.parquet').dropna()
fig=plt.figure(figsize=(W,58*MM)); gs=gridspec.GridSpec(1,4,wspace=0.62,left=0.07,right=0.975,top=0.86,bottom=0.25)
axA,axB,axC,axD=[fig.add_subplot(gs[0,i]) for i in range(4)]; lim=[0.08,15]
gg=g[g.body.isin(ch.body)&(g.n_planes>=32)].pivot(index='body',columns='slab',values='diam_um').dropna()
axA.scatter(gg['A'],gg['C'],s=3,color=C['grey'],alpha=0.5,lw=0,rasterized=True); axA.plot(lim,lim,'k--',lw=0.5); axA.set_xscale('log'); axA.set_yscale('log'); plain_log(axA); axA.set_xlim(lim); axA.set_ylim(lim)
axA.set_xlabel('diameter, anterior slab (µm)'); axA.set_ylabel('diameter, posterior slab (µm)'); axA.text(0.04,0.96,f'n = {len(gg)}\nρ = {stats.spearmanr(gg.A,gg.C)[0]:.2f}',transform=axA.transAxes,va='top')
ct=g[g.body.isin(ch.body)].cos_theta.dropna(); axB.hist(ct,bins=np.linspace(0.5,1,51),color=C['DN_light'],lw=0); axB.set_xlabel('cos θ (axon versus slab normal)'); axB.set_ylabel('axon–slab pairs')
axB.text(0.04,0.96,f'median {ct.median():.3f}',transform=axB.transAxes,va='top')
lim2=[0.01,120]; axC.scatter(mv.area_um2,mv.area_mesh_um2,s=5,color=C['grey'],alpha=0.7,lw=0,rasterized=True); axC.plot(lim2,lim2,'k--',lw=0.5); axC.set_xscale('log'); axC.set_yscale('log'); plain_log(axC); axC.set_xlim(lim2); axC.set_ylim(lim2)
axC.set_xlabel('voxel area (µm²)'); axC.set_ylabel('mesh-section area (µm²)'); axC.text(0.04,0.96,f'n = {len(mv)}\nmedian ratio {np.median(mv.area_mesh_um2/mv.area_um2):.2f}',transform=axC.transAxes,va='top')
mm=ch.merge(dl[['bodyId','r_skel_neck']],on='bodyId'); axD.scatter(mm.r_skel_neck*2*0.008,mm.diam_um,s=3,color=C['grey'],alpha=0.4,lw=0,rasterized=True); axD.plot(lim,lim,'k--',lw=0.5)
axD.set_xscale('log'); axD.set_yscale('log'); plain_log(axD); axD.set_xlim(lim); axD.set_ylim(lim); axD.set_xlabel('skeleton diameter (µm)'); axD.set_ylabel('segmentation diameter (µm)')
axD.text(0.04,0.96,f'n = {len(mm)}\nρ = {stats.spearmanr(mm.r_skel_neck,mm.diam_um)[0]:.2f}',transform=axD.transAxes,va='top')
panel_labels(fig,[axA,axB,axC,axD],'abcd',rows=['abcd'])
export(fig,'FigureS1',axes=[axA,axB,axC,axD],ids=list('abcd'))
