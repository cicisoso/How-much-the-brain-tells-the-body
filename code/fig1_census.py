import os
DATA_DIR=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
from figlib import *
# publication settings (also set in figlib): Arial/Helvetica sans-serif, editable text, alignment gate before export
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans']
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})
try:
    from audit_panel_alignment import require_matplotlib_panel_alignment  # optional QA gate (nature-figure toolkit)
except ImportError:
    pass
import matplotlib.gridspec as gridspec
from scipy import stats
D=DATA_DIR; R='results'
ann=pd.read_feather(f'{D}/body-annotations-male-cns-v1.0-minconf-0.5.feather')
cal=pd.read_parquet(f'{R}/caliber_per_body.parquet'); ch=pd.read_parquet(f'{R}/channel_axons_male.parquet')
neck=pd.read_csv(f'{D}/neck_neurons.csv'); rng=np.random.default_rng(1)
pts=[]
for b in rng.choice(neck.bodyId.values,900,replace=False):
    d=pd.read_csv(f'{D}/swc/{b}.swc',sep=r'\s+',comment='#',header=None,usecols=[2,3,4],names=['x','y','z']); pts.append(d.values[::3])
P=np.vstack(pts)*0.008
z=np.load(f'{D}/neck_slabs/plane_B_z26976.npz'); pl=z['plane']; H=np.load(f'{R}/plane_B_hullmask.npz')['hull']
ids=np.unique(pl); cat=np.full(len(ids),5,dtype=np.uint8)
cls=dict(zip(cal[cal.neck_listed].body.astype(np.int64),cal[cal.neck_listed].superclass)); sup=dict(zip(ann.bodyId.astype(np.int64),ann.superclass.astype(str)))
code={'descending_neuron':0,'ascending_neuron':1,'sensory_ascending':2,'sensory_descending':3}
for i,b in enumerate(ids):
    if b==0: cat[i]=6
    elif b in cls and cls[b] in code: cat[i]=code[cls[b]]
    elif b in sup and ('motor' in sup[b] or 'efferent' in sup[b]): cat[i]=4
img=cat[np.searchsorted(ids,pl)]
pal=np.array([matplotlib.colors.to_rgb(c) for c in [C['DN'],C['AN'],C['SA'],C['SD'],C['EFF'],C['UNANN'],C['ECS']]])
rgb=pal[img]; rgb[~H]=np.array([1,1,1])*0.97
ys,xs=np.where(H); x0,x1,y0,y1=xs.min()-40,xs.max()+40,ys.min()-40,ys.max()+40
crop=rgb[y0:y1,x0:x1]; ext=[0,(x1-x0)*0.016,(y1-y0)*0.016,0]
lst=cal[cal.neck_listed&cal.type.notna()&cal.somaSide.isin(['L','R'])]
t=lst.groupby(['type','somaSide']).diam_um.median().unstack().dropna(); rho=stats.spearmanr(t.L,t.R)[0]
mv=pd.read_parquet(f'{R}/mip_validation.parquet')
fig=plt.figure(figsize=(180*MM,104*MM))
gs=gridspec.GridSpec(2,3,width_ratios=[1.15,1.9,1.2],height_ratios=[1,1],wspace=0.55,hspace=0.6,left=0.075,right=0.99,top=0.9,bottom=0.145)
axA=fig.add_subplot(gs[:,0]); axB=fig.add_subplot(gs[:,1]); axC=fig.add_subplot(gs[0,2]); axD=fig.add_subplot(gs[1,2])
axA.hexbin(P[:,0],P[:,2],gridsize=(60,140),bins='log',cmap='Greys',mincnt=1,linewidths=0.1,rasterized=True)
zn=54000*0.008; axA.axhline(zn,color=C['AN'],lw=0.8,ls='--'); axA.text(P[:,0].max()+8,zn,'neck\nplane',color=C['AN'],fontsize=8,va='center')
axA.invert_yaxis(); axA.set_aspect('equal'); axA.set_xlabel('x (µm)'); axA.set_ylabel('anterior–posterior (µm)')
axA.text(P[:,0].max()+8,P[:,2].min()+120,'brain',ha='left',fontsize=8,color=C['grey']); axA.text(P[:,0].max()+8,P[:,2].max()-200,'VNC',ha='left',fontsize=8,color=C['grey'])
axA.set_title('3,721 neck-crossing\nneurons')
axB.imshow(crop,extent=ext,interpolation='nearest',rasterized=True); axB.set_xticks([]); axB.set_yticks([])
for s in axB.spines.values(): s.set_visible(False)
axB.plot([ext[1]-14,ext[1]-4],[ext[2]-3,ext[2]-3],color='k',lw=1.5); axB.text(ext[1]-9,ext[2]-5,'10 µm',ha='center',va='bottom',fontsize=8)
from matplotlib.patches import Patch
handles=[Patch(color=C['DN'],label='Descending (DN)'),Patch(color=C['AN'],label='Ascending (AN)'),Patch(color=C['SA'],label='Sensory ascending (SA)'),Patch(color=C['SD'],label='Sensory descending (SD)'),Patch(color=C['EFF'],label='Efferent (exits via nerve)'),Patch(color=C['UNANN'],label='Glia / unannotated')]
axB.legend(handles=handles,loc='upper center',bbox_to_anchor=(0.5,-0.01),ncol=3,fontsize=8,handlelength=1,columnspacing=1.0,labelspacing=0.3)
axB.set_title('Neck connective cross-section (16 nm EM segmentation)')
axC.scatter(t.L,t.R,s=4,color=C['grey'],alpha=0.6,lw=0,rasterized=True); lim=[0.08,15]; axC.plot(lim,lim,color='k',lw=0.5,ls='--')
axC.set_xscale('log'); axC.set_yscale('log'); plain_log(axC); axC.set_xlim(lim); axC.set_ylim(lim); axC.set_xlabel('diameter, left (µm)'); axC.set_ylabel('diameter, right (µm)')
axC.text(0.04,0.96,f'{len(t)} cell types\nρ = {rho:.2f}',transform=axC.transAxes,va='top',fontsize=8)
axD.scatter(mv.area_um2_raw,mv.area_um2_mip0,s=4,color=C['grey'],alpha=0.6,lw=0,rasterized=True); lim2=[0.008,100]; axD.plot(lim2,lim2,color='k',lw=0.5,ls='--')
axD.set_xscale('log'); axD.set_yscale('log'); plain_log(axD); axD.set_xlim(lim2); axD.set_ylim(lim2); axD.set_xlabel('area at 16 nm (µm²)'); axD.set_ylabel('area at 8 nm (µm²)')
axD.text(0.04,0.96,f'{len(mv)} axons\nmedian ratio {mv.ratio.median():.2f}',transform=axD.transAxes,va='top',fontsize=8)
panel_labels_tight(fig,[axA,axB,axC,axD],'abcd')
base='figures/Figure1'
require_matplotlib_panel_alignment(fig,json_out=base+'.alignment.json',overlay_svg=base+'.alignment.svg',tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=False,strict=True,axes=[axA,axB,axC,axD],panel_ids=['a','b','c','d'],column_groups=[['c','d']],
     exemptions=[{'panels':['a','b'],'checks':['row','column','panel-width','horizontal-gutter','vertical-gutter'],'reason':'a and b are full-height image/schematic panels with intentional unequal widths'}])
fig.savefig(base+'.pdf'); fig.savefig(base+'.svg'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300); print('saved',base)
