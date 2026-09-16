from figlib import *
# publication settings (also set in figlib): Arial/Helvetica sans-serif, editable text, alignment gate before export
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans']
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})
try:
    from audit_panel_alignment import require_matplotlib_panel_alignment  # optional QA gate (nature-figure toolkit)
except ImportError:
    pass
import matplotlib.gridspec as gridspec
R='results'; ch=pd.read_parquet(f'{R}/channel_axons_male.parquet'); dl=pd.read_parquet(f'{R}/channel_axons_male_delay.parquet')
fig=plt.figure(figsize=(180*MM,72*MM))
gs=gridspec.GridSpec(1,4,width_ratios=[1,1,1,1],wspace=1.05,left=0.165,right=0.985,top=0.86,bottom=0.33)
axA,axB,axC,axD=[fig.add_subplot(gs[0,i]) for i in range(4)]
def share_bars(ax,df,key,order,labels,title):
    g=df.groupby(key).agg(n=('body','size'),area=('area_um2','sum')).reindex(order).fillna(0)
    cnt=g.n/g.n.sum(); ar=g.area/g.area.sum(); y=np.arange(len(order))[::-1]
    ax.barh(y+0.18,cnt,height=0.34,color=C['light'],edgecolor=C['grey'],lw=0.4,label='axon count'); ax.barh(y-0.18,ar,height=0.34,color={'DN':C['DN'],'AN':C['AN'],'SA':C['SA']}[title[:2]],label='cross-section area')
    ax.set_yticks(y); ax.set_yticklabels(labels); ax.set_xlabel('share'); ax.set_title(title); ax.set_xlim(0,max(cnt.max(),ar.max())*1.15)
dn=ch[ch.superclass=='descending_neuron']
order=['xn','xl','lt','ut','fl','nt','it','ht','ad','other']; dn=dn.assign(sub=dn.subclass.where(dn.subclass.isin(order[:-1]),'other'))
lab=['multiple neuropils','all leg neuropils','lower tectulum','upper tectulum','front leg','neck tectulum','interm. tectulum','haltere tectulum','abdomen','other']
share_bars(axA,dn,'sub',order,lab,'DN by VNC target')
an=ch[ch.superclass=='ascending_neuron']; an=an.assign(seg=an.somaNeuromere.where(an.somaNeuromere.isin(['T1','T2','T3']),'abdomen'))
share_bars(axB,an,'seg',['T1','T2','T3','abdomen'],['T1 (front legs)','T2 (wings, mid legs)','T3 (halteres, hind legs)','abdomen'],'AN by segment of origin')
sa=ch[ch.superclass=='sensory_ascending']; sa=sa.assign(org=sa.subclass.where(sa.subclass.isin(['campaniform sensilla','haltere','leg bristle','chordotonal organ','abdomen']),'other'))
share_bars(axC,sa,'org',['campaniform sensilla','haltere','leg bristle','chordotonal organ','abdomen','other'],['wing campaniform','haltere','leg bristle','chordotonal','abdomen','other'],'SA by sense organ')
h,l=axC.get_legend_handles_labels(); fig.legend(h,l,loc='lower center',bbox_to_anchor=(0.42,0.0),ncol=2,fontsize=8,handlelength=1.2,columnspacing=1.5,frameon=False)
d=dl[(dl.end_brain=='branch')&(dl.end_vnc=='branch')]
bins=np.linspace(0,2.5,40)
for sc in ['descending_neuron','ascending_neuron']:
    x=d[d.superclass==sc].latency_ms; axD.hist(x,bins=bins,histtype='step',color=CLASS_COLOR[sc],lw=1.0,label={'descending_neuron':'DN','ascending_neuron':'AN'}[sc]); axD.axvline(x.median(),color=CLASS_COLOR[sc],lw=0.6,ls=':')
gf=dl[dl.type=='DNp01'].latency_ms.max(); axD.plot([gf,gf],[10,212],color=C['grey'],lw=0.5,clip_on=False); axD.text(gf,218,'GF',ha='center',va='bottom',fontsize=8,color=C['DN'],clip_on=False)
axD.set_xlabel('conduction time (ms)'); axD.set_ylabel('axons'); axD.set_ylim(0,210); axD.legend(fontsize=8,handlelength=1.2,loc='upper right',labelspacing=0.2)
panel_labels_tight(fig,[axA,axB,axC,axD],'abcd')
base='figures/Figure3'; require_matplotlib_panel_alignment(fig,json_out=base+'.alignment.json',overlay_svg=base+'.alignment.svg',tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=False,strict=True)
fig.savefig(base+'.pdf'); fig.savefig(base+'.svg'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300); print('saved',base)
