"""Figure 2. Descending axons are fewer but hold most of the cross-section; allocation by target, origin and sense organ."""
from figlib_rsif import *
import matplotlib.gridspec as gridspec
ch=pd.read_parquet(f'{R}/channel_axons_male.parquet'); B=json.load(open(f'{R}/budget_male.json'))
fig=plt.figure(figsize=(W,128*MM))
gs=gridspec.GridSpec(2,3,width_ratios=[1.25,1,1],height_ratios=[1,1.05],wspace=0.62,hspace=0.75,left=0.165,right=0.985,top=0.9,bottom=0.15)
axA,axB,axC=[fig.add_subplot(gs[0,i]) for i in range(3)]; axD,axE,axF=[fig.add_subplot(gs[1,i]) for i in range(3)]
# (a) diameter distributions
bins=np.logspace(np.log10(0.08),np.log10(12),45)
for sc in ['ascending_neuron','sensory_ascending','descending_neuron']:
    d=ch[ch.superclass==sc].diam_um; short={'ascending_neuron':'AN','sensory_ascending':'SA','descending_neuron':'DN'}[sc]
    axA.hist(d,bins=bins,histtype='step',lw=1.0,color=CLASS_COLOR[sc],label=f'{short} (n = {num(len(d))})'); axA.axvline(d.median(),ymax=0.62,color=CLASS_COLOR[sc],lw=0.6,ls=':')
axA.set_xscale('log'); plain_log(axA,'x'); axA.set_xlabel('axon diameter at the neck (µm)'); axA.set_ylabel('axons'); axA.set_ylim(0,275)
axA.legend(loc='upper left',handlelength=1.0,borderaxespad=0.1,labelspacing=0.25)
top=ch.sort_values('diam_um',ascending=False).drop_duplicates('type').head(2)
for _,r in top.iterrows(): axA.annotate(r.type.replace('DNp01','DNp01 (GF)'),xy=(r.diam_um,2),xytext=(r.diam_um,95 if r.type=='DNp01' else 50),ha='center',rotation=90,va='bottom',arrowprops=dict(arrowstyle='-',lw=0.4,color=C['grey']))
# (b) Lorenz curve of area
a=np.sort(ch.area_um2.values)[::-1]; cs=np.cumsum(a)/a.sum(); x=np.arange(1,len(a)+1)/len(a)
axB.plot(x,cs,color=C['DN'],lw=1.2); axB.plot([0,1],[0,1],color=C['grey'],lw=0.5,ls='--')
n50=B['n_axons_for_half_area']; axB.plot([n50/len(a)],[0.5],'o',ms=3,color=C['AN'])
axB.annotate(f'{n50} axons hold\nhalf of the area',xy=(n50/len(a),0.5),xytext=(0.34,0.18),arrowprops=dict(arrowstyle='-',lw=0.4,color=C['grey']))
axB.set_xlabel('fraction of axons (thickest first)'); axB.set_ylabel('cumulative area fraction'); axB.set_xlim(0,1); axB.set_ylim(0,1); axB.text(0.6,0.05,f'Gini = {B["gini_area"]:.2f}')
# (c) shares by count, summed diameter and summed area
dn=ch[ch.dir=='down']; up=ch[ch.dir=='up']
vd=np.array([len(dn),dn.diam_um.sum(),dn.area_um2.sum()]); vu=np.array([len(up),up.diam_um.sum(),up.area_um2.sum()]); fd=vd/(vd+vu)
xx=np.arange(3); axC.bar(xx,fd,color=C['DN'],width=0.6,label='brain → body'); axC.bar(xx,1-fd,bottom=fd,color=C['AN'],width=0.6,label='body → brain')
for i,f in enumerate(fd): axC.text(i,f/2,f'{100*f:.0f}%',ha='center',va='center',color='white'); axC.text(i,f+(1-f)/2,f'{100*(1-f):.0f}%',ha='center',va='center',color='white')
axC.set_xticks(xx); axC.set_xticklabels(['count','Σ$\mathit{d}$','Σ$\mathit{d}$²']); axC.set_ylabel('share of the channel'); axC.set_ylim(0,1)
axC.legend(loc='upper center',bbox_to_anchor=(0.5,-0.12),ncol=1,handlelength=1,borderaxespad=0.1,labelspacing=0.2)
# (d)-(f) count and area shares by target, origin and sense organ
def share_bars(ax,df,key,order,labels,title,color):
    g=df.groupby(key).agg(n=('body','size'),area=('area_um2','sum')).reindex(order).fillna(0)
    cnt=g.n/g.n.sum(); ar=g.area/g.area.sum(); y=np.arange(len(order))[::-1]
    ax.barh(y+0.18,cnt,height=0.34,color=C['light'],edgecolor=C['grey'],lw=0.4,label='axon count'); ax.barh(y-0.18,ar,height=0.34,color=color,label='cross-sectional area')
    ax.set_yticks(y); ax.set_yticklabels(labels); ax.set_xlabel('share'); ax.set_title(title); ax.set_xlim(0,max(cnt.max(),ar.max())*1.15)
dd=ch[ch.superclass=='descending_neuron']; order=['xn','xl','lt','ut','fl','nt','it','ht','ad','other']; dd=dd.assign(sub=dd.subclass.where(dd.subclass.isin(order[:-1]),'other'))
share_bars(axD,dd,'sub',order,['multiple neuropils','all leg neuropils','lower tectulum','upper tectulum','front leg','neck tectulum','intermediate tectulum','haltere tectulum','abdomen','other'],'DNs by VNC target',C['DN'])
an=ch[ch.superclass=='ascending_neuron']; an=an.assign(seg=an.somaNeuromere.where(an.somaNeuromere.isin(['T1','T2','T3']),'abdomen'))
share_bars(axE,an,'seg',['T1','T2','T3','abdomen'],['T1 (front legs)','T2 (wings, mid legs)','T3 (halteres, hind legs)','abdomen'],'ANs by segment of origin',C['AN'])
sa=ch[ch.superclass=='sensory_ascending']; sa=sa.assign(org=sa.subclass.where(sa.subclass.isin(['campaniform sensilla','haltere','leg bristle','chordotonal organ','abdomen']),'other'))
axF.set_xticks([0,0.2,0.4,0.6]); share_bars(axF,sa,'org',['campaniform sensilla','haltere','leg bristle','chordotonal organ','abdomen','other'],['wing campaniform','haltere','leg bristle','chordotonal','abdomen','other'],'SAs by sense organ',C['SA'])
from matplotlib.patches import Patch
fig.legend(handles=[Patch(facecolor=C['light'],edgecolor=C['grey'],lw=0.4,label='axon count'),Patch(color=C['DN'],label='area (DN)'),Patch(color=C['AN'],label='area (AN)'),Patch(color=C['SA'],label='area (SA)')],
           loc='lower center',bbox_to_anchor=(0.575,0.0),ncol=4,handlelength=1.2,columnspacing=1.5)
panel_labels(fig,[axA,axB,axC,axD,axE,axF],'abcdef',rows=['abc','def'])
export(fig,'Figure2',axes=[axA,axB,axC,axD,axE,axF],ids=list('abcdef'),
       exemptions=[{'panels':['a'],'checks':['panel-width'],'reason':'distribution panel intentionally wider'},
                   {'panels':['d','e','f'],'checks':['column','panel-width','horizontal-gutter'],'reason':'horizontal bar panels carry category labels of different lengths left of the axes'}])
