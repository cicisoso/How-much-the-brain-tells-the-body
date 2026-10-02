"""Figure 5. The organization of the channel is shared between a male and a female connectome."""
from figlib_rsif import *
import matplotlib.gridspec as gridspec
m=pd.read_parquet(f'{R}/channel_axons_male.parquet'); f=pd.read_parquet(f'{R}/channel_axons_female.parquet'); j=pd.read_parquet(f'{R}/type_matched_relative.parquet'); J=json.load(open(f'{R}/compare_sex2.json'))
fig=plt.figure(figsize=(W,60*MM))
gs=gridspec.GridSpec(1,4,width_ratios=[1,1,1.35,1],wspace=0.72,left=0.075,right=0.985,top=0.8,bottom=0.24)
axA,axB,axC,axD=[fig.add_subplot(gs[0,i]) for i in range(4)]
axA.scatter(j.male_abs,j.female_abs,s=4,color=C['grey'],alpha=0.6,lw=0,rasterized=True); lim=[0.08,15]; axA.plot(lim,lim,'k--',lw=0.5)
k=2**J['matched_abs_median_log2_f_over_m']; xs=np.array([lim[0]/k*1.02,lim[1]*0.98]); axA.plot(xs,k*xs,color=C['female'],lw=0.8)
axA.set_xscale('log'); axA.set_yscale('log'); plain_log(axA); axA.set_xlim(lim); axA.set_ylim(lim); axA.set_xlabel('male diameter (µm)'); axA.set_ylabel('female diameter (µm)')
axA.text(0.98,0.04,f'n = {len(j)}\nρ = {J["matched_abs_spearman"]:.2f}\nslope {k:.2f}',transform=axA.transAxes,ha='right',va='bottom')
xx=np.arange(3); sm=[J['male_share_count'],J['male_share_sumd'],J['male_share_area']]; sf=[J['female_share_count'],J['female_share_sumd'],J['female_share_area']]
axB.bar(xx-0.18,sm,width=0.34,color=C['male'],label='male (MaleCNS)'); axB.bar(xx+0.18,sf,width=0.34,color=C['female'],label='female (BANC)'); axB.axhline(0.5,color=C['grey'],lw=0.5,ls='--')
axB.set_xticks(xx); axB.set_xticklabels(['count','Σ$\\mathit{d}$','Σ$\\mathit{d}$²']); axB.set_ylabel('descending share of channel'); axB.set_ylim(0,0.8)
axB.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=1,handlelength=1,borderaxespad=0.1)
dm=m[m.superclass=='descending_neuron']; df_=f[f.superclass=='descending_neuron']
data=[dm[dm.dimcat=='isomorphic'].rel,dm[dm.dimcat=='male-specific'].rel,df_[df_.dimcat=='isomorphic'].rel,df_[df_.dimcat=='female-specific'].rel]
pos=[0.5,1.5,2.5,3.5]; cols=[C['light'],C['male'],C['light'],C['female']]
bp=axC.boxplot(data,positions=pos,widths=0.3,whis=(5,95),showfliers=False,patch_artist=True,medianprops=dict(color='k',lw=0.8),whiskerprops=dict(lw=0.5),capprops=dict(lw=0.5))
for patch,c in zip(bp['boxes'],cols): patch.set_facecolor(c); patch.set_edgecolor('k'); patch.set_linewidth(0.5)
rng=np.random.default_rng(0)
for x0,d,c in zip([1.5,3.5],[data[1],data[3]],[C['male'],C['female']]): axC.scatter(x0+rng.uniform(-0.08,0.08,len(d)),d,s=6,color=c,edgecolor='k',lw=0.3,zorder=3)
axC.set_yscale('log'); plain_log(axC,'y'); axC.set_ylim(0.12,14); axC.set_xticks(pos); axC.set_xticklabels(['shared','sex-\nspecific','shared','sex-\nspecific']); axC.set_xlim(0,4)
axC.set_ylabel('relative calibre of DNs'); axC.axhline(1,color=C['grey'],lw=0.5,ls='--')
axC.text(1.0,10.5,'male',ha='center',color=C['male']); axC.text(3.0,10.5,'female',ha='center',color=C['female'])
d=j[j.dim=='sexually dimorphic']
axD.scatter(d.male,d.female,s=6,color=C['grey'],alpha=0.8,lw=0,rasterized=True); lim2=[0.15,6]; axD.plot(lim2,lim2,'k--',lw=0.5); axD.set_xscale('log'); axD.set_yscale('log'); plain_log(axD); axD.set_xlim(lim2); axD.set_ylim(lim2)
axD.set_xlabel('relative calibre, male'); axD.set_ylabel('relative calibre, female')
slots={'DNb05':4.8,'DNg37':3.5,'DNde002':2.55,'aSP22':1.85}
for t,r in d.sort_values('log2_rel_m_over_f',ascending=False).head(4).iterrows():
    axD.annotate(t,xy=(r.male,r.female),xytext=(0.19,slots.get(t,2.0)),ha='left',va='center',arrowprops=dict(arrowstyle='-',lw=0.3,color=C['grey'],shrinkA=9,shrinkB=2,relpos=(1,0.5)))
axD.text(0.98,0.04,f'n = {len(d)}',transform=axD.transAxes,ha='right',va='bottom')
panel_labels(fig,[axA,axB,axC,axD],'abcd',rows=['abcd'])
export(fig,'Figure5',axes=[axA,axB,axC,axD],ids=list('abcd'),exemptions=[{'panels':['c'],'checks':['panel-width'],'reason':'box-plot panel intentionally wider to separate category labels'}])
