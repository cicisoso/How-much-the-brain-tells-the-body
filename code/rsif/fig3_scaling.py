"""Figure 3. Synaptic output scales as d^1.4, and the exponent splits into cable length and release-site density."""
from figlib_rsif import *
import matplotlib.gridspec as gridspec
m=pd.read_parquet(f'{R}/channel_axons_male_arbor.parquet'); J=json.load(open(f'{R}/connectivity_scaling.json')); A=json.load(open(f'{R}/arbor_scaling.json'))
m=m[m.diam_um>=0.1]
dn=m[m.superclass=='descending_neuron']; an=m[m.superclass=='ascending_neuron']; sa=m[m.superclass=='sensory_ascending']
fig=plt.figure(figsize=(W,118*MM))
gs=gridspec.GridSpec(2,3,wspace=0.55,hspace=0.62,left=0.1,right=0.985,top=0.9,bottom=0.1)
ax=[fig.add_subplot(gs[i,j]) for i in range(2) for j in range(3)]
xlim=[0.09,12]
def scat(a,d,col,color,label,key):
    a.scatter(d.diam_um,d[col].clip(lower=1),s=3,color=C['grey'],alpha=0.45,lw=0,rasterized=True)
    fit_line(a,d.diam_um.values,d[col].values.astype(float),color,xlim=[d.diam_um.min(),d.diam_um.max()])
    j=J[key]; a.text(0.97,0.03,f'{label}, n = {num(j["n"])}\nslope {j["slope"]:.2f} {ci(*j["slope_ci"])}\nρ = {j["rho"]:.2f}'.replace('-','−'),transform=a.transAxes,va='bottom',ha='right')
    a.set_xscale('log'); a.set_yscale('log'); plain_log(a); a.set_xlim(xlim); a.set_ylim(8,3e5); a.set_xlabel('neck diameter (µm)')
scat(ax[0],dn,'n_out_vnc',C['DN'],'DN','descending_neuron:n_out_vnc'); ax[0].set_ylabel('output synapses in the VNC')
scat(ax[1],an,'n_out_brain',C['AN'],'AN','ascending_neuron:n_out_brain'); ax[1].set_ylabel('output synapses in the brain')
scat(ax[2],sa,'n_out_brain',C['SA'],'SA','sensory_ascending:n_out_brain'); ax[2].set_ylabel('output synapses in the brain')
# (d) cable built in the target region, (e) release sites per micrometre of that cable
dd=dn[(dn.n_out_vnc>0)&(dn.cable_vnc_um>0)&(dn.tbar_vnc>0)]; aa=an[(an.n_out_brain>0)&(an.cable_brain_um>0)&(an.tbar_brain>0)]
for a,ylab,fD,fA,key in [(ax[3],'cable in target region (µm)',lambda d: d.cable_vnc_um,lambda d: d.cable_brain_um,'L'),
                         (ax[4],'release sites per µm of cable',lambda d: d.tbar_vnc/d.cable_vnc_um,lambda d: d.tbar_brain/d.cable_brain_um,'T_per_L')]:
    for d,f,color,lab,sc in [(dd,fD,C['DN'],'DN','descending_neuron'),(aa,fA,C['AN'],'AN','ascending_neuron')]:
        y=f(d).values; a.scatter(d.diam_um,y,s=2.5,color=color,alpha=0.22,lw=0,rasterized=True); fit_line(a,d.diam_um.values,y,color,xlim=[d.diam_um.min(),d.diam_um.max()])
        s=A[sc]['decomposition'][key]['slope']; a.plot([],[],color=color,lw=1.2,label=f'{lab}: slope {s["coef"]:.2f} {ci(*s["ci"])}')
    a.set_xscale('log'); a.set_yscale('log'); plain_log(a); a.set_xlim(xlim); a.set_xlabel('neck diameter (µm)'); a.set_ylabel(ylab)
    a.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),handlelength=1.2,labelspacing=0.2,borderaxespad=0.1)
ax[3].set_ylim(5,2e4); ax[4].set_ylim(0.02,5)
# (f) the exponent of synaptic output is the sum of three exponents
a=ax[5]; comps=[('L','cable','#3B3B3B'),('T_per_L','sites/µm','#8C8C8C'),('S_per_T','contacts/site','#D0D0D0')]
cls=[('descending_neuron','DN'),('ascending_neuron','AN'),('sensory_ascending','SA')]; w=0.2
for i,(sc,lab) in enumerate(cls):
    dec=A[sc]['decomposition']
    for k,(key,name,col) in enumerate(comps):
        s=dec[key]['slope']; x0=i+(k-1.5)*w
        a.bar(x0,s['coef'],width=w,color=col,edgecolor='k',lw=0.3,label=name if i==0 else None)
        a.errorbar(x0,s['coef'],yerr=[[s['coef']-s['ci'][0]],[s['ci'][1]-s['coef']]],color='k',lw=0.6,capsize=1.2)
    jj=J[{'descending_neuron':'descending_neuron:n_out_vnc','ascending_neuron':'ascending_neuron:n_out_brain','sensory_ascending':'sensory_ascending:n_out_brain'}[sc]]
    s={'coef':jj['slope'],'ci':jj['slope_ci']}; x0=i+1.5*w   # output exponent as in panels (a)-(c) and the text
    a.errorbar(x0,s['coef'],yerr=[[s['coef']-s['ci'][0]],[s['ci'][1]-s['coef']]],fmt='D',ms=3.5,color=CLASS_COLOR[sc],mec='k',mew=0.4,lw=0.8,capsize=1.5,label='output' if i==0 else None)
a.axhline(0,color='k',lw=0.5); a.set_xticks(range(3)); a.set_xticklabels([l for _,l in cls]); a.set_xlim(-0.55,2.55); a.set_ylim(-0.75,1.65)
a.set_ylabel('exponent of neck diameter')
a.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,handlelength=1.0,columnspacing=0.8,labelspacing=0.2,borderaxespad=0.1)
panel_labels(fig,ax,'abcdef',rows=['abc','def'])
export(fig,'Figure3',axes=ax,ids=list('abcdef'))
