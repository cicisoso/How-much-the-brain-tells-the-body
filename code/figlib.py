import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# the optional panel-alignment gate comes from the nature-figure toolkit; a no-op fallback is defined below
try:
    try:
    from audit_panel_alignment import require_matplotlib_panel_alignment  # optional QA gate (nature-figure toolkit)
except ImportError:
    pass
except Exception:
    def require_matplotlib_panel_alignment(*a,**k): return None
# Communications Biology figure lettering: Arial/Helvetica 8-12 pt, minimal size variation, lowercase panel letters
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans','Liberation Sans']
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
plt.rcParams.update({'font.size':8,'axes.labelsize':8,'axes.titlesize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,
    'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6,'xtick.major.size':2.5,'ytick.major.size':2.5,'legend.frameon':False,'lines.linewidth':1.0,'mathtext.default':'regular'})
MM=1/25.4; W=180*MM
# semantic colours: descending = blue family; ascending = orange family; sensory ascending = teal; sensory descending = violet
C={'DN':'#0F4D92','DN_light':'#7FA6D6','AN':'#D9782D','AN_light':'#F0B98A','SA':'#2E8B8B','SD':'#7B4F9E','EFF':'#4D4D4D','UNANN':'#C9C9C9','ECS':'#FFFFFF','male':'#0F4D92','female':'#B8386F','grey':'#767676','light':'#E6E6E6'}
CLASS_LABEL={'descending_neuron':'Descending (DN)','ascending_neuron':'Ascending (AN)','sensory_ascending':'Sensory ascending (SA)','sensory_descending':'Sensory descending (SD)'}
CLASS_COLOR={'descending_neuron':C['DN'],'ascending_neuron':C['AN'],'sensory_ascending':C['SA'],'sensory_descending':C['SD']}
def panel_labels_tight(fig,axes,letters,dx_pt=0,dy_pt=2):
    fig.canvas.draw(); r=fig.canvas.get_renderer()
    for ax,l in zip(axes,letters):
        bb=ax.get_tightbbox(r).transformed(fig.transFigure.inverted())
        fig.text(max(0.004,bb.x0+dx_pt/72/fig.get_figwidth()),min(1-12/72/fig.get_figheight(),bb.y1+dy_pt/72/fig.get_figheight()),l,fontsize=10,fontweight='bold',ha='left',va='bottom')
def export(fig,base):
    fig.savefig(base+'.pdf'); fig.savefig(base+'.svg'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300); print('saved',base)
def plain_log(ax,axis='both'):
    fmt=matplotlib.ticker.FuncFormatter(lambda v,p: f'{v:g}' if v<1e4 else f'{v/1000:g}k')
    if axis in ('x','both'): ax.xaxis.set_major_formatter(fmt)
    if axis in ('y','both'): ax.yaxis.set_major_formatter(fmt)
def fit_line(ax,x,y,color,xlim=None,lw=1.2):
    ok=(x>0)&(y>0); lx=np.log10(x[ok]); ly=np.log10(y[ok]); s,i=np.polyfit(lx,ly,1)
    xx=np.array(xlim if xlim else [x[ok].min(),x[ok].max()]); ax.plot(xx,10**(i+s*np.log10(xx)),color=color,lw=lw,zorder=4); return s
