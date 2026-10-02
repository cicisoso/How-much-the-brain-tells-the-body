"""Shared figure style for the J. R. Soc. Interface submission.

Royal Society figure guidance (royalsociety.org/journals/authors/author-guidelines, read 2 October 2026): Times New Roman,
9 pt preferred and never below 7.5 pt; a thin space between values and units and as the thousands separator (no commas);
a full space either side of '='; panel labels italic letters in roman brackets, '(a)', placed consistently at the top
left; at least 300 dpi (TIFF/PNG/EPS/JPG accepted). Figures are drawn at the double-column width of 174 mm."""
import os, sys, json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
if os.environ.get('NATURE_FIGURE_SCRIPTS'): sys.path.insert(0,os.environ['NATURE_FIGURE_SCRIPTS'])   # optional QA toolkit folder
try:
    from audit_panel_alignment import require_matplotlib_panel_alignment
except Exception:                                   # the alignment gate is optional outside the authors' machine
    def require_matplotlib_panel_alignment(*a,**k): return None

import tempfile
# run from the repository root: results/ and figures/rsif/ are relative to the working directory, public data from FLYBRAIN_DATA
R='results'; D=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data')); FIG=os.path.join('figures','rsif'); QA=os.path.join(tempfile.gettempdir(),'rsif_figure_qa')
os.makedirs(FIG,exist_ok=True); os.makedirs(QA,exist_ok=True)

FS=8                                                # 8 pt everywhere (>= 7.5 pt floor); panel labels 9 pt
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Times','DejaVu Serif'],
    'mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold','mathtext.default':'regular',
    'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42,
    'font.size':FS,'axes.labelsize':FS,'axes.titlesize':FS,'xtick.labelsize':FS,'ytick.labelsize':FS,'legend.fontsize':FS,
    'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6,
    'xtick.minor.width':0.4,'ytick.minor.width':0.4,'xtick.major.size':2.5,'ytick.major.size':2.5,'legend.frameon':False,'lines.linewidth':1.0})
MM=1/25.4; W=174*MM
THIN=' '
C={'DN':'#0F4D92','DN_light':'#7FA6D6','AN':'#D9782D','AN_light':'#F0B98A','SA':'#2E8B8B','SD':'#7B4F9E','EFF':'#4D4D4D','UNANN':'#C9C9C9','ECS':'#FFFFFF',
   'male':'#0F4D92','female':'#B8386F','grey':'#767676','light':'#E6E6E6'}
CLASS_COLOR={'descending_neuron':C['DN'],'ascending_neuron':C['AN'],'sensory_ascending':C['SA'],'sensory_descending':C['SD']}

def num(v,dec=None):
    """Royal Society number style: four-digit numbers closed up, thin-space thousands separators from five digits."""
    s=f'{v:.{dec}f}' if dec is not None else (f'{int(round(v))}' if float(v).is_integer() else f'{v:g}')
    ip,_,fp=s.partition('.'); neg=ip.startswith('-'); ip=ip.lstrip('-')
    if len(ip)>4: ip=f'{int(ip):,}'.replace(',',THIN)
    return ('−' if neg else '')+ip+('.'+fp if fp else '')

def ci(lo,hi,dec=2): return f'[{lo:.{dec}f}, {hi:.{dec}f}]'.replace('-','−')

def plain_log(ax,axis='both'):
    fmt=matplotlib.ticker.FuncFormatter(lambda v,p: num(v) if v>=1 else f'{v:g}')
    if axis in ('x','both'): ax.xaxis.set_major_formatter(fmt)
    if axis in ('y','both'): ax.yaxis.set_major_formatter(fmt)

def fit_line(ax,x,y,color,xlim=None,lw=1.2):
    ok=(x>0)&(y>0); lx=np.log10(x[ok]); ly=np.log10(y[ok]); s,i=np.polyfit(lx,ly,1)
    xx=np.array(xlim if xlim else [x[ok].min(),x[ok].max()]); ax.plot(xx,10**(i+s*np.log10(xx)),color=color,lw=lw,zorder=4); return s

def panel_labels(fig,axes,letters,dx_pt=0,dy_pt=2,rows=None):
    """Italic letters in roman brackets, '(a)', at the top left of each panel's tight bounding box. Panels listed together
    in `rows` (strings of letters, e.g. ['abc','def']) share one label height so labels sit in the same position."""
    fig.canvas.draw(); r=fig.canvas.get_renderer()
    bbs={l:ax.get_tightbbox(r).transformed(fig.transFigure.inverted()) for ax,l in zip(axes,letters)}
    top={l:bb.y1 for l,bb in bbs.items()}
    for grp in (rows or []):
        y=max(bbs[l].y1 for l in grp)
        for l in grp: top[l]=y
    for l,bb in bbs.items():
        fig.text(max(0.004,bb.x0+dx_pt/72/fig.get_figwidth()),min(1-11/72/fig.get_figheight(),top[l]+dy_pt/72/fig.get_figheight()),
                 r'($\mathit{%s}$)'%l,fontsize=9,ha='left',va='bottom')

def export(fig,name,axes=None,ids=None,**gate):
    base=os.path.join(FIG,name)
    if axes is not None:
        require_matplotlib_panel_alignment(fig,json_out=os.path.join(QA,name+'.alignment.json'),overlay_svg=os.path.join(QA,name+'.alignment.svg'),
            tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=False,strict=True,axes=axes,panel_ids=ids,**gate)
    fig.savefig(base+'.pdf'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300)
    print('saved',base)

def boxes(ax,data,pos,color,width=0.32,alpha=0.35,pts=True,s=2.5,palpha=0.6,rng=None):
    bp=ax.boxplot(data,positions=pos,widths=width,whis=(5,95),showfliers=False,patch_artist=True,medianprops=dict(color='k',lw=0.8),
                  whiskerprops=dict(lw=0.5),capprops=dict(lw=0.5),boxprops=dict(lw=0.5))
    for p in bp['boxes']: p.set_facecolor(color); p.set_alpha(alpha)
    if pts:
        rng=rng or np.random.default_rng(0)
        for x0,dd in zip(pos,data): ax.scatter(x0+rng.uniform(-0.1,0.1,len(dd)),dd,s=s,color=color,alpha=palpha,lw=0,rasterized=True,zorder=3)
    return bp
