"""Electronic supplementary material data files for the J. R. Soc. Interface submission:
Data S1 (per-axon tables for both connectomes, including the new arbor and release-site columns) and
Data S2 (numerical values behind every figure panel and table)."""
import os, json, numpy as np, pandas as pd
# run from the repository root: results/ in, data_tables/ out, public data from FLYBRAIN_DATA
HERE='data_tables'; R='results'
D=os.path.expanduser(os.environ.get('FLYBRAIN_DATA','../data'))
CODE='https://github.com/cicisoso/How-much-the-brain-tells-the-body (Zenodo https://doi.org/10.5281/zenodo.22758901)'

m=pd.read_parquet(f'{R}/channel_axons_male_arbor.parquet'); dl=pd.read_parquet(f'{R}/channel_axons_male_delay.parquet')
m=m.merge(dl[['body','len_to_brain_um','len_to_vnc_um','end_brain','end_vnc','v_ms','latency_ms']],on='body',how='left')
male=pd.DataFrame({'bodyId':m.body,'superclass':m.superclass,'type':m.type,'instance':m.instance,'somaSide':m.somaSide,'subclass_DN_target':m.subclass,'somaNeuromere_AN_origin':m.somaNeuromere,'dimorphism':m.dimorphism,
 'predicted_neurotransmitter':m.consensus_nt,'neurons_in_cell_type':m.type_n,'neck_area_um2':m.area_um2,'neck_diameter_um':m.diam_um,'diameter_cv_across_slabs':m.diam_cv,'n_slabs':m.n_slabs,
 'output_synapses_total':m.n_out,'output_synapses_VNC':m.n_out_vnc,'output_synapses_brain':m.n_out_brain,'postsynaptic_partners':m.n_out_partners,'postsynaptic_partners_ge5':m.n_out_partners5,
 'input_synapses_total':m.n_in,'input_synapses_VNC':m.n_in_vnc,'input_synapses_brain':m.n_in_brain,'presynaptic_partners':m.n_in_partners,'presynaptic_partners_ge5':m.n_in_partners5,
 'release_sites_VNC':m.tbar_vnc.astype(int),'release_sites_brain':m.tbar_brain.astype(int),
 'cable_VNC_um':m.cable_vnc_um.round(1),'cable_brain_um':m.cable_brain_um.round(1),'cable_connective_um':m.cable_neck_um.round(1),'cable_total_um':m.cable_total_um.round(1),
 'terminal_tips_VNC':m.tips_vnc,'terminal_tips_brain':m.tips_brain,'branch_points_VNC':m.bp_vnc,'branch_points_brain':m.bp_brain,
 'mean_skeleton_diameter_VNC_um':m.skel_diam_vnc_um.round(3),'mean_skeleton_diameter_brain_um':m.skel_diam_brain_um.round(3),
 'path_neck_to_brain_arbor_um':m.len_to_brain_um,'path_neck_to_vnc_arbor_um':m.len_to_vnc_um,'walk_end_brain':m.end_brain,'walk_end_vnc':m.end_vnc,'model_velocity_m_per_s':m.v_ms,'model_arbor_to_arbor_delay_ms':m.latency_ms})
f=pd.read_parquet(f'{R}/channel_axons_female.parquet')
female=pd.DataFrame({'root_888':f.root_888.astype(str),'superclass':f.superclass,'cell_type':f.cell_type,'malecns_cell_type':f.malecns_cell_type,'side':f.side,'dimorphism':f.dim,'n_planes':f.n_planes,'neck_area_um2':f.area_um2,'neck_diameter_um':f.diam_um,'relative_calibre':f.rel})
readme=pd.DataFrame({'Data S1. Per-axon calibres, classes, targets, synapses, release sites, arbors and delays for every neck-crossing neuron of the MaleCNS and the BANC connectomes':[
 'Sheet "MaleCNS": one row per brain-body axon of the male connectome (MaleCNS v1.0; 3698 channel axons). Identifiers and annotations (bodyId, superclass, type, instance, somaSide; subclass_DN_target: VNC-target code of DNs, xn multiple neuropils, xl all leg neuropils, lt lower tectulum, ut upper tectulum, fl front leg, nt neck tectulum, it intermediate tectulum, ht haltere tectulum, ad abdomen; somaNeuromere_AN_origin; dimorphism; predicted_neurotransmitter, consensus prediction of the release; neurons_in_cell_type).',
 'Calibre: neck_area_um2 (median cross-sectional area over three 16 nm slabs, tilt-corrected), neck_diameter_um (diameter of the circle of equal area), diameter_cv_across_slabs, n_slabs.',
 'Synapses (partner table, confidence >= 0.5; each row of the table is one presynaptic-postsynaptic contact): output and input contacts in total, in the VNC and in the brain (assigned by the neuropil of the postsynaptic site); numbers of postsynaptic and presynaptic partners, all and those connected by >= 5 synapses; release_sites_VNC and release_sites_brain, the number of distinct presynaptic locations (T-bars) of the axon in each region.',
 'Arbors (published skeletons): cable length in the VNC (z > 56 500 voxels), in the brain (z < 51 500 voxels), in the connective and in total; terminal tips and branch points in each region; mean skeleton diameter in each region (skeleton radii are quantized with a floor of 0.26 um and were used only for rank correlations).',
 'Delays: path lengths from the neck plane to the first branch point on each side; walk_end (branch or terminus); model_velocity_m_per_s (0.7*sqrt(d)); model_arbor_to_arbor_delay_ms.',
 'Sheet "BANC": one row per brain-body axon of the female connectome (BANC v888) with a mesh section at the curated neck plane: root_888, superclass, cell_type, malecns_cell_type (matched male type), side, dimorphism, n_planes, neck_area_um2 (median mesh-section area over three planes), neck_diameter_um, relative_calibre (diameter divided by the median of isomorphic axons of the same class in the same dataset).',
 'Source datasets: MaleCNS v1.0 (https://male-cns.janelia.org, CC BY 4.0) and BANC v888 (gs://lee-lab_brain-and-nerve-cord-fly-connectome). Code: '+CODE+'.']})
with pd.ExcelWriter(f'{HERE}/Data_S1_per_axon_tables.xlsx',engine='openpyxl') as w:
    readme.to_excel(w,sheet_name='README',index=False); male.to_excel(w,sheet_name='MaleCNS',index=False); female.to_excel(w,sheet_name='BANC',index=False)
male.to_csv(f'{HERE}/DataS1_MaleCNS_neck_axons.csv',index=False); female.to_csv(f'{HERE}/DataS1_BANC_neck_axons.csv',index=False)
print('Data_S1_per_axon_tables.xlsx', male.shape, female.shape)

# ---------------- Data S2: one sheet per figure panel or table
src={}; ann=pd.read_feather(f'{D}/body-annotations-male-cns-v1.0-minconf-0.5.feather'); cal=pd.read_parquet(f'{R}/caliber_per_body.parquet'); ch=pd.read_parquet(f'{R}/channel_axons_male.parquet')
lst=cal[cal.neck_listed&cal.type.notna()&cal.somaSide.isin(['L','R'])]; t=lst.groupby(['type','somaSide']).diam_um.median().unstack().dropna()
src['Fig1c']=t.reset_index().rename(columns={'L':'diameter_left_um','R':'diameter_right_um'})
mv=pd.read_parquet(f'{R}/mip_validation.parquet'); src['Fig1d']=mv.reset_index()[['body','area_um2_raw','area_um2_mip0']].rename(columns={'body':'bodyId','area_um2_raw':'area_16nm_um2','area_um2_mip0':'area_8nm_um2'})
src['Fig2a']=ch[['body','superclass','diam_um']].rename(columns={'body':'bodyId'})
a=np.sort(ch.area_um2.values)[::-1]; src['Fig2b']=pd.DataFrame({'fraction_of_axons':np.arange(1,len(a)+1)/len(a),'cumulative_area_fraction':np.cumsum(a)/a.sum()})
dn=ch[ch.dir=='down']; up=ch[ch.dir=='up']
src['Fig2c']=pd.DataFrame({'measure':['axon count','summed diameter','summed area'],'descending':[len(dn),dn.diam_um.sum(),dn.area_um2.sum()],'ascending':[len(up),up.diam_um.sum(),up.area_um2.sum()]})
def shares(df,key): g=df.groupby(key).agg(n=('body','size'),area_um2=('area_um2','sum')); g['count_share']=g.n/g.n.sum(); g['area_share']=g.area_um2/g.area_um2.sum(); return g.reset_index()
src['Fig2d']=shares(ch[ch.superclass=='descending_neuron'].fillna({'subclass':'other'}),'subclass'); src['Fig2e']=shares(ch[ch.superclass=='ascending_neuron'].fillna({'somaNeuromere':'other'}),'somaNeuromere'); src['Fig2f']=shares(ch[ch.superclass=='sensory_ascending'].fillna({'subclass':'other'}),'subclass')
mm=m[m.diam_um>=0.1]
src['Fig3a']=mm[mm.superclass=='descending_neuron'][['body','type','diam_um','n_out_vnc']].rename(columns={'body':'bodyId'})
src['Fig3b']=mm[mm.superclass=='ascending_neuron'][['body','type','diam_um','n_out_brain']].rename(columns={'body':'bodyId'})
src['Fig3c']=mm[mm.superclass=='sensory_ascending'][['body','type','diam_um','n_out_brain']].rename(columns={'body':'bodyId'})
de=mm[mm.superclass.isin(['descending_neuron','ascending_neuron'])].copy()
de['target_cable_um']=np.where(de.superclass=='descending_neuron',de.cable_vnc_um,de.cable_brain_um); de['target_release_sites']=np.where(de.superclass=='descending_neuron',de.tbar_vnc,de.tbar_brain)
de['target_output_synapses']=np.where(de.superclass=='descending_neuron',de.n_out_vnc,de.n_out_brain); de['release_sites_per_um']=de.target_release_sites/de.target_cable_um
src['Fig3d_e']=de[['body','superclass','diam_um','target_cable_um','target_release_sites','release_sites_per_um','target_output_synapses']].rename(columns={'body':'bodyId'})
A=json.load(open(f'{R}/arbor_scaling.json')); rows=[]
for sc in ['descending_neuron','ascending_neuron','sensory_ascending']:
    for k,lab in [('S','output'),('L','target cable'),('T_per_L','release sites per um'),('S_per_T','contacts per release site')]:
        s=A[sc]['decomposition'][k]['slope']; rows.append({'class':sc,'quantity':lab,'exponent':s['coef'],'ci_low':s['ci'][0],'ci_high':s['ci'][1],'n':A[sc]['decomposition'][k]['n']})
src['Fig3f']=pd.DataFrame(rows)
src['Fig4a']=mm[mm.superclass.isin(['descending_neuron','ascending_neuron'])][['body','superclass','diam_um','n_out_partners5']].rename(columns={'body':'bodyId'})
src['Fig4b']=pd.concat([pd.read_csv(f'{R}/type_population_{sc}.csv').assign(superclass=sc) for sc in ['descending_neuron','ascending_neuron']])
src['Fig4c']=mm[mm.superclass.isin(['descending_neuron','ascending_neuron'])&mm.consensus_nt.isin(['acetylcholine','gaba','glutamate'])][['body','superclass','consensus_nt','diam_um']].rename(columns={'body':'bodyId'})
j=pd.read_parquet(f'{R}/type_matched_relative.parquet'); j.index.name='type'; src['Fig5a']=j.reset_index()[['type','male_abs','female_abs']].rename(columns={'male_abs':'male_diameter_um','female_abs':'female_diameter_um'})
J=json.load(open(f'{R}/compare_sex2.json')); src['Fig5b']=pd.DataFrame({'measure':['axon count','summed diameter','summed area'],'male_descending_share':[J['male_share_count'],J['male_share_sumd'],J['male_share_area']],'female_descending_share':[J['female_share_count'],J['female_share_sumd'],J['female_share_area']]})
fm=pd.read_parquet(f'{R}/channel_axons_male.parquet'); ff=pd.read_parquet(f'{R}/channel_axons_female.parquet')
src['Fig5c']=pd.concat([fm[fm.superclass=='descending_neuron'][['body','dimcat','rel']].assign(sex='male').rename(columns={'body':'id'}),ff[ff.superclass=='descending_neuron'][['root_888','dimcat','rel']].assign(sex='female').rename(columns={'root_888':'id'})])
src['Fig5c']['id']=src['Fig5c'].id.astype(str)
src['Fig5d']=j[j.dim=='sexually dimorphic'].reset_index()[['type','male','female']].rename(columns={'male':'relative_calibre_male','female':'relative_calibre_female'})
dl2=dl[(dl.end_brain=='branch')&(dl.end_vnc=='branch')&dl.superclass.isin(['descending_neuron','ascending_neuron'])]; src['Fig6a']=dl2[['body','superclass','diam_um','L_um','v_ms','latency_ms']].rename(columns={'body':'bodyId'})
K=json.load(open(f'{R}/capacity_sexes.json')); src['Fig6b']=pd.DataFrame({'beta':K['betas'],'ratio_male':K['male']['ratio_curve'],'ratio_female':K['female']['ratio_curve']})
order=m.sort_values('diam_um',ascending=False); ks=np.concatenate([[0],np.unique(np.round(np.logspace(0,3,40)).astype(int))]); rows=[]
for k in ks:
    rest=order.iloc[k:]; dd=rest[rest.dir=='down']; uu=rest[rest.dir=='up']; rows.append({'thickest_removed':int(k),'capacity_ratio_beta1':dd.diam_um.sum()/uu.diam_um.sum(),'descending_area_share':dd.area_um2.sum()/rest.area_um2.sum()})
src['Fig6c']=pd.DataFrame(rows)
src['Table1']=pd.read_csv(f'{R}/budget_table.csv')
g=pd.read_parquet(f'{R}/caliber_per_slab.parquet'); gg=g[g.body.isin(ch.body)&(g.n_planes>=32)].pivot(index='body',columns='slab',values='diam_um').dropna(); src['FigS1a']=gg.reset_index().rename(columns={'A':'diameter_anterior_um','B':'diameter_middle_um','C':'diameter_posterior_um'})
src['FigS1b']=g[g.body.isin(ch.body)][['body','slab','cos_theta']].dropna()
mvv=pd.read_parquet(f'{R}/male_mesh_validation.parquet').dropna(); src['FigS1c']=mvv[['body','area_um2','area_mesh_um2']]
dlb=pd.read_parquet(f'{D}/neck_delay_budget.parquet'); src['FigS1d']=ch.merge(dlb[['bodyId','r_skel_neck']],on='bodyId')[['body','r_skel_neck','diam_um']].assign(skeleton_diameter_um=lambda x: x.r_skel_neck*2*0.008)
rng=np.random.default_rng(0); di=dn.diam_um.values; ui=up.diam_um.values; src['FigS2a']=pd.DataFrame({'bootstrap_ratio_beta1':[rng.choice(di,len(di)).sum()/rng.choice(ui,len(ui)).sum() for i in range(2000)]})
src['FigS2b']=mm[mm.superclass.isin(['descending_neuron','ascending_neuron'])][['body','superclass','diam_um','area_um2','n_out_vnc','n_out_brain']].rename(columns={'body':'bodyId'})
src['FigS2c']=src['Fig3d_e']
src['FigS2d']=mm[mm.superclass=='sensory_ascending'][['body','diam_um','cable_brain_um','n_out_brain']].rename(columns={'body':'bodyId'})
with pd.ExcelWriter(f'{HERE}/Data_S2_figure_source_data.xlsx',engine='openpyxl') as w:
    pd.DataFrame({'Data S2':['Numerical values behind every panel of figures 1-6, table 1 and electronic supplementary material figures S1-S2 of "Scaling laws of axon calibre and the capacity of the brain-body channel in Drosophila". Sheet names give the figure and panel. Figures 1a and 1b are images (skeleton density and a segmentation plane) without numerical source data beyond the public datasets. Code: '+CODE+'.']}).to_excel(w,sheet_name='README',index=False)
    for k,v in src.items(): v.to_excel(w,sheet_name=k,index=False)
print('Data_S2_figure_source_data.xlsx sheets:',list(src))
