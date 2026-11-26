"""Usage (archive root): python3 scripts/05_figures.py  -> figures/*.pdf and 600-dpi *.tiff"""
import numpy as np,pandas as pd,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':8,'pdf.fonttype':42,'font.family':'DejaVu Sans'})
S=['NC1','NC2','NC3','NM1','NM2','NM3']; CN,CM='#0072B2','#D55E00'
def save(f,n):
    f.savefig(f'figures/{n}.pdf',bbox_inches='tight'); f.savefig(f'figures/{n}.tiff',dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'}); plt.close(f)
# 1 marker heatmap
L=pd.read_csv('results/scores/log2CPM_tested_genes.csv')
panel={'HSC / MPP':['Hlf','Mecom','Mpl','Kit','Cd34'],'Myeloid progenitor':['Mpo','Elane','Prtn3','Ctsg','Prss57','Ms4a3'],
 'Granulocyte':['Ly6g','Cebpe','Camp','Ltf','S100a8','S100a9','Retnlg'],'Monocyte':['Csf1r','Ly6c2','Adgre1','F13a1','Ms4a6c'],
 'Erythroid':['Gata1','Klf1','Gypa','Alas2','Slc4a1','Rhd'],'B cell':['Cd19','Cd79a','Pax5','Ms4a1','Ighd','Fcer2a'],
 'T cell':['Cd3e','Cd3d','Cd3g','Trac','Cd5'],'Mast cell':['Cma1','Mcpt4','Cpa3','Ms4a2'],'Basophil':['Mcpt8','Prss34','Fcer1a']}
rows=[];lab=[];missing=[]
for g,gs in panel.items():
    for x in gs:
        r=L[L.symbol==x]
        if len(r)!=1: missing.append((g,x)); continue
        v=r[S].values[0].astype(float); rows.append((v-v.mean())/v.std(ddof=1)); lab.append((g,x))
M=np.array(rows)
f,a=plt.subplots(figsize=(4.2,0.15*len(lab)+1.2))
im=a.imshow(M,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
a.set_xticks(range(6),S); a.set_yticks(range(len(lab)),[x for _,x in lab],style='italic'); a.axvline(2.5,color='k',lw=1.2)
prev=None
for i,(g,_) in enumerate(lab):
    if g!=prev:
        if prev: a.axhline(i-.5,color='k',lw=.8)
        n=sum(1 for h,_ in lab if h==g); a.text(5.8,i+(n-1)/2,g,va='center'); prev=g
f.colorbar(im,ax=a,orientation='horizontal',fraction=.025,pad=.05,label='row z-score of log2 CPM')
a.set_title('Selected lineage markers (display panel)'); save(f,'Figure1_marker_heatmap')
pd.DataFrame([(g,x,'shown') for g,x in lab]+[(g,x,'not tested (low expression or absent)') for g,x in missing],columns=['group','symbol','status']).to_csv('results/scores/heatmap_panel_genes.csv',index=False)
# 2 scores
s=pd.read_csv('results/scores/lineage_scores_by_sample.csv'); s=s[s.method=='mean_z'].set_index('signature')
c=pd.read_csv('results/camera_roast/lineage_sets_camera_fry.csv').set_index('signature')
order=['HSC_MPP','Myeloid_progenitor_vs_mature [post hoc]','Granulocyte_neutrophil','Monocyte','Erythroid','B_cell','T_cell','Mast_cell','Basophil']
names=['HSC / MPP','Myeloid progenitor (post hoc)','Granulocyte','Monocyte','Erythroid','B cell','T cell','Mast cell','Basophil']
f,ax=plt.subplots(3,3,figsize=(7,6.6),sharey=True)
for a,k,n in zip(ax.flat,order,names):
    r=s.loc[k]
    for j,(g,col) in enumerate([('NC',CN),('NM',CM)]):
        v=[r[f'{g}{i}'] for i in (1,2,3)]; a.scatter(np.array([j-.12,j,j+.12]),v,color=col,s=22); a.hlines(np.mean(v),j-.25,j+.25,color=col)
    a.set_xticks([0,1],['NPM1::CCDC28A','NPM1::MLF1'],fontsize=6.5); a.axhline(0,color='.7',lw=.5)
    a.set_title(f"{n} ({int(r.genes_tested)} genes)\nWelch p={r.welch_p:.3f}; fry FDR={c.loc[k].fry_FDR:.3f}",fontsize=7)
for a in ax[:,0]: a.set_ylabel('signature score (mean z)')
f.suptitle('Lineage signature scores (expression scores, not cell proportions)'); f.tight_layout(); save(f,'Figure2_lineage_scores')
# 3 marker enrichment
E=pd.read_csv('results/marker_enrichment/DE_lineage_marker_enrichment.csv')
keys=[k for k in order if k!='Myeloid_progenitor_vs_mature [post hoc]']+['Myeloid_progenitor_vs_mature [post hoc]']
f,ax=plt.subplots(1,2,figsize=(7,3.4),sharey=True)
for a,level,t in zip(ax,['mouse genes','mapped human orthologs'],['522 concordant genes (mouse)','398 mapped human orthologs']):
    for off,d,col,lab_ in [(-.18,'up_in_NM',CM,'higher in MLF1'),(.18,'up_in_NC',CN,'higher in CCDC28A')]:
        e=E[(E.level==level)&(E.direction==d)].set_index('signature').loc[keys]
        pc=100*e.overlap/e.markers_in_background.replace(0,np.nan)
        a.barh(np.arange(len(keys))+off,pc,.34,color=col,label=lab_)
        for i,(v,r) in enumerate(zip(pc,e.itertuples())):
            if r.overlap: a.text(v+1,i+off,f"{r.overlap}/{r.markers_in_background}"+('*' if r.FDR<.05 else ''),va='center',fontsize=6)
    a.set_title(t); a.set_xlabel('% of lineage markers in the DE list'); a.set_xlim(0,120)
ax[0].set_yticks(range(len(keys)),[names[order.index(k)] for k in keys]); ax[0].invert_yaxis(); ax[1].legend(frameon=False,loc='lower right')
f.text(.5,-.02,'* hypergeometric FDR < 0.05',ha='center'); f.tight_layout(); save(f,'Figure3_marker_enrichment')
# 4 sensitivity
o=pd.read_csv('results/sensitivity/ORA_summary.csv'); o=o[(o['list']=='concordant')&(o.direction=='up_in_NM')].set_index('collection')
g=pd.read_csv('results/sensitivity/GSEA_summary.csv'); cols=['GOBP','GOCC','GOMF','REACTOME']
f,ax=plt.subplots(1,2,figsize=(7,3))
x=np.arange(4); a=ax[0]
a.bar(x-.2,o.loc[cols].significant_full+.0,.38,color='.65',label='all genes'); a.bar(x+.2,o.loc[cols].significant_without_markers,.38,color='#009E73',label='markers removed')
r=o.loc[cols]; a.errorbar(x+.2,r.random_removal_median,yerr=[r.random_removal_median-r.random_removal_2_5pct,r.random_removal_97_5pct-r.random_removal_median],fmt='o',color='k',ms=3,capsize=3,label='random removal (median, 95% range)')
a.set_yscale('symlog'); a.set_xticks(x,cols); a.set_ylabel('terms with FDR < 0.05'); a.set_title('A  ORA: genes higher in MLF1',loc='left'); a.legend(frameon=False,fontsize=6)
t=g.groupby('collection')[['significant_full','significant_without_markers']].sum().loc[cols]; a=ax[1]
a.bar(x-.2,t.significant_full,.38,color='.65'); a.bar(x+.2,t.significant_without_markers,.38,color='#009E73'); a.set_yscale('symlog'); a.set_xticks(x,cols)
a.set_title('B  Preranked GSEA: both directions',loc='left'); f.tight_layout(); save(f,'Figure4_marker_exclusion_sensitivity')
# 5 correlation-aware
C=pd.read_csv('results/camera_roast/GO_Reactome_camera_fry_all_genes.csv')
f,ax=plt.subplots(1,3,figsize=(7,2.5),sharey=False)
for a,col,t in zip(ax,['camera_preset_FDR','camera_FDR','fry_FDR'],['CAMERA, preset correlation 0.01','CAMERA, estimated correlation','fry (rotation, self-contained)']):
    v=C[col].values; a.hist(v,bins=np.linspace(0,1,41),color='.5'); a.axvline(.05,color=CM,lw=1); a.set_title(f'{t}\nFDR<0.05: {(v<.05).sum()}; min {v.min():.3g}',fontsize=7); a.set_xlabel('FDR')
ax[0].set_ylabel('gene sets'); f.tight_layout(); save(f,'Figure5_correlation_aware_tests')
# 6 TF
f,ax=plt.subplots(1,2,figsize=(7,4.2))
for a,k in zip(ax,[10,20]):
    t=pd.read_csv(f'results/tf/TF_resolved_min{k}.csv').head(21).iloc[::-1]
    a.barh(t.TF+np.where(t.renamed,' †',''),t.NM_minus_NC,color=[CM if q<.05 else '.7' for q in t.FDR]); a.set_title(f"{'A' if k==10 else 'B'}  ≥{k} targets: top 21 by p",loc='left'); a.set_xlabel('score difference, MLF1 − CCDC28A')
f.text(.5,-.02,'orange: limma FDR < 0.05;  † name translated from a UniProt accession or non-mouse symbol.  Inferred activity, not measured.',ha='center',fontsize=7)
f.tight_layout(); save(f,'Figure6_TF_resolved')
print('missing panel genes',missing)
