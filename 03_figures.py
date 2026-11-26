"""Usage (from archive root): python3 scripts/03_figures.py  -> figures/Figure1-5 .pdf and 600-dpi .tiff"""
import textwrap,numpy as np,pandas as pd,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.cluster.hierarchy import linkage,dendrogram,leaves_list
from scipy.spatial.distance import squareform
plt.rcParams.update({'font.size':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','axes.titlesize':9})
NM,NC,GR='#D55E00','#0072B2','#BBBBBB'; six=['NM1','NM2','NM3','NC1','NC2','NC3']
def save(f,name):
    f.savefig(f'figures/{name}.pdf',bbox_inches='tight')
    f.savefig(f'figures/{name}.tiff',dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'}); plt.close(f)
def letter(ax,l): ax.set_title(l,loc='left')
d=pd.read_csv('results/de/NM_vs_NC_all_genes_annotated.csv'); lc=pd.read_csv('results/qc/log2cpm.csv',index_col=0)

# Figure 1
f=plt.figure(figsize=(7.2,6.6)); gs=f.add_gridspec(2,1,height_ratios=[1,1.1],hspace=.35)
ax=f.add_subplot(gs[0]); ax.axis('off'); ax.set_xlim(0,10); ax.set_ylim(-.3,6); letter(ax,'A  Revised analytical workflow')
def box(x,y,w,h,s,c='#F2F2F2'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.08',fc=c,ec='#444444',lw=.8)); ax.text(x+w/2,y+h/2,s,ha='center',va='center',fontsize=7.5)
def arrow(a,b): ax.annotate('',b,a,arrowprops=dict(arrowstyle='->',lw=.8))
box(3,5,4,.7,'GSE274558 raw counts\n9 samples, 55,402 genes'); box(.2,3.5,4,.8,'Control vs fusion (9 samples)\nconfounded with inferred sex: exploratory only','#FBE3D6')
box(5.8,3.5,4,.8,'NM vs NC (6 fusion samples)\ninferred sex constant','#DCEBF7'); arrow((4,4.9),(2.4,4.4)); arrow((6,4.9),(7.6,4.4))
box(5.8,2.3,4,.7,'QC; DESeq2 and edgeR; filter sensitivity'); box(5.8,1.1,4,.7,'Ensembl 116 orthologs; GO, Reactome, GSEA'); box(5.8,-.05,4,.7,'CollecTRI-based TF activity (inferred)')
for a,b in [(3.4,3.08),(2.2,1.88),(1.0,.73)]: arrow((7.8,a),(7.8,b))
full=pd.read_csv('inputs/GSE274558_raw_counts.txt.gz',sep='\t',index_col=0)
ids={'Xist':'086503','Kdm5d':'056673','Ddx3y':'069045','Uty':'068457','Eif2s3y':'069049'}
cpm=full/full.sum()*1e6; m=np.log2(cpm.loc[['ENSMUSG00000'+i for i in ids.values()]]+1)
ax=f.add_subplot(gs[1]); im=ax.imshow(m.values,aspect='auto',cmap='viridis')
ax.set_xticks(range(9),[c.replace('WT','WT') for c in m.columns]); ax.set_yticks(range(5),list(ids),style='italic')
for i in range(5):
    for j in range(9): ax.text(j,i,f'{cpm.loc[m.index[i]].iloc[j]:.1f}',ha='center',va='center',fontsize=6.5,color='white' if m.values[i,j]<m.values.max()*.6 else 'black')
for x in (2.5,5.5): ax.axvline(x,color='white',lw=3)
letter(ax,'B  Sex-associated genes (cell labels: CPM)'); f.colorbar(im,ax=ax,label='log2(CPM+1)',shrink=.8)
f.text(.12,.015,'WT1-3 are file labels (publication describes vector-transduced controls; not verified). Inferred sex: WT male, NM and NC female.\nDocumented sex not determined. ~ sex + group on all nine samples: rank 3 of 4 columns (not estimable).',fontsize=7)
save(f,'Figure1_workflow_sex_audit')

# Figure 2
f,A=plt.subplots(2,2,figsize=(7.2,6.8)); f.subplots_adjust(hspace=.45,wspace=.4)
p=pd.read_csv('results/qc/pca.csv'); ax=A[0,0]
for _,r in p.iterrows():
    nm=r['sample'].startswith('NM'); ax.scatter(r.PC1,r.PC2,c=NM if nm else NC,marker='o' if nm else 's',s=40); ax.annotate(r['sample'],(r.PC1,r.PC2),xytext=(3,3),textcoords='offset points',fontsize=7)
ax.set_xlabel(f'PC1 ({p.PC1_var[0]:.1%})'); ax.set_ylabel(f'PC2 ({p.PC2_var[0]:.1%})'); letter(ax,'A  PCA')
c=pd.read_csv('results/qc/spearman.csv',index_col=0); Z=linkage(squareform(1-c.values,checks=False),'average'); o=leaves_list(Z); ax=A[0,1]
im=ax.imshow(c.values[np.ix_(o,o)],cmap='viridis'); ax.set_xticks(range(6),c.index[o],rotation=90); ax.set_yticks(range(6),c.index[o]); f.colorbar(im,ax=ax,label='Spearman r',shrink=.8)
ax.set_title('B  Correlation, clustered',loc='left',pad=24); top=ax.inset_axes([0,1.01,1,.13]); dendrogram(Z,ax=top,no_labels=True,color_threshold=0,above_threshold_color='black'); top.axis('off')
ax=A[1,0]; y=-np.log10(d.edgeR_FDR.clip(lower=1e-300))
for mask,col,lab in [(~d.concordant,GR,'other'),(d.concordant&(d.edgeR_logFC>0),NM,'higher in NM'),(d.concordant&(d.edgeR_logFC<0),NC,'higher in NC')]:
    ax.scatter(d.edgeR_logFC[mask],y[mask],s=3,c=col,label=f'{lab} ({mask.sum():,})',rasterized=True)
ax.axhline(-np.log10(.05),ls='--',lw=.6,c='k'); ax.set_xlabel('edgeR log2 fold change (NM / NC)'); ax.set_ylabel('-log10 FDR (edgeR)'); ax.legend(frameon=False,fontsize=6.5); letter(ax,'C  Differential expression')
s=pd.read_csv('results/de/summary_counts.csv'); s=s[s['filter'].str.startswith('filterBy')]; ax=A[1,1]; x=np.arange(3)
for k,(col,lab,colr) in enumerate([('sig','FDR < 0.05','#555555'),('sig_lfc1','and |log2FC| >= 1','#AAAAAA')]):
    b=ax.bar(x+(k-.5)*.38,s[col],.38,label=lab,color=colr); ax.bar_label(b,fontsize=6.5)
ax.set_xticks(x,['edgeR','DESeq2','Concordant']); ax.set_ylabel('Genes'); ax.legend(frameon=False,fontsize=6.5); ax.set_ylim(0,s.sig.max()*1.2); letter(ax,'D  Significant genes')
save(f,'Figure2_QC_differential_expression')

# Figure 3
f=plt.figure(figsize=(7.2,8.6)); gs=f.add_gridspec(1,2,width_ratios=[1.15,1],wspace=.55); ax=f.add_subplot(gs[0]); ax.set_box_aspect(1)
for mask,col,lab in [(~d.concordant,GR,'not concordant'),(d.concordant,'#222222','concordant')]:
    ax.scatter(d.edgeR_logFC[mask],d.DESeq2_log2FC[mask],s=3,c=col,label=lab,rasterized=True)
lim=np.nanmax(np.abs(d[['edgeR_logFC','DESeq2_log2FC']].values))*1.05; ax.plot([-lim,lim],[-lim,lim],lw=.6,c=NM); ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
r=d[['edgeR_logFC','DESeq2_log2FC']].corr().iloc[0,1]; ax.set_xlabel('edgeR log2FC'); ax.set_ylabel('DESeq2 log2FC'); ax.legend(frameon=False,loc='upper left'); letter(ax,f'A  Method concordance (r = {r:.3f})')
k=d[d.concordant_lfc1].sort_values(['edgeR_FDR','gene']); g=pd.concat([k[k.edgeR_logFC>0].head(25),k[k.edgeR_logFC<0].head(25)])
M=lc.loc[g.gene,six].values; M=(M-M.mean(1,keepdims=True))/M.std(1,ddof=1,keepdims=True)
ax=f.add_subplot(gs[1]); im=ax.imshow(M,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2); ax.set_xticks(range(6),six,rotation=90)
ax.set_yticks(range(len(g)),g.mouse_symbol.fillna(g.gene),fontsize=5.8,style='italic'); ax.axhline(24.5,color='black',lw=1); ax.axvline(2.5,color='black',lw=1)
letter(ax,'B  Top 25 per direction'); f.colorbar(im,ax=ax,label='Row z-score',shrink=.4)
save(f,'Figure3_concordance_heatmap')

# Figure 4
o=pd.read_csv('results/enrichment/ORA_all_terms.csv'); o=o[(o['list']=='concordant_FDR')&(o.FDR<.05)]
G=pd.read_csv('results/enrichment/GSEA_all_terms.csv'); G=G[G['FDR q-val']<.05]
def name(t,n=48):
    t=t.split('_',1)[1].replace('_',' ').capitalize(); return textwrap.shorten(t,n,placeholder='...')
colls=['GOBP','GOMF','GOCC','REACTOME']
f,A=plt.subplots(4,2,figsize=(10.5,11),gridspec_kw=dict(wspace=1.25,hspace=.55))
for i,c in enumerate(colls):
    ax=A[i,0]; sub=pd.concat([o[(o.collection==c)&(o.direction==dr)].nsmallest(4,'FDR') for dr in ('up_in_NM','up_in_NC')])
    if len(sub):
        v=np.where(sub.direction=='up_in_NM',1,-1)*-np.log10(sub.FDR); ax.barh(range(len(sub)),v,color=np.where(v>0,NM,NC))
        ax.set_yticks(range(len(sub)),[name(t) for t in sub.term],fontsize=6.3); ax.invert_yaxis(); ax.axvline(0,c='k',lw=.6)
    else: ax.text(.5,.5,'No term with FDR < 0.05',ha='center',transform=ax.transAxes)
    ax.set_xlabel('signed -log10 FDR (left: higher in NC)'); letter(ax,f'{"ABCD"[i]}  {c} over-representation')
    ax=A[i,1]; s=G[G.collection==c]; sub=pd.concat([s[s.NES>0].nlargest(4,'NES'),s[s.NES<0].nsmallest(4,'NES')])
    if len(sub):
        ax.barh(range(len(sub)),sub.NES,color=np.where(sub.NES>0,NM,NC)); ax.set_yticks(range(len(sub)),[name(t) for t in sub.Term],fontsize=6.3); ax.invert_yaxis(); ax.axvline(0,c='k',lw=.6)
    else: ax.text(.5,.5,'No set with FDR q < 0.05',ha='center',transform=ax.transAxes)
    ax.set_xlabel('Normalized enrichment score (positive: NM)'); letter(ax,f'{"EFGH"[i]}  {c} ranked GSEA')
f.text(.1,.045,'Orange: higher in NM. Blue: higher in NC. Up to four terms per direction, FDR < 0.05. MSigDB 7.5.1, human one-to-one orthologs.',fontsize=7)
save(f,'Figure4_enrichment')

# Figure 5
a=pd.read_csv('results/tf/TF_NM_vs_NC_min10.csv').set_index('TF'); b=pd.read_csv('results/tf/TF_NM_vs_NC_min20.csv').set_index('TF')
sel=[t for t in a.index.union(b.index) if ((t in a.index and a.FDR[t]<.05) or (t in b.index and b.FDR[t]<.05)) and a.name_is_mouse_symbol[t]]
sel=sorted(sel,key=lambda t:a.NM_minus_NC[t])
f,ax=plt.subplots(figsize=(6,1.8+.24*len(sel)))
for y_,t in enumerate(sel):
    ax.plot([0,a.NM_minus_NC[t]],[y_,y_],c=GR,lw=.8,zorder=1)
    for tab,mk,off,lab in [(a,'o',.13,'min 10'),(b,'^',-.13,'min 20')]:
        if t in tab.index:
            col=NM if tab.NM_minus_NC[t]>0 else NC
            ax.scatter(tab.NM_minus_NC[t],y_+off,marker=mk,s=28,facecolors=col if tab.FDR[t]<.05 else 'white',edgecolors=col,zorder=2)
ax.set_yticks(range(len(sel)),[f'{t} ({a.targets_used[t]})' for t in sel],style='italic'); ax.axvline(0,c='k',lw=.7); ax.set_xlabel('Activity score difference, NM minus NC')
h=[plt.Line2D([],[],ls='',marker=m_,mfc=c_,mec='k',label=l) for m_,c_,l in [('o','k','min 10, FDR<0.05'),('o','white','min 10, not significant'),('^','k','min 20, FDR<0.05'),('^','white','min 20, not significant')]]
ax.legend(handles=h,frameon=False,fontsize=6.5,loc='lower right'); letter(ax,'Inferred transcription factor activity')
f.text(.02,-.02,'Number of targets in brackets. Inferred from target expression, not measured activity. TF names that are not\nmouse gene symbols are omitted here and listed in the result tables.',fontsize=6.5)
save(f,'Figure5_TF_activity')
print('TFs plotted',len(sel))
