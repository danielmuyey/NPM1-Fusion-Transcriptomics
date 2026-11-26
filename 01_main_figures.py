"""Usage: python3 scripts/01_main_figures.py <folder holding the four analysis archives>
Draws Figures 1-5 and the graphical abstract (PDF + 600-dpi TIFF) from the audited result tables."""
import sys,numpy as np,pandas as pd,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage,leaves_list
R=sys.argv[1]+'/'; D=R+'GSE274558_Direct_Fusion_Comparison/'; L=R+'GSE274558_Lineage_and_Cell_Composition_Analysis/results/'; X=R+'AML_Cross_Dataset_Integration_Analysis/results/'
plt.rcParams.update({'font.family':'Liberation Sans','font.size':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':8.5,'axes.linewidth':.6,'legend.frameon':False})
plt.rcParams.update({'mathtext.fontset':'custom','mathtext.rm':'Liberation Sans','mathtext.it':'Liberation Sans:italic','mathtext.bf':'Liberation Sans:bold','mathtext.default':'regular','mathtext.cal':'Liberation Sans'})
CM,CC,GR,INK='#eb6834','#2a78d6','#c8c7c2','#333333'
NM,NC=['NM1','NM2','NM3'],['NC1','NC2','NC3']; S=NM+NC
NAME=lambda s:s.replace('NM','MLF1 ').replace('NC','CCDC28A ').replace('WT','VEC ')
def save(f,n,d='figures'):
    f.savefig(f'{d}/{n}.pdf',bbox_inches='tight',metadata={'CreationDate':None}); f.savefig(f'{d}/{n}.tiff',dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'})
    f.savefig(f'metadata/preview_{n}.png',dpi=150,bbox_inches='tight'); plt.close(f)
def title(a,l,t): a.set_title(f'{l}  {t}',loc='left',fontweight='bold')
def strip(a,vals,ylabel):
    for i,(g,c,m) in enumerate([(NM,CM,'o'),(NC,CC,'s')]):
        v=np.array([vals[s] for s in g],float); a.scatter(i+np.array([-.12,0,.12]),v,color=c,marker=m,s=22,zorder=3,edgecolor='white',lw=.5); a.hlines(v.mean(),i-.28,i+.28,color=INK,lw=1)
    a.set_xticks([0,1],['MLF1','CCDC28A'],fontsize=6.5); a.set_xlim(-.5,1.5); a.set_ylabel(ylabel,fontsize=7)
def box(a,x,y,w,h,s,c='#f2f2f0',fs=7.5):
    a.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.08',fc=c,ec='#555555',lw=.7)); a.text(x+w/2,y+h/2,s,ha='center',va='center',fontsize=fs)
def arrow(a,p,q): a.annotate('',q,p,arrowprops=dict(arrowstyle='->',lw=.8,color=INK))
d=pd.read_csv(D+'results/de/NM_vs_NC_all_genes_annotated.csv'); lc=pd.read_csv(D+'results/qc/log2cpm.csv',index_col=0)

# Figure 1
f=plt.figure(figsize=(7.1,5.4)); gs=f.add_gridspec(2,1,height_ratios=[1,.9],hspace=.3)
a=f.add_subplot(gs[0]); a.axis('off'); a.set_xlim(0,10); a.set_ylim(0,5); title(a,'A','Study design')
box(a,2.6,4,4.8,.7,'GSE274558: GFP-positive marrow before leukemia\n9 libraries (3 vector, 3 MLF1, 3 CCDC28A)')
box(a,.1,2.2,4.2,1,'Fusion vs vector control\nExcluded from conclusions:\ninferred sex coincides with group','#f6e3da')
box(a,5.7,2.2,4.2,1,'MLF1 vs CCDC28A (3 vs 3)\nPrimary analysis','#dde9f8')
arrow(a,(3.6,3.9),(2.4,3.3)); arrow(a,(6.4,3.9),(7.6,3.3))
box(a,5.7,.3,4.2,1.2,'Differential expression; HOX genes\nLineage markers; gene sets\nFOXM1 target score (exploratory)'); arrow(a,(7.8,2.1),(7.8,1.6))
box(a,.1,.3,4.2,1.2,'External context (negative or untestable)\nGSE308825: drug exposure, 2 xenografts\nGSE248597: spatial, 8 biopsies'); arrow(a,(5.6,.9),(4.4,.9))
x=pd.read_csv(D+'inputs/GSE274558_raw_counts.txt.gz',sep='\t',index_col=0)
ids={'Xist':'086503','Kdm5d':'056673','Ddx3y':'069045','Uty':'068457','Eif2s3y':'069049'}
cols=['WT1','WT2','WT3']+S; cpm=(x/x.sum()*1e6).loc[['ENSMUSG00000'+i for i in ids.values()],cols]; m=np.log2(cpm+1)
a=f.add_subplot(gs[1]); im=a.imshow(m.values,aspect='auto',cmap='Blues'); a.set_xticks(range(9),[NAME(c) for c in cols],fontsize=6.5,rotation=35,ha='right'); a.set_yticks(range(5),list(ids),style='italic')
for i in range(5):
    for j in range(9): a.text(j,i,f'{cpm.values[i,j]:.1f}',ha='center',va='center',fontsize=6,color='white' if m.values[i,j]>m.values.max()*.55 else INK)
for v in (2.5,5.5): a.axvline(v,color='white',lw=3)
title(a,'B','Sex-associated transcripts (counts per million)'); f.colorbar(im,ax=a,label='log$_2$(CPM+1)',shrink=.8); save(f,'Figure1')

# Figure 2
f=plt.figure(figsize=(7.1,6.6),layout='constrained'); gs=f.add_gridspec(2,2,width_ratios=[1.25,1])
p=pd.read_csv(D+'results/qc/pca.csv'); a=f.add_subplot(gs[0,0])
for g,c,mk,lab in [(NM,CM,'o','MLF1'),(NC,CC,'s','CCDC28A')]:
    q=p[p['sample'].isin(g)]; a.scatter(q.PC1,q.PC2,color=c,marker=mk,s=36,label=lab)
    for r in q.itertuples(): a.annotate(r.sample[-1],(r.PC1,r.PC2),xytext=(4,3),textcoords='offset points',fontsize=6)
a.set(xlabel=f'PC1 ({p.PC1_var.iloc[0]*100:.1f}%)',ylabel=f'PC2 ({p.PC2_var.iloc[0]*100:.1f}%)'); a.legend(); title(a,'A','Principal components')
a=f.add_subplot(gs[1,0]); y=-np.log10(d.edgeR_p); c=d.concordant
a.scatter(d.edgeR_logFC[~c],y[~c],s=3,color=GR,rasterized=True,label='Not concordant')
for sel,col,lab in [(c&(d.edgeR_logFC>0),CM,'Higher in MLF1 (349)'),(c&(d.edgeR_logFC<0),CC,'Higher in CCDC28A (173)')]: a.scatter(d.edgeR_logFC[sel],y[sel],s=5,color=col,label=lab,rasterized=True)
a.set(xlabel='log$_2$ fold change (MLF1 vs CCDC28A)',ylabel='-log$_{10}$ P (edgeR)'); a.set_ylim(top=y.max()*1.4); a.legend(fontsize=6,loc='upper left'); title(a,'B','Differential expression')
cd=d[d.concordant]; top=pd.concat([cd[cd.edgeR_logFC>0].nsmallest(20,'edgeR_FDR'),cd[cd.edgeR_logFC<0].nsmallest(20,'edgeR_FDR')])
M=lc.loc[top.gene,S].values; M=(M-M.mean(1,keepdims=True))/M.std(1,ddof=1,keepdims=True)
a=f.add_subplot(gs[:,1]); im=a.imshow(M,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2); a.set_xticks(range(6),[NAME(s) for s in S],rotation=90,fontsize=6.5)
a.set_yticks(range(len(top)),top.mouse_symbol.fillna(top.gene),style='italic',fontsize=5.8); a.axvline(2.5,color=INK,lw=1); a.axhline(19.5,color=INK,lw=1)
title(a,'C','Leading genes'); f.colorbar(im,ax=a,label='row z-score',orientation='horizontal',shrink=.7); save(f,'Figure2')

# Figure 3
hox=d[d.mouse_symbol.fillna('').str.match(r'^Hox[ab]\d+$')&d.human_symbol.notna()].copy()
hox['k']=hox.mouse_symbol.str[3]+hox.mouse_symbol.str[4:].str.zfill(2); hox=hox.sort_values('k'); assert len(hox)==17
f=plt.figure(figsize=(7.1,3.9),layout='constrained'); gs=f.add_gridspec(2,3,width_ratios=[1.5,1,1])
a=f.add_subplot(gs[:,0])
for i,r in enumerate(hox.itertuples()):
    a.hlines(i,r.edgeR_logFC,0,color=GR,lw=1); a.scatter(r.edgeR_logFC,i,s=30,facecolor=CC if r.concordant else 'white',edgecolor=CC,zorder=3)
a.set_yticks(range(17),hox.mouse_symbol,style='italic'); a.invert_yaxis(); [t.set_fontweight('bold') for t in a.get_yticklabels() if t.get_text()=='Hoxb6']; a.axvline(0,color=INK,lw=.7); a.set_xlabel('log$_2$ fold change (MLF1 vs CCDC28A)')
a.legend(handles=[Line2D([],[],marker='o',ls='',mfc=CC,mec=CC,label='FDR<0.05, both methods'),Line2D([],[],marker='o',ls='',mfc='white',mec=CC,label='Not met')],fontsize=6,loc='upper center',bbox_to_anchor=(.5,-.13),ncol=2); title(a,'A','Tested Hox genes')
full=x[S]; lcpm=np.log2(full/full.sum()*1e6+1)
for k,g in enumerate(['Hoxa13','Hoxb6','Hoxb8','Hoxb9']):
    a=f.add_subplot(gs[k//2,1+k%2]); strip(a,lcpm.loc[d[d.mouse_symbol==g].gene.iloc[0]],'log$_2$(CPM+1)'); a.set_ylim(bottom=-.15)
    a.set_title(('B  ' if k==0 else '')+g,loc='left',fontstyle='italic',fontweight='bold')
save(f,'Figure3')

# Figure 4
sc=pd.read_csv(L+'scores/lineage_scores_by_sample.csv'); sc=sc[sc.method=='mean_z'].set_index('signature')
fry=pd.read_csv(L+'camera_roast/lineage_sets_camera_fry.csv').set_index('signature')
lin=[('Erythroid','Erythroid'),('B_cell','B cell'),('Basophil','Basophil'),('Mast_cell','Mast cell'),('T_cell','T cell'),('Monocyte','Monocyte'),('Granulocyte_neutrophil','Granulocyte'),('HSC_MPP','HSC/MPP')]
f=plt.figure(figsize=(7.1,6.2),layout='constrained'); sf=f.subfigures(2,1,height_ratios=[1.35,1])
ax=sf[0].subplots(2,4,sharey=True); sf[0].suptitle('A  Marker scores per sample',x=.01,ha='left',fontweight='bold',fontsize=8.5)
for a,(k,n) in zip(ax.flat,lin):
    strip(a,sc.loc[k],'mean z-score' if a in ax[:,0] else ''); a.set_title(f'{n} ({int(sc.genes_tested[k])} genes)\nfry FDR {fry.fry_FDR[k]:.3f}',fontsize=7)
E=pd.read_csv(L+'marker_enrichment/DE_lineage_marker_enrichment.csv'); E=E[E.level=='mouse genes']
a=sf[1].subplots()
for off,dr,col,lab in [(-.2,'up_in_NM',CM,'Higher in MLF1'),(.2,'up_in_NC',CC,'Higher in CCDC28A')]:
    q=E[E.direction==dr].set_index('signature').loc[[k for k,_ in lin]]; v=100*q.overlap/q.markers_in_background
    a.bar(np.arange(8)+off,v,width=.36,color=col,label=lab)
    for i,(val,fd,o,n) in enumerate(zip(v,q.FDR,q.overlap,q.markers_in_background)): a.text(i+off,val+1,f'{o}/{n}'+('*' if fd<.05 else ''),ha='center',fontsize=5.8,rotation=90,va='bottom')
a.set_xticks(range(8),[n for _,n in lin],fontsize=6.5); a.set_ylabel('markers among concordant genes (%)'); a.set_ylim(0,130); a.set_yticks([0,25,50,75,100]); a.legend(fontsize=6.5); title(a,'B','Markers within the differential gene list'); save(f,'Figure4')

# Figure 5
tf=pd.read_csv(L+'tf/TF_resolved_min10.csv').sort_values('NM_minus_NC').reset_index(drop=True); assert len(tf)==354
v=pd.read_csv(X+'foxm1_GSE274558_sample_values.csv').set_index('measure')
f=plt.figure(figsize=(7.1,2.9),layout='constrained'); gs=f.add_gridspec(1,3,width_ratios=[1.9,1,1])
a=f.add_subplot(gs[0]); sig=tf.FDR<.05
a.scatter(tf.index[~sig],tf.NM_minus_NC[~sig],s=5,color=GR,label='FDR≥0.05'); a.scatter(tf.index[sig],tf.NM_minus_NC[sig],s=14,color=INK,label=f'FDR<0.05 ({sig.sum()})')
i=tf.index[tf.TF=='Foxm1'][0]; a.annotate('Foxm1',(i,tf.NM_minus_NC[i]),xytext=(i-120,tf.NM_minus_NC[i]+.15),style='italic',arrowprops=dict(arrowstyle='-',lw=.6))
a.axhline(0,color=INK,lw=.6); a.set(xlabel='354 regulators, ranked',ylabel='score difference (MLF1 - CCDC28A)'); a.legend(fontsize=6,loc='upper left'); title(a,'A','Regulator target scores')
a=f.add_subplot(gs[1]); strip(a,v.loc['Foxm1_mRNA_log2CPM'],'Foxm1 log$_2$ CPM'); a.set_ylim(0,6.5); title(a,'B','Foxm1 mRNA')
a=f.add_subplot(gs[2]); strip(a,v.loc['Foxm1_targets_signed'],'mean signed z-score'); a.axhline(0,color=GR,lw=.7); title(a,'C','Target score'); save(f,'Figure5')

# Graphical abstract
f,a=plt.subplots(figsize=(7,3.6)); a.axis('off'); a.set_xlim(0,10); a.set_ylim(0,5.2)
a.text(5,5,'NPM1::MLF1 vs NPM1::CCDC28A in preleukemic mouse marrow (3 vs 3)',ha='center',fontsize=10.5,fontweight='bold')
a.text(5,4.55,'Public RNA-sequencing reanalysis: 522 concordant differential genes',ha='center',fontsize=9)
box(a,.15,2.2,4.55,1.8,'Higher with NPM1::MLF1 (349 genes)\n\nErythroid and B-cell marker genes\nFOXM1 target score (exploratory;\nFoxm1 mRNA unchanged)','#fbe6dd',8.5)
box(a,5.3,2.2,4.55,1.8,'Higher with NPM1::CCDC28A (173 genes)\n\nHoxb6: retained under all expression filters\nHoxa13, Hoxb8, Hoxb9: filter-sensitive\nAll four: zero or near-zero MLF1 counts','#dde9f8',8.5)
box(a,.15,.2,9.7,1.4,'Not shown by these data: effects relative to controls (sex confounded); composition-independent regulation;\nsample-level pathway change; association with drug exposure; spatial localization; independent validation',fs=7.8)
save(f,'Graphical_abstract')

# Journal upload: write TIFFs as 600-dpi RGB (matplotlib writes RGBA)
import glob
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
for _t in glob.glob('figures/*.tiff'):
    Image.open(_t).convert('RGB').save(_t,compression='tiff_lzw',dpi=(600,600))
