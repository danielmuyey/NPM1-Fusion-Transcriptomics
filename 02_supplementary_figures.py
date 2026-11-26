"""Usage: python3 scripts/02_supplementary_figures.py <folder holding the four analysis archives>
S1 and S4 are drawn here. S2, S3 and S5-S11 are the audited plotting scripts rerun on a temporary copy with an Arial-metric font."""
import sys,os,shutil,subprocess,tempfile,re,numpy as np,pandas as pd,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R=os.path.abspath(sys.argv[1])+'/'; D=R+'GSE274558_Direct_Fusion_Comparison/'; Ld=R+'GSE274558_Lineage_and_Cell_Composition_Analysis/'; Xd=R+'AML_Cross_Dataset_Integration_Analysis/'
OUT=os.path.abspath('supplementary_figures'); PRE=os.path.abspath('metadata')
FONT={'font.family':'Liberation Sans','font.size':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False}
plt.rcParams.update(FONT); CM,CC,GR='#eb6834','#2a78d6','#c8c7c2'; S=['NM1','NM2','NM3','NC1','NC2','NC3']; lab=[s.replace('NM','MLF1 ').replace('NC','CCDC28A ') for s in S]
plt.rcParams.update({'mathtext.fontset':'custom','mathtext.rm':'Liberation Sans','mathtext.it':'Liberation Sans:italic','mathtext.bf':'Liberation Sans:bold','mathtext.default':'regular','mathtext.cal':'Liberation Sans'})
def save(f,n):
    f.savefig(f'{OUT}/{n}.pdf',bbox_inches='tight',metadata={'CreationDate':None}); f.savefig(f'{OUT}/{n}.tiff',dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'})
    f.savefig(f'{PRE}/preview_{n}.png',dpi=110,bbox_inches='tight'); plt.close(f)
# S1
cr=pd.read_csv(D+'results/qc/spearman.csv',index_col=0).loc[S,S]; d=pd.read_csv(D+'results/de/NM_vs_NC_all_genes_annotated.csv')
f,A=plt.subplots(1,2,figsize=(7.1,3.1),layout='constrained'); a=A[0]; im=a.imshow(cr.values,cmap='Blues'); a.set_xticks(range(6),lab,rotation=90); a.set_yticks(range(6),lab)
for i in range(6):
    for j in range(6): a.text(j,i,f'{cr.values[i,j]:.3f}',ha='center',va='center',fontsize=5.5,color='white' if cr.values[i,j]>(cr.values.min()+cr.values.max())/2 else 'black')
a.set_title('A  Spearman correlation',loc='left',fontweight='bold'); f.colorbar(im,ax=a,shrink=.8)
a=A[1]; c=d.concordant; a.scatter(d.edgeR_logFC[~c],d.DESeq2_log2FC[~c],s=2,color=GR,rasterized=True,label='Not concordant'); a.scatter(d.edgeR_logFC[c],d.DESeq2_log2FC[c],s=4,color='#333333',rasterized=True,label='Concordant')
a.set(xlabel='edgeR log$_2$ fold change',ylabel='DESeq2 log$_2$ fold change'); a.legend(); a.set_title('B  Agreement between methods',loc='left',fontweight='bold'); save(f,'FigureS1')
# S4 (cited after S3)
o=pd.read_csv(D+'results/enrichment/ORA_all_terms.csv'); o=o[(o['list']=='concordant_FDR')&(o.FDR<.05)]
g=pd.read_csv(Ld+'results/camera_roast/GO_Reactome_camera_fry_all_genes.csv')
f=plt.figure(figsize=(7.1,6.6),layout='constrained'); gs=f.add_gridspec(2,2,height_ratios=[1.5,1])
a=f.add_subplot(gs[0,:]); rows=pd.concat([o[o.direction=='up_in_NM'].nsmallest(10,'FDR'),o[o.direction=='up_in_NC'].nsmallest(3,'FDR')])
a.barh(range(len(rows)),-np.log10(rows.FDR),color=[CM if x=='up_in_NM' else CC for x in rows.direction],height=.7)
a.set_yticks(range(len(rows)),[re.sub(r'^(GOBP|GOCC|GOMF|REACTOME)_',r'\1: ',t).replace('_',' ').lower()[:80] for t in rows.term],fontsize=6); a.invert_yaxis(); a.set_xlabel('-log$_{10}$ FDR (orange, higher in MLF1; blue, higher in CCDC28A)')
a.set_title('A  Over-represented terms (descriptive, gene level)',loc='left',fontweight='bold')
cols=[c for c in g.columns if c.endswith('FDR')]; want=[('camera_FDR','CAMERA'),('fry_FDR','fry')]
for k,(col,t) in enumerate(want):
    a=f.add_subplot(gs[1,k]); assert col in g.columns,(col,cols); a.hist(g[col],bins=np.linspace(0,1,41),color='#8a8984'); a.axvline(.05,color='black',lw=.8,ls='--')
    a.set(xlabel='FDR',ylabel='gene sets'); a.set_title(('B  ' if k==0 else 'C  ')+f'{t}\nminimum FDR {g[col].min():.3f}',loc='left',fontweight='bold',fontsize=8)
assert len(g)==6611; save(f,'FigureS4')
# S2, S3, S5-S11: rerun audited plotting code on copies
def rerun(src,script,mapping):
    with tempfile.TemporaryDirectory() as t:
        w=t+'/w'; shutil.copytree(src,w); code=open(w+'/'+script).read().replace('DejaVu Sans','Liberation Sans')
        code=code.replace("def save(","plt.rcParams.update({'font.family':'Liberation Sans'})\ndef save(",1); open(w+'/'+script,'w').write(code)
        shutil.rmtree(w+'/figures'); os.makedirs(w+'/figures'); subprocess.run([sys.executable,script],cwd=w,check=True,stdout=subprocess.DEVNULL)
        for old,new in mapping.items():
            for e in ('pdf','tiff'): shutil.copy(f'{w}/figures/{old}.{e}',f'{OUT}/{new}.{e}')
            from PIL import Image; Image.MAX_IMAGE_PIXELS=None; im=Image.open(f'{OUT}/{new}.tiff'); im.thumbnail((1100,1400)); im.convert('RGB').save(f'{PRE}/preview_{new}.png')
rerun(Ld,'scripts/05_figures.py',{'Figure1_marker_heatmap':'FigureS2','Figure4_marker_exclusion_sensitivity':'FigureS3','Figure6_TF_resolved':'FigureS5'})
rerun(Xd,'scripts/07_figures.py',{'Fig1_GSE308825_QC':'FigureS6','Fig2_signature_drug_exposure':'FigureS7','Fig3_FOXM1_across_datasets':'FigureS8','Fig4_spatial_marker_sets':'FigureS9','Fig5_spatial_signature_scores':'FigureS10','Fig6_evidence_matrix':'FigureS11'})

# Journal upload: write TIFFs as 600-dpi RGB (matplotlib writes RGBA)
import glob
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
for _t in glob.glob('supplementary_figures/*.tiff'):
    Image.open(_t).convert('RGB').save(_t,compression='tiff_lzw',dpi=(600,600))
