"""Usage (from archive root): python3 scripts/02_enrichment.py
ORA (hypergeometric, BH within collection x direction x list) and preranked GSEA (gseapy)."""
import numpy as np,pandas as pd,gseapy
from scipy.stats import hypergeom,false_discovery_control as bh
t=pd.read_csv('results/de/NM_vs_NC_all_genes_annotated.csv')
bg=t.dropna(subset=['human_id']); B=set(bg.human_id); sym=dict(zip(bg.human_id,bg.human_symbol))
sets={}
for line in open('inputs/msigdb_7.5.1_human_GO_Reactome_ensembl.gmt'):
    f=line.rstrip('\n').split('\t'); s=set(f[2:])&B
    if 10<=len(s)<=500: sets.setdefault(f[1],{})[f[0]]=s
rows=[]
for lst,col in [('concordant_FDR','concordant'),('concordant_FDR_lfc1','concordant_lfc1')]:
    for d,sign in [('up_in_NM',1),('up_in_NC',-1)]:
        q=set(bg[bg[col]&(np.sign(bg.edgeR_logFC)==sign)].human_id)
        for coll,S in sets.items():
            r=[]
            for name,g in S.items():
                k=len(q&g)
                r.append(dict(list=lst,direction=d,collection=coll,term=name,overlap=k,query_size=len(q),set_size=len(g),
                    background=len(B),fold=(k/len(q))/(len(g)/len(B)) if q else 0,p=hypergeom.sf(k-1,len(B),len(g),len(q)),
                    genes=';'.join(sorted(sym[i] for i in q&g))))
            r=pd.DataFrame(r); r['FDR']=bh(r.p); rows.append(r)
ora=pd.concat(rows); ora.to_csv('results/enrichment/ORA_all_terms.csv',index=False)
ora[ora.FDR<.05].sort_values(['list','direction','collection','FDR']).to_csv('results/enrichment/ORA_significant.csv',index=False)
c=ora.assign(sig=ora.FDR<.05).groupby(['list','direction','collection']).agg(terms_tested=('term','size'),significant=('sig','sum'),query_size=('query_size','first'))
c.to_csv('results/enrichment/ORA_counts.csv'); print(c)
rnk=bg.assign(score=np.sign(bg.edgeR_logFC)*np.sqrt(bg.edgeR_F)).sort_values(['score','human_id'],ascending=[False,True])[['human_id','score']]
rnk.to_csv('results/enrichment/GSEA_ranked_list.tsv',sep='\t',index=False,header=False)
print('tied scores:',int(rnk.score.duplicated().sum()))
G=[]
for coll,S in sets.items():
    res=gseapy.prerank(rnk=rnk,gene_sets={k:sorted(v) for k,v in S.items()},min_size=10,max_size=500,
        permutation_num=1000,seed=1,threads=2,outdir=None,verbose=False).res2d.copy()
    res.insert(0,'collection',coll); G.append(res); print(coll,len(res),flush=True)
G=pd.concat(G).drop(columns=['Name']); G['Lead_genes']=G.Lead_genes.map(lambda x:';'.join(sym.get(i,i) for i in str(x).split(';')))
for k in ['NES','FDR q-val','NOM p-val','ES']: G[k]=pd.to_numeric(G[k])
G.to_csv('results/enrichment/GSEA_all_terms.csv',index=False)
g=G.assign(NM=(G['FDR q-val']<.05)&(G.NES>0),NC=(G['FDR q-val']<.05)&(G.NES<0)).groupby('collection').agg(tested=('Term','size'),sig_up_NM=('NM','sum'),sig_up_NC=('NC','sum'))
g.to_csv('results/enrichment/GSEA_counts.csv'); print(g)

# Themes: lists fixed before any gene or term names from this comparison were inspected
themes={'HOX/MEIS':['Hoxa5','Hoxa7','Hoxa9','Hoxa10','Hoxb2','Hoxb4','Hoxb5','Meis1','Pbx3'],
 'XPO1-associated':['Xpo1','Nup98','Nup214','Ran','Rangap1','Ranbp1'],
 'Myeloid differentiation':['Cebpa','Cebpe','Spi1','Mpo','Elane','Itgam','Cd14','Ly6g','Kit','Cd34','Csf1r'],
 'Metabolism':['Hk2','Pfkp','Ldha','Idh1','Idh2','Cpt1a','Fasn','Mdh2'],
 'Immune regulation':['H2-Ab1','Cd74','Ciita','Cd274','Tlr2','Tlr4','Il1b','Tnf','Stat1','Irf8']}
rows=[]
for th,genes in themes.items():
    for g_ in genes:
        r=t[t.mouse_symbol==g_]
        rows.append(dict(theme=th,mouse_symbol=g_,tested=not r.empty,**({} if r.empty else r.iloc[0][['edgeR_logFC','edgeR_FDR','DESeq2_log2FC','DESeq2_FDR','concordant','concordant_lfc1']].to_dict())))
pd.DataFrame(rows).to_csv('results/enrichment/theme_genes.csv',index=False)
terms=['GOBP_MYELOID_CELL_DIFFERENTIATION','GOBP_PROTEIN_EXPORT_FROM_NUCLEUS','REACTOME_ACTIVATION_OF_ANTERIOR_HOX_GENES_IN_HINDBRAIN_DEVELOPMENT_DURING_EARLY_EMBRYOGENESIS',
 'GOBP_OXIDATIVE_PHOSPHORYLATION','GOBP_GLYCOLYTIC_PROCESS','GOBP_FATTY_ACID_BETA_OXIDATION','GOBP_INNATE_IMMUNE_RESPONSE','GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION','REACTOME_INTERFERON_SIGNALING']
th=G[G.Term.isin(terms)][['collection','Term','NES','NOM p-val','FDR q-val']]; th.to_csv('results/enrichment/theme_terms_GSEA.csv',index=False)
print('theme terms missing:',set(terms)-set(th.Term)); print(th.to_string(index=False))
open('results/sessionInfo_python.txt','w').write(f'gseapy {gseapy.__version__}\nnumpy {np.__version__}\npandas {pd.__version__}\n')
