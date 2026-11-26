"""Usage (archive root): python3 scripts/02_sensitivity_enrichment.py  (after 01_lineage.R)
Exploratory lineage-marker exclusion: ORA + preranked GSEA without markers, with a random-removal control for ORA."""
import numpy as np,pandas as pd,gseapy
from scipy.stats import hypergeom,false_discovery_control as bh
t=pd.read_csv('inputs/NM_vs_NC_all_genes_annotated.csv'); bg=t.dropna(subset=['human_id']).reset_index(drop=True)
ex=pd.read_csv('results/marker_enrichment/excluded_lineage_markers.csv')
excl=set(ex.human_id.dropna())&set(bg.human_id)
# also drop any row whose mouse gene is a marker
excl|=set(bg[bg.gene.isin(ex.gene)].human_id)
ids=bg.human_id.values; pos={g:i for i,g in enumerate(ids)}; n=len(ids)
sets={}
for line in open('inputs/msigdb_7.5.1_human_GO_Reactome_ensembl.gmt'):
    f=line.rstrip('\n').split('\t'); s=set(f[2:])&set(ids)
    if 10<=len(s)<=500: sets.setdefault(f[1],{})[f[0]]=s
M={c:(list(S),np.array([[g in S[k] for g in ids] for k in S],dtype=np.int32)) for c,S in sets.items()}
isx=np.isin(ids,list(excl))
def ora(query,keep):
    out={}
    for c,(names,A) in M.items():
        q=(query&keep).astype(np.int32); K=A@keep.astype(np.int32); k=A@q
        p=hypergeom.sf(k-1,keep.sum(),K,q.sum()); p=np.where(K>0,p,1.0)
        out[c]=(names,k,K,int(q.sum()),int(keep.sum()),p,bh(p))
    return out
rng=np.random.default_rng(1); rows=[]; summ=[]
for lst in ['concordant','concordant_lfc1']:
  for d,sign in [('up_in_NM',1),('up_in_NC',-1)]:
    query=(bg[lst].values.astype(bool))&(np.sign(bg.edgeR_logFC.values)==sign)
    full=ora(query,np.ones(n,bool)); red=ora(query,~isx)
    nq=int((query&isx).sum()); nb=int((~query&isx).sum())
    qi=np.flatnonzero(query); oi=np.flatnonzero(~query); rand={c:[] for c in M}
    for it in range(200):
        keep=np.ones(n,bool); keep[rng.choice(qi,nq,replace=False)]=False; keep[rng.choice(oi,nb,replace=False)]=False
        r=ora(query,keep)
        for c in M: rand[c].append(int((r[c][6]<.05).sum()))
    for c in M:
        names,k,K,Q,B,p,f=red[c]; f0=full[c][6]
        for i,nm in enumerate(names):
            rows.append((lst,d,c,nm,int(full[c][1][i]),float(f0[i]),int(k[i]),int(K[i]),Q,B,float(p[i]),float(f[i])))
        a=f0<.05; b=f<.05
        summ.append(dict(list=lst,direction=d,collection=c,query_full=int(query.sum()),query_markers_removed=nq,significant_full=int(a.sum()),
            significant_without_markers=int(b.sum()),retained=int((a&b).sum()),lost=int((a&~b).sum()),gained=int((~a&b).sum()),
            random_removal_median=float(np.median(rand[c])),random_removal_2_5pct=float(np.percentile(rand[c],2.5)),
            random_removal_97_5pct=float(np.percentile(rand[c],97.5))))
R=pd.DataFrame(rows,columns='list direction collection term overlap_full FDR_full overlap_nomarkers set_size_nomarkers query_nomarkers background_nomarkers p_nomarkers FDR_nomarkers'.split())
R.to_csv('results/sensitivity/ORA_full_vs_marker_excluded.csv',index=False)
S=pd.DataFrame(summ); S.to_csv('results/sensitivity/ORA_summary.csv',index=False); print(S.to_string(index=False),flush=True)
bg2=bg[~isx]
rnk=bg2.assign(score=np.sign(bg2.edgeR_logFC)*np.sqrt(bg2.edgeR_F)).sort_values(['score','human_id'],ascending=[False,True])[['human_id','score']]
G=[]
for c,S_ in sets.items():
    gs={k:sorted(v-excl) for k,v in S_.items()}
    res=gseapy.prerank(rnk=rnk,gene_sets=gs,min_size=10,max_size=500,permutation_num=1000,seed=1,threads=2,outdir=None,verbose=False).res2d.copy()
    res.insert(0,'collection',c); G.append(res); print(c,len(res),flush=True)
G=pd.concat(G).drop(columns=['Name'])
for k in ['NES','FDR q-val','NOM p-val','ES']: G[k]=pd.to_numeric(G[k])
G.to_csv('results/sensitivity/GSEA_marker_excluded.csv',index=False)
F=pd.read_csv('inputs/GSEA_all_terms.csv')
m=F[['collection','Term','NES','FDR q-val']].merge(G[['collection','Term','NES','FDR q-val']],on=['collection','Term'],how='left',suffixes=('_full','_nomarkers'))
m.to_csv('results/sensitivity/GSEA_full_vs_marker_excluded.csv',index=False)
out=[]
for (c),g in m.groupby('collection'):
    for d,s in [('NM',1),('NC',-1)]:
        a=(g['FDR q-val_full']<.05)&(np.sign(g.NES_full)==s); b=(g['FDR q-val_nomarkers']<.05)&(np.sign(g.NES_nomarkers)==s)
        out.append(dict(collection=c,higher_in=d,significant_full=int(a.sum()),significant_without_markers=int(b.sum()),retained=int((a&b).sum()),lost=int((a&~b).sum()),gained=int((~a&b).sum())))
out=pd.DataFrame(out); out.to_csv('results/sensitivity/GSEA_summary.csv',index=False); print(out.to_string(index=False))
import scipy
open('results/sessionInfo_python.txt','w').write(f'gseapy {gseapy.__version__}\nnumpy {np.__version__}\npandas {pd.__version__}\nscipy {scipy.__version__}\n')
