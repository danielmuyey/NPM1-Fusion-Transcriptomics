"""Usage (archive root): python3 scripts/01_freeze_signatures.py
Freeze the direct-comparison signatures. Direction: positive log2FC = higher in NPM1::MLF1 (NM) than NPM1::CCDC28A (NC)."""
import pandas as pd,hashlib
t=pd.read_csv('inputs/NM_vs_NC_all_genes_annotated.csv')
cols=['gene','mouse_symbol','human_id','human_symbol','edgeR_logFC','edgeR_FDR','DESeq2_log2FC','DESeq2_FDR','concordant_lfc1']
c=t[t.concordant].copy(); c['direction']=c.edgeR_logFC.map(lambda x:'higher_in_NM' if x>0 else 'higher_in_NC')
c=c.sort_values(['direction','edgeR_FDR','gene'])[cols+['direction']].rename(columns={'gene':'mouse_ensembl_id'})
h=c.dropna(subset=['human_id'])
assert len(c)==522 and len(h)==398
# exclusion documentation
A=pd.read_csv('inputs/derived_lineage_markers.csv'); B=pd.read_csv('inputs/mMCPcounter_haematopoietic_markers.csv')
x=pd.read_csv('inputs/excluded_lineage_markers.csv')
reason={}
primary=A[((A.source=='ImmGen')|((A.source=='MouseRNAseq')&(A.lineage=='Erythroid')))&A.gene.notna()]
for r in primary.itertuples(): reason.setdefault(r.gene,[]).append(f'{r.source}:{r.lineage}')
for r in B.dropna(subset=['gene']).itertuples(): reason.setdefault(r.gene,[]).append(f'mMCPcounter:{r.lineage}')
assert set(reason)==set(x.gene)
ex=c[c.mouse_ensembl_id.isin(x.gene)].copy()
ex['excluded_because']=ex.mouse_ensembl_id.map(lambda g:'; '.join(sorted(set(reason[g]))))
cm=c[~c.mouse_ensembl_id.isin(x.gene)]; hm=cm.dropna(subset=['human_id'])
out={'fusion_522_mouse':c,'fusion_398_human':h,'marker_excluded_mouse':cm,'marker_excluded_human':hm,'excluded_genes':ex}
for n,d in out.items(): d.to_csv(f'signatures/{n}.csv',index=False)
all_=x.copy(); all_['excluded_because']=all_.gene.map(lambda g:'; '.join(sorted(set(reason[g])))); all_['in_522']=all_.gene.isin(c.mouse_ensembl_id)
all_.to_csv('signatures/all_excluded_lineage_markers.csv',index=False)
with open('signatures/MANIFEST.tsv','w') as f:
    f.write('file\tgenes\thigher_in_NM\thigher_in_NC\tsha256\n')
    for n,d in out.items():
        b=open(f'signatures/{n}.csv','rb').read()
        f.write(f'{n}.csv\t{len(d)}\t{(d.direction=="higher_in_NM").sum()}\t{(d.direction=="higher_in_NC").sum()}\t{hashlib.sha256(b).hexdigest()}\n')
print(open('signatures/MANIFEST.tsv').read())
# note: post hoc progenitor genes are NOT excluded (not in the predefined set)
post=A[A.lineage=='Myeloid_progenitor_vs_mature']
print('post hoc progenitor genes retained in marker-excluded signature:',sorted(post[post.gene.isin(cm.mouse_ensembl_id)].symbol))
print('FOXM1:',t[t.mouse_symbol=='Foxm1'][['edgeR_logFC','edgeR_FDR','DESeq2_log2FC','DESeq2_FDR','concordant']].to_dict('records'))
