"""Usage (output folder): python3 scripts/03_supplementary_tables.py <folder holding the four analysis archives>
Writes Supplementary Tables S1-S13 (CSV + one workbook) with a column dictionary."""
import sys,pandas as pd,openpyxl
from openpyxl.styles import Font
R=sys.argv[1]+'/'; A=R+'GSE274558_Sex_Confounding_Audit/'; D=R+'GSE274558_Direct_Fusion_Comparison/results/'; L=R+'GSE274558_Lineage_and_Cell_Composition_Analysis/results/'; X=R+'AML_Cross_Dataset_Integration_Analysis/'
r=pd.read_csv
de=r(D+'de/NM_vs_NC_all_genes_annotated.csv')
cols=['gene','mouse_symbol','human_id','human_symbol']+[c for c in de.columns if c not in('gene','mouse_symbol','human_id','human_symbol')]
de=de[cols].rename(columns={'gene':'mouse_ensembl_id'})
T=[('S1','Sample sex inference (9 libraries)',[('',r(A+'mouse/results_marker_counts/sample_sex_calls.csv'))]),
('S2','Differential expression, all tested genes',[('',de)]),
('S3','Concordant genes (522)',[('',de[de.concordant])]),
('S4','Human ortholog gene sets (398 and 336)',[('398',r(X+'signatures/fusion_398_human.csv')),('336',r(X+'signatures/marker_excluded_human.csv'))]),
('S5','Filter sensitivity',[('counts',r(D+'de/summary_counts.csv')),('overlap',r(D+'de/filter_sensitivity_overlap.csv'))]),
('S6','Hox, XPO1-related and other preselected genes',[('themes',r(D+'enrichment/theme_genes.csv')),('allHox',de[de.mouse_symbol.fillna('').str.startswith('Hox')])]),
('S7','Lineage-associated expression',[('scores',r(L+'scores/lineage_scores_by_sample.csv')),('tests',r(L+'camera_roast/lineage_sets_camera_fry.csv')),('enrich',r(L+'marker_enrichment/DE_lineage_marker_enrichment.csv')),('ORAsens',r(L+'sensitivity/ORA_summary.csv')),('GSEAsens',r(L+'sensitivity/GSEA_summary.csv'))]),
('S8','Enrichment, descriptive (gene level)',[('ORA',r(D+'enrichment/ORA_significant.csv')),('GSEA',r(D+'enrichment/GSEA_all_terms.csv'))]),
('S9','Correlation-aware gene-set tests, sample level',[('',r(L+'camera_roast/GO_Reactome_camera_fry_all_genes.csv'))]),
('S10','Regulator target scores',[('min10',r(L+'tf/TF_resolved_min10.csv')),('min20',r(L+'tf/TF_resolved_min20.csv'))]),
('S11','FOXM1, mouse data',[('mRNA',r(X+'results/foxm1_GSE274558_mRNA_DE.csv')),('samples',r(X+'results/foxm1_GSE274558_sample_values.csv')),('settests',r(X+'results/foxm1_GSE274558_set_tests.csv')),('targets',r(X+'results/foxm1_GSE274558_targets.csv'))]),
('S12','External datasets',[('drugDE',r(X+'results/GSE308825/DE_summary.csv')),('QC',r(X+'results/GSE308825/sample_qc.csv')),('sets',r(X+'results/GSE308825/signature_camera_fry.csv')),('sign',r(X+'results/GSE308825/sign_agreement.csv')),('FOXM1',r(X+'results/GSE308825/foxm1_set_tests.csv')),('FOXM1mRNA',r(X+'results/GSE308825/foxm1_mRNA_DE.csv')),('spatial',r(X+'results/GSE248597/descriptive_by_region_type.csv')),('coverage',r(X+'results/GSE248597/set_coverage.csv'))]),
('S13','Evidence matrix',[('',r(X+'results/evidence_matrix.csv'))])]
DICT=[('Direction','In GSE274558 tables a positive effect means higher with NPM1::MLF1 (NM); negative means higher with NPM1::CCDC28A (NC).'),
('NM / NC / WT','NPM1::MLF1 / NPM1::CCDC28A / empty-vector control sample labels.'),
('edgeR_logFC, DESeq2_log2FC','log2 fold change, NM vs NC.'),('*_p / p','Unadjusted p value.'),('*_FDR / FDR / padj / adj.P.Val','Benjamini-Hochberg adjusted value within the analysis named by the table (not across the study).'),
('concordant / concordant_lfc1','FDR<0.05 in both methods, same direction / and |log2 FC|>=1 in both.'),
('S8, gene-level tables','Gene-level over-representation and gene-permutation GSEA: descriptive only; not sample-level evidence.'),
('S9, camera_preset columns','Preset inter-gene correlation 0.01; shown for transparency, not used for inference.'),
('S7','Scores are expression summaries, not cell proportions. "[post hoc]" sets were defined after the first run. Rows flagged for fewer than 10 genes are unreliable.'),
('S10, renamed','Regulator name translated from a UniProt accession or non-mouse symbol. Activity is inferred, not measured. Regulon file rebuilt from run.R logic; the original source table and its download date were not retained.'),
('S11','FOXM1 was selected from the same samples; not independent evidence.'),('S12','Drug exposure is not clinical response. Spatial values are descriptive; regions are not independent. S12_sign: gene-level sign agreement, added during analysis; ignores gene correlation.'),
('Gene-set files','MSigDB 7.5.1 membership files are not redistributed here; term names and results only.')]
wb=openpyxl.Workbook(); w=wb.active; w.title='README'; w.append(['Item','Meaning']); [w.append(list(x)) for x in DICT]; w.append([]); w.append(['Sheet','Content','Rows','Columns'])
for t,title,parts in T:
    for suf,df in parts:
        name=t+('_'+suf if suf else ''); df.to_csv(f'supplementary_tables/{name}.csv',index=False)
        w.append([name,title,len(df),'; '.join(map(str,df.columns))])
        s=wb.create_sheet(name[:31]); s.append([str(c) for c in df.columns])
        for row in df.itertuples(index=False): s.append([None if (isinstance(v,float) and v!=v) else (v.item() if hasattr(v,'item') else v) for v in row])
        for c in s[1]: c.font=Font(bold=True)
        s.freeze_panes='A2'
for c in w[1]: c.font=Font(bold=True)
w.column_dimensions['A'].width=28; w.column_dimensions['B'].width=90
wb.save('Supplementary_Tables.xlsx'); print(len(wb.sheetnames)-1,'data sheets')
