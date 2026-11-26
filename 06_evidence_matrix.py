# Usage (archive root): python3 scripts/06_evidence_matrix.py
# Every cell answers one stated claim (metadata/evidence_rules.md). Built from the result tables, not typed by hand.
import pandas as pd
D = ['GSE274558 mouse', 'GSE308825 AML903', 'GSE308825 AML001', 'GSE248597 spatial']
cam = pd.read_csv('results/GSE308825/signature_camera_fry.csv').query("analysis=='primary'")
fox = pd.read_csv('results/GSE308825/foxm1_set_tests.csv').query("analysis=='primary'")
mr = pd.read_csv('results/GSE308825/foxm1_mRNA_DE.csv').query("analysis=='primary'")
ag = pd.read_csv('results/GSE308825/sign_agreement.csv'); ora = pd.read_csv('results/GSE308825/overlap_ORA.csv')
mm = pd.read_csv('results/foxm1_GSE274558_set_tests.csv').set_index('set'); md = pd.read_csv('results/foxm1_GSE274558_mRNA_DE.csv').iloc[0]
ms = pd.read_csv('results/foxm1_GSE274558_summary.csv').set_index('measure')
rows = []
def add(claim, dataset, status, detail): rows.append(dict(claim=claim, dataset=dataset, status=status, detail=detail))
def st(ok): return 'supportive' if ok else 'non-significant'
for t, lab in [('Quiz', 'quizartinib'), ('Palb', 'palbociclib'), ('Combo', 'combination')]:
    c = f'Fusion-signature gene set shifts with {lab} exposure (CAMERA)'
    add(c, D[0], 'unavailable', 'no drug exposure'); add(c, D[3], 'unavailable', 'no drug exposure')
    for m, d in [('AML903', D[1]), ('AML001', D[2])]:
        x = cam[cam.contrast == f'{m}_{t}']
        if x.empty: add(c, d, 'unavailable', 'one sample' if t == 'Palb' else 'arm absent')
        else: add(c, d, st((x.camera_FDR < .05).any()), f'min FDR {x.camera_FDR.min():.2f}')
    c = f'Gene-level sign agreement with {lab} exceeds or falls below matched random (added; ignores gene correlation)'
    add(c, D[0], 'unavailable', 'no drug exposure'); add(c, D[3], 'unavailable', 'no drug exposure')
    for m, d in [('AML903', D[1]), ('AML001', D[2])]:
        x = ag[ag.contrast == f'{m}_{t}']
        if x.empty: add(c, d, 'unavailable', 'one sample' if t == 'Palb' else 'arm absent'); continue
        n = int((x.empirical_FDR < .05).sum())
        txt = '; '.join(f'{r.signature} {r.same_sign_fraction:.2f} (random {r.random_median_fraction:.2f})' for r in x.itertuples())
        add(c, d, 'supportive' if n == len(x) else 'non-significant', ('only one of two signatures passes: ' if 0 < n < len(x) else '') + txt)
c = 'Signature genes overlap drug-associated gene list beyond chance'
add(c, D[0], 'unavailable', 'no drug exposure'); add(c, D[3], 'unavailable', 'no drug exposure')
for m, d in [('AML903', D[1]), ('AML001', D[2])]:
    x = ora[ora.contrast.str.startswith(m)]; assert (x.drug_genes == 0).all()
    add(c, d, 'unavailable', 'no gene at FDR<0.05')
c = 'FOXM1 mRNA differs between the compared conditions'
add(c, D[0], st(md.edgeR_FDR < .05 and md.DESeq2_FDR < .05), f'log2FC {md.edgeR_logFC:.2f}, FDR {md.edgeR_FDR:.2f}')
for m, d in [('AML903', D[1]), ('AML001', D[2])]:
    x = mr[mr.contrast.str.startswith(m)]
    add(c, d, st((x.FDR_genomewide < .05).any()), '; '.join(f'{r.contrast.split("_")[1]} {r.logFC:.2f} (FDR {r.FDR_genomewide:.2f})' for r in x.itertuples()))
add(c, D[3], 'unavailable', 'no biopsy IDs')
for s, c in [('targets', 'FOXM1 target set differs between conditions'), ('proliferation', 'Proliferation marker set differs between conditions')]:
    if s == 'targets':
        t = mm.loc['targets']; add(c, D[0], st(t.camera_p < .05 and t.fry_p < .05), f'higher in NM; CAMERA p {t.camera_p:.3f}; selected post hoc, not independent')
        key = 'FOXM1_targets_signed'
    else:
        p = ms.loc['proliferation_set', 'welch_p']; add(c, D[0], st(p < .05), f'Welch p {p:.2f}'); key = 'proliferation_set'
    for m, d in [('AML903', D[1]), ('AML001', D[2])]:
        x = fox[(fox.set == key) & fox.contrast.str.startswith(m)]
        add(c, d, st((x.camera_FDR < .05).any()), f'min CAMERA FDR {x.camera_FDR.min():.2f}')
    add(c, D[3], 'unavailable', 'no biopsy IDs')
for c, why in [('Fusion signature or FOXM1 activity differs by spatial context', ['no spatial data'] * 3 + ['no biopsy IDs; descriptive only']),
               ('Immune organisation differs between region types (tested)', ['no spatial data'] * 3 + ['no biopsy IDs; descriptive only']),
               ('Difference is independent of cell composition', ['bulk, composition confounded', 'bulk', 'bulk', 'no biopsy IDs']),
               ('Replicated in human AML with the same NPM1 fusions', ['mouse model', 'genotype not established', 'genotype not established', 'genotype not given']),
               ('Predicts clinical therapeutic response', ['no treatment', 'drug exposure in xenograft', 'drug exposure in xenograft', 'treatment-naive']),
               ('Mechanism (causal role) demonstrated', ['no perturbation of FOXM1', 'no perturbation of FOXM1', 'no perturbation of FOXM1', 'observational'])]:
    for d, w in zip(D, why): add(c, d, 'unavailable', w)
e = pd.DataFrame(rows); assert len(e) == e.claim.nunique() * 4 and not e.duplicated(['claim', 'dataset']).any()
e.to_csv('results/evidence_matrix_long.csv', index=False)
e.assign(cell=e.status + ': ' + e.detail).pivot(index='claim', columns='dataset', values='cell').reindex(index=e.claim.unique(), columns=D).to_csv('results/evidence_matrix.csv')
print(e.status.value_counts()); print(e[e.status != 'unavailable'].to_string())
