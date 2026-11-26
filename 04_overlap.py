# Usage (archive root): python3 scripts/04_overlap.py
# Signature x drug-exposure overlap. Planned ORA, plus threshold-free sign agreement. Overlap is not validation.
import pandas as pd, numpy as np
from scipy.stats import hypergeom, spearmanr
rng = np.random.default_rng(20261126); N = 10000
def bh(p):
    p = np.asarray(p, float); q = np.full(p.shape, np.nan); i = np.flatnonzero(~np.isnan(p)); o = i[np.argsort(p[i])]
    q[o] = np.minimum(1, np.minimum.accumulate((p[o] * len(o) / np.arange(1, len(o) + 1))[::-1])[::-1]); return q
sigs = {'S398': pd.read_csv('signatures/fusion_398_human.csv'), 'S336': pd.read_csv('signatures/marker_excluded_human.csv')}
a = pd.read_csv('inputs/NM_vs_NC_all_genes_annotated.csv')
# human genes tested in GSE274558: one-to-one confidence-1 orthologs
o = pd.read_csv('inputs/Mouse_Human_Orthologs_Ensembl116.tsv', sep='\t', quoting=3); o.columns = ['mid','msym','hid','hsym','type','conf']
one = o[(o.type == 'ortholog_one2one') & (o.conf == 1)].drop_duplicates('mid')
bg = set(one[one.mid.isin(a.gene)].hsym.dropna())
contrasts = ['AML903_Quiz', 'AML001_Quiz', 'AML001_Palb', 'AML001_Combo']
de = {c: pd.read_csv(f'results/GSE308825/DE_primary_{c}.csv.gz').set_index('gene') for c in contrasts}
ora, agree, uni = [], [], []
for c, d in de.items():
    u = d[d.index.isin(bg)].copy(); u['bin'] = pd.qcut(u.AveExpr.rank(method='first'), 20, labels=False)
    uni.append(dict(contrast=c, GSE274558_one2one_background=len(bg), eligible_universe=len(u)))
    nbin = np.bincount(u['bin'], minlength=20)
    pos = np.bincount(u['bin'], weights=(u.logFC > 0), minlength=20).astype(int)
    for sn, s in sigs.items():
        s = s[s.human_symbol.isin(u.index)]
        for rule, fc in [('FDR<0.05 and |log2FC|>=1', 1), ('FDR<0.05', 0)]:
            for sd in ['higher_in_NM', 'higher_in_NC']:
                g = s[s.direction == sd].human_symbol
                for dd, sel in [('up', (u['adj.P.Val'] < .05) & (u.logFC > fc)), ('down', (u['adj.P.Val'] < .05) & (u.logFC < -fc))]:
                    K, n, k = int(sel.sum()), len(g), int(sel.reindex(g).sum())
                    r = dict(contrast=c, signature=sn, signature_direction=sd, drug_direction=dd, rule=rule, universe=len(u),
                             signature_genes=n, drug_genes=K, overlap=k, expected=n * K / len(u), status='tested' if K else 'not testable: no drug genes',
                             hypergeometric_p=np.nan, random_matched_p=np.nan)
                    if K:
                        r['hypergeometric_p'] = hypergeom.sf(k - 1, len(u), K, n)
                        hit = np.bincount(u['bin'], weights=sel, minlength=20).astype(int); take = np.bincount(u.loc[g, 'bin'], minlength=20)
                        null = sum(rng.hypergeometric(hit[b], nbin[b] - hit[b], take[b], N) for b in range(20) if take[b])
                        r['random_matched_p'] = (1 + (null >= k).sum()) / (N + 1)
                    ora.append(r)
        # sign agreement: drug logFC sign equals NM-NC sign; null = expression-matched genes with same direction mix
        lfc = d.loc[s.human_symbol, 'logFC'].values; want = np.where(s.direction == 'higher_in_NM', 1, -1)
        k = int((np.sign(lfc) == want).sum()); null = np.zeros(N, int)
        for sign, sd in [(1, 'higher_in_NM'), (-1, 'higher_in_NC')]:
            take = np.bincount(u.loc[s[s.direction == sd].human_symbol, 'bin'], minlength=20)
            good = pos if sign == 1 else np.bincount(u['bin'], weights=(u.logFC < 0), minlength=20).astype(int)
            for b in range(20):
                if take[b]: null += rng.hypergeometric(good[b], nbin[b] - good[b], take[b], N)
        p = min(1, 2 * (1 + min((null >= k).sum(), (null <= k).sum())) / (N + 1))
        agree.append({'contrast': c, 'signature': sn, 'genes': len(s), 'same_sign': k, 'same_sign_fraction': k / len(s),
                      'random_median_fraction': np.median(null) / len(s), 'random_lo95': np.quantile(null, .025) / len(s), 'random_hi95': np.quantile(null, .975) / len(s),
                      'empirical_p_two_sided': p, 'spearman_rho_with_NM_vs_NC': spearmanr(s.edgeR_logFC, lfc)[0]})
ora = pd.DataFrame(ora)
for rule in ora.rule.unique():
    i = ora.rule == rule; ora.loc[i, 'hypergeometric_FDR'] = bh(ora.loc[i, 'hypergeometric_p']); ora.loc[i, 'random_matched_FDR'] = bh(ora.loc[i, 'random_matched_p'])
agree = pd.DataFrame(agree); agree['empirical_FDR'] = bh(agree.empirical_p_two_sided)
# same drug, two models
cross = []
for sn, s in sigs.items():
    g = [x for x in s.human_symbol if x in de['AML903_Quiz'].index and x in de['AML001_Quiz'].index]
    x, y = de['AML903_Quiz'].loc[g, 'logFC'], de['AML001_Quiz'].loc[g, 'logFC']
    shared = de['AML903_Quiz'].index.intersection(de['AML001_Quiz'].index)
    X, Y = de['AML903_Quiz'].loc[shared, 'logFC'], de['AML001_Quiz'].loc[shared, 'logFC']
    cross.append(dict(signature=sn, genes_in_both=len(g), same_sign_fraction=float((np.sign(x) == np.sign(y)).mean()), spearman_rho=spearmanr(x, y)[0],
                      all_genes_same_sign_fraction=float((np.sign(X) == np.sign(Y)).mean()), all_genes_spearman_rho=spearmanr(X, Y)[0]))
ora.to_csv('results/GSE308825/overlap_ORA.csv', index=False); agree.to_csv('results/GSE308825/sign_agreement.csv', index=False)
pd.DataFrame(uni).to_csv('results/GSE308825/overlap_universe.csv', index=False)
pd.DataFrame(cross).to_csv('results/GSE308825/quizartinib_between_model_agreement.csv', index=False)
print(pd.DataFrame(uni)); print(ora.status.value_counts()); print(agree.round(3).to_string()); print(pd.DataFrame(cross).round(3).to_string())

# POST HOC diagnostic: does the AML001 sign agreement hold within each replicate (library set A = R1, B = R2)?
x = pd.read_csv('results/GSE308825/logCPM_all_samples.csv').set_index('gene'); rep = []
for t in ['Quiz', 'Palb', 'Combo']:
    for r in ['R1', 'R2']:
        for sn, s in sigs.items():
            s2 = s[s.human_symbol.isin(de[f'AML001_{t}'].index)]
            diff = x.loc[s2.human_symbol, f'AML001_{t}_{r}'].values - x.loc[s2.human_symbol, f'AML001_DMSO_{r}'].values
            want = np.where(s2.direction == 'higher_in_NM', 1, -1)
            allg = x.loc[de[f'AML001_{t}'].index, f'AML001_{t}_{r}'] - x.loc[de[f'AML001_{t}'].index, f'AML001_DMSO_{r}']
            up = float((allg > 0).mean()); expected = float(np.where(want == 1, up, 1 - up).mean())
            rep.append(dict(treatment=t, replicate=r, signature=sn, genes=len(s2), same_sign_fraction=float((np.sign(diff) == want).mean()),
                            expected_from_all_genes_unmatched=expected))
pd.DataFrame(rep).to_csv('results/GSE308825/posthoc_sign_agreement_by_replicate.csv', index=False)
print(pd.DataFrame(rep).round(3).to_string())
for f in ['replicate_blocked']:
    for c in contrasts[1:]:
        d = pd.read_csv(f'results/GSE308825/DE_{f}_{c}.csv.gz').set_index('gene')
        for sn, s in sigs.items():
            s2 = s[s.human_symbol.isin(d.index)]; w = np.where(s2.direction == 'higher_in_NM', 1, -1)
            print(f, c, sn, round(float((np.sign(d.loc[s2.human_symbol, 'logFC']) == w).mean()), 3))
