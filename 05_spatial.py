# Usage (archive root): python3 -I scripts/05_spatial.py
# GSE248597 GeoMx: DESCRIPTIVE ONLY. Regions are nested in biopsies and no biopsy identifier is available,
# so no between-region-type tests are made. The random-set percentile refers to gene sampling, not to patients.
import pandas as pd, numpy as np, openpyxl, os
rng = np.random.default_rng(20261126); os.makedirs('results/GSE248597', exist_ok=True)
rows = list(openpyxl.load_workbook('inputs/external/GSE248597_Q3.xlsx', read_only=True)['Sheet1'].iter_rows(values_only=True))
reg = pd.read_csv('metadata/GSE248597_regions.csv')
assert [h.replace('LYMPHOID', '') for h in rows[0][1:]] == reg.region_type.tolist()
X = pd.DataFrame([r[1:] for r in rows[1:]], index=[r[0] for r in rows[1:]], columns=reg.gsm).astype(float)
assert X.index.is_unique and not X.isna().any().any()
L = np.log2(X + 1)
R = X.rank(axis=0, method='average') / len(X)          # within-region percentile rank
immune = {'T_cell': 'CD3D CD3E CD3G TRAC', 'CD8_T': 'CD8A CD8B', 'B_cell': 'CD19 MS4A1 CD79A CD79B', 'Plasma_Ig': 'JCHAIN MZB1 IGHG1 IGKC',
          'NK_cytotoxic': 'NCAM1 NKG7 KLRD1 GNLY', 'Monocyte_macrophage': 'CD68 CD14 LYZ CD163', 'Myeloid_progenitor_blast': 'CD34 MPO KIT',
          'Erythroid': 'GYPA HBB ALAS2', 'Lymphoid_aggregate_chemokines': 'CCL19 CCL21 CXCL13 LTB'}
sets = {k: [(g, 1) for g in v.split()] for k, v in immune.items()}
for sn, f in [('S398', 'fusion_398_human'), ('S336', 'marker_excluded_human')]:
    s = pd.read_csv(f'signatures/{f}.csv')
    for d in ['higher_in_NM', 'higher_in_NC']:
        sets[f'{sn}_{d}'] = [(g, 1) for g in s[s.direction == d].human_symbol]
fx = pd.read_csv('signatures/foxm1_targets_human.csv').dropna(subset=['human_symbol'])
prolif = 'MKI67 TOP2A PCNA MCM2 CCNB1 CDK1 BUB1 AURKB'.split()
sets['FOXM1_targets_signed'] = list(zip(fx.human_symbol, fx.sign))
sets['FOXM1_targets_without_proliferation_set'] = [(g, w) for g, w in sets['FOXM1_targets_signed'] if g not in prolif]
sets['proliferation_set'] = [(g, 1) for g in prolif]
mean_expr = L.mean(axis=1); bins = pd.qcut(mean_expr.rank(method='first'), 20, labels=False).values
pools = [np.flatnonzero(bins == b) for b in range(20)]
types = ['LYMPHOIDAGGREGATE', 'MIXED1', 'MIXED2', 'CONTROL']
lab = reg.region_type.replace('AGGREGATE', 'LYMPHOIDAGGREGATE').values
Rv = R.values; typemean = np.stack([Rv[:, lab == t].mean(axis=1) for t in types])     # type x gene
scores, cov, ref = {}, [], []
for name, items in sets.items():
    got = [(g, w) for g, w in items if g in X.index]
    cov.append(dict(set=name, genes_in_set=len(items), genes_measured=len(got), measured=';'.join(g for g, _ in got) if len(items) < 60 else '',
                    absent=';'.join(g for g, _ in items if g not in X.index) if len(items) < 60 else ''))
    j = X.index.get_indexer([g for g, _ in got]); w = np.array([w for _, w in got])
    scores[name] = np.where(w[:, None] > 0, Rv[j], 1 - Rv[j]).mean(axis=0)
    null = np.empty((1000, len(types)))
    for it in range(1000):
        k = np.empty(len(j), int)
        for b in np.unique(bins[j]):
            m = bins[j] == b; k[m] = rng.choice(pools[b], m.sum(), replace=False)
        null[it] = np.where(w[None, :] > 0, typemean[:, k], 1 - typemean[:, k]).mean(axis=1)
    for i, t in enumerate(types):
        obs = scores[name][lab == t].mean()
        ref.append(dict(set=name, region_type=t, genes=len(got), mean_score=obs, random_median=np.median(null[:, i]),
                        random_lo95=np.quantile(null[:, i], .025), random_hi95=np.quantile(null[:, i], .975), percentile_in_random=100 * (null[:, i] < obs).mean()))
S = pd.DataFrame(scores, index=reg.gsm)
S.insert(0, 'FOXM1_log2_Q3', L.loc['FOXM1'].values)
full = reg.set_index('gsm').join(S); full.to_csv('results/GSE248597/region_scores.csv')
pd.DataFrame(cov).to_csv('results/GSE248597/set_coverage.csv', index=False)
pd.DataFrame(ref).to_csv('results/GSE248597/random_set_reference.csv', index=False)
long = full.melt(id_vars=['region_type', 'slide_barcode'], value_vars=list(S.columns), var_name='measure')
def desc(g): return g.value.agg(regions='size', median='median', mean='mean', q25=lambda v: v.quantile(.25), q75=lambda v: v.quantile(.75)).reset_index()
desc(long.groupby(['measure', 'region_type'])).to_csv('results/GSE248597/descriptive_by_region_type.csv', index=False)
desc(long.groupby(['measure', 'slide_barcode', 'region_type'])).to_csv('results/GSE248597/descriptive_by_slide_and_type.csv', index=False)
print(pd.DataFrame(cov)[['set', 'genes_in_set', 'genes_measured', 'absent']].to_string())
print(full.groupby('region_type')[list(S.columns)].mean().T.round(3).to_string())
print(pd.crosstab(reg.slide_barcode, reg.region_type))
r = pd.DataFrame(ref); print(r[~r.set.isin(immune)][['set', 'region_type', 'mean_score', 'random_median', 'percentile_in_random']].round(3).to_string())
w = full.groupby(['slide_barcode', 'region_type'])[list(S.columns)].mean().T.round(3); print(w.to_string())
