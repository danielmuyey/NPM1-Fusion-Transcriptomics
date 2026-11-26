# Usage (archive root): python3 scripts/07_figures.py   -> figures/*.pdf and *.tiff (600 dpi, LZW)
import pandas as pd, numpy as np, matplotlib, os
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
os.makedirs('figures', exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.titlesize': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'axes.edgecolor': '#777777', 'axes.linewidth': .6, 'legend.frameon': False})
BLUE, ORANGE, VIOLET, GREY = '#2a78d6', '#eb6834', '#4a3aa7', '#8a8984'
TC = {'DMSO': GREY, 'Quiz': BLUE, 'Palb': ORANGE, 'Combo': VIOLET}; MK = {'DMSO': 'o', 'Quiz': 's', 'Palb': '^', 'Combo': 'D'}
def save(f, name):
    f.savefig(f'figures/{name}.pdf', bbox_inches='tight', metadata={'CreationDate': None})
    f.savefig(f'figures/{name}.tiff', dpi=600, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'}); plt.close(f)
handles = [Line2D([], [], marker=MK[t], ls='', color=TC[t], label=t) for t in TC]
# 1 QC
q = pd.read_csv('results/GSE308825/sample_qc.csv'); q = q.sort_values(['model', 'treatment', 'replicate']); pc = pd.read_csv('results/GSE308825/pca.csv')
f, ax = plt.subplots(2, 2, figsize=(7.2, 6), layout='constrained')
for a, col, lab, div in [(ax[0, 0], 'total_counts', 'assigned counts (millions)', 1e6), (ax[0, 1], 'genes_detected', 'genes detected (thousands)', 1e3)]:
    a.barh(range(len(q)), q[col] / div, color=[TC[t] for t in q.treatment], height=.7); a.set_yticks(range(len(q)), q.title); a.invert_yaxis(); a.set_xlabel(lab)
for a, m in zip(ax[1], ['AML903', 'AML001']):
    x = pc[pc.model == m]
    for r in x.itertuples():
        _, t, rp = r.sample.split('_'); a.scatter(r.PC1, r.PC2, color=TC[t], marker=MK[t], s=30); a.annotate(rp, (r.PC1, r.PC2), fontsize=6, xytext=(3, 2), textcoords='offset points')
    a.set(title=m, xlabel=f'PC1 ({x.PC1_pct.iloc[0]:.0f}%)', ylabel=f'PC2 ({x.PC2_pct.iloc[0]:.0f}%)')
for a, t in zip(ax.ravel(), ['A  Assigned counts', 'B  Genes detected', 'C', 'D']): a.set_title(t, loc='left', fontweight='bold')
f.legend(handles=handles, ncol=4, loc='outside lower center'); f.suptitle('GSE308825: library quality and sample structure', fontsize=9)
save(f, 'Fig1_GSE308825_QC')
# 2 signatures vs drug exposure
ag = pd.read_csv('results/GSE308825/sign_agreement.csv'); cam = pd.read_csv('results/GSE308825/signature_camera_fry.csv').query("analysis=='primary'")
cons = ['AML903_Quiz', 'AML001_Quiz', 'AML001_Palb', 'AML001_Combo']
f, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), layout='constrained', gridspec_kw={'width_ratios': [1.2, 1]})
a = ax[0]
for i, c in enumerate(cons):
    for off, s, m in [(-.16, 'S398', 'o'), (.16, 'S336', 's')]:
        r = ag[(ag.contrast == c) & (ag.signature == s)].iloc[0]; y = i + off
        a.plot([r.random_lo95, r.random_hi95], [y, y], color='#c3c2b7', lw=5, solid_capstyle='butt')
        a.plot(r.random_median_fraction, y, '|', color='#52514e', ms=7); a.plot(r.same_sign_fraction, y, m, color=BLUE, ms=5)
        a.annotate(f'FDR {r.empirical_FDR:.3f}', (.66, y), fontsize=5.5, va='center')
a.set(yticks=range(4), yticklabels=cons, xlim=(.38, .72), xlabel='fraction of signature genes with same sign\nas NM vs NC (grey: matched random 95%)', title='A  Gene-level sign agreement')
a.invert_yaxis()
a.legend(handles=[Line2D([], [], marker='o', ls='', color=BLUE, label='398'), Line2D([], [], marker='s', ls='', color=BLUE, label='336')], loc='upper center', bbox_to_anchor=(.5, -.32), ncol=2, fontsize=6)
b = ax[1]
M = np.array([[cam[(cam.contrast == c) & (cam.signature == s) & (cam.signature_direction == d)].camera_FDR.iloc[0] for s in ['S398', 'S336'] for d in ['higher_in_NM', 'higher_in_NC']] for c in cons])
b.imshow(M, cmap=LinearSegmentedColormap.from_list('b', ['#256abf', '#f0efec']), vmin=0, vmax=1, aspect='auto')
for (i, j), v in np.ndenumerate(M): b.text(j, i, f'{v:.2f}', ha='center', va='center', color='black' if v > .5 else 'white')
b.set(xticks=range(4), xticklabels=['398\nNM-up', '398\nNC-up', '336\nNM-up', '336\nNC-up'], yticks=range(4), yticklabels=cons, title='B  CAMERA FDR (none < 0.05)' if (M >= .05).all() else 'B  CAMERA FDR')
save(f, 'Fig2_signature_drug_exposure')
# 3 FOXM1
mv = pd.read_csv('results/foxm1_GSE274558_sample_values.csv').set_index('measure'); pv = pd.read_csv('results/GSE308825/foxm1_sample_values.csv')
sp = pd.read_csv('results/GSE248597/region_scores.csv'); types = ['LYMPHOIDAGGREGATE', 'MIXED1', 'MIXED2', 'CONTROL']; tl = ['aggregate', 'MIXED1', 'MIXED2', 'control']
f, ax = plt.subplots(2, 4, figsize=(7.2, 4.6), layout='constrained', gridspec_kw={'width_ratios': [.7, 1, 1.2, 1.3]})
rng = np.random.default_rng(1)
for row, (mk, pk, sk, lab) in enumerate([('Foxm1_mRNA_log2CPM', 'FOXM1_logCPM', 'FOXM1_log2_Q3', 'FOXM1 mRNA (log2)'), ('Foxm1_targets_signed', 'FOXM1_target_score', 'FOXM1_targets_signed', 'FOXM1 target score')]):
    a = ax[row, 0]
    for i, g in enumerate(['NC', 'NM']): a.scatter([i - .1, i, i + .1], mv.loc[mk, [f'{g}{k}' for k in (1, 2, 3)]], color=[BLUE, ORANGE][i], marker=['s', 'o'][i], s=18)
    a.set(xticks=[0, 1], xticklabels=['NPM1::\nCCDC28A', 'NPM1::\nMLF1'], xlim=(-.5, 1.5), ylabel=lab)
    for a, m, tr in [(ax[row, 1], 'AML903', ['DMSO', 'Quiz', 'Palb']), (ax[row, 2], 'AML001', ['DMSO', 'Quiz', 'Palb', 'Combo'])]:
        for i, t in enumerate(tr):
            v = pv[(pv.model == m) & (pv.treatment == t)][pk]; a.scatter(i + np.linspace(-.12, .12, len(v)) if len(v) > 1 else [i], v, color=TC[t], marker=MK[t], s=18)
        a.set(xticks=range(len(tr)), xticklabels=tr, xlim=(-.5, len(tr) - .5))
    a = ax[row, 3]
    for i, t in enumerate(types):
        v = sp[sp.region_type.str.replace('^AGGREGATE$', 'LYMPHOIDAGGREGATE', regex=True) == t][sk]; a.scatter(i + rng.uniform(-.22, .22, len(v)), v, color=GREY, s=4, alpha=.7, linewidths=0)
        a.plot([i - .3, i + .3], [v.median()] * 2, color='black', lw=1.2)
    a.set(xticks=range(4), xticklabels=[f'{x}' for x in tl]); a.tick_params(axis='x', rotation=45)
for a, t in zip(ax[:, 0], ['A', 'B']): a.text(-.75, 1.06, t, transform=a.transAxes, fontweight='bold', fontsize=9)
for a, t in zip(ax[0], ['GSE274558 (mouse)', 'AML903 xenograft', 'AML001 xenograft', 'GSE248597 regions']): a.set_title(t)
f.suptitle('FOXM1, exploratory. Scales differ between datasets and are not comparable. Spatial regions are not independent.', fontsize=7.5)
save(f, 'Fig3_FOXM1_across_datasets')
# 4 spatial heatmap
d = pd.read_csv('results/GSE248597/descriptive_by_slide_and_type.csv', dtype={'slide_barcode': str}); a_ = pd.read_csv('results/GSE248597/descriptive_by_region_type.csv')
names = ['T_cell', 'CD8_T', 'B_cell', 'Plasma_Ig', 'NK_cytotoxic', 'Monocyte_macrophage', 'Myeloid_progenitor_blast', 'Erythroid', 'Lymphoid_aggregate_chemokines']
cov = pd.read_csv('results/GSE248597/set_coverage.csv').set_index('set')
norm = lambda s: s.replace({'AGGREGATE': 'LYMPHOIDAGGREGATE'}); a_['region_type'] = norm(a_.region_type); d['region_type'] = norm(d.region_type)
cmap = LinearSegmentedColormap.from_list('s', ['#f4f8fd', '#9ec5f4', '#2a78d6', '#0d366b'])
f, ax = plt.subplots(1, 3, figsize=(7.2, 3.4), layout='constrained', sharey=True, gridspec_kw={'width_ratios': [4, 3, 4]})
for a, tab, tt, title in [(ax[0], a_, types, 'A  All regions'), (ax[1], d[d.slide_barcode.str.endswith('752')], types[:3], 'B  Slide ...752'), (ax[2], d[d.slide_barcode.str.endswith('859')], types, 'C  Slide ...859')]:
    piv = tab.pivot(index='measure', columns='region_type', values='mean').reindex(index=names, columns=tt)
    n = tab[tab.measure == names[0]].set_index('region_type').regions
    im = a.imshow(piv.values, cmap=cmap, vmin=0, vmax=1, aspect='auto')
    for (i, j), v in np.ndenumerate(piv.values): a.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=5.5, color='white' if v > .55 else 'black')
    a.set(title=title, xticks=range(len(tt)), xticklabels=[f'{tl[types.index(t)]} ({int(n[t])})' for t in tt]); a.tick_params(axis='x', rotation=60)
ax[0].set_yticks(range(len(names)), [f"{x.replace('_', ' ')} ({cov.genes_measured[x]}/{cov.genes_in_set[x]})" for x in names])
f.colorbar(im, ax=ax, label='mean within-region percentile rank', shrink=.7)
f.suptitle('GSE248597 immune and lineage marker sets: descriptive only (no biopsy identifiers)', fontsize=8)
save(f, 'Fig4_spatial_marker_sets')
# 5 spatial signatures vs random sets
r = pd.read_csv('results/GSE248597/random_set_reference.csv')
sel = ['S398_higher_in_NM', 'S398_higher_in_NC', 'S336_higher_in_NM', 'S336_higher_in_NC', 'FOXM1_targets_signed', 'proliferation_set']
f, ax = plt.subplots(2, 3, figsize=(7.2, 4.2), layout='constrained', sharex=True)
for a, s, letter in zip(ax.ravel(), sel, 'ABCDEF'):
    x = r[r['set'] == s].set_index('region_type').loc[types]
    for i, row in enumerate(x.itertuples()):
        a.plot([i, i], [row.random_lo95, row.random_hi95], color='#c3c2b7', lw=6, solid_capstyle='butt'); a.plot(i, row.mean_score, 'o', color=BLUE, ms=4)
    a.set(title=f"{letter}  {s.replace('_', ' ')} (n={x.genes.iloc[0]})", xticks=range(4), xticklabels=tl); a.tick_params(axis='x', rotation=45); a.title.set_fontsize(6.5)
f.supylabel('mean score (dot); expression-matched random sets, 95% (grey)', fontsize=7)
f.suptitle('Spatial scores against random gene sets (reference over genes, not over patients)', fontsize=8)
save(f, 'Fig5_spatial_signature_scores')
# 6 evidence matrix
e = pd.read_csv('results/evidence_matrix_long.csv'); claims = list(e.claim.unique()); D = ['GSE274558 mouse', 'GSE308825 AML903', 'GSE308825 AML001', 'GSE248597 spatial']
col = {'supportive': '#9ec5f4', 'non-significant': '#e1e0d9', 'discordant': '#f5b79c', 'unavailable': '#ffffff'}
lab = {'supportive': 'supportive', 'non-significant': 'non-sig.', 'discordant': 'discordant', 'unavailable': 'unavailable'}
f, a = plt.subplots(figsize=(7.2, .36 * len(claims) + 1))
import textwrap
for r_ in e.itertuples():
    i, j = claims.index(r_.claim), D.index(r_.dataset)
    a.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=col[r_.status], edgecolor='#c3c2b7', lw=.6))
    a.text(j + .5, i + .5, lab[r_.status] + ('*' if r_.status == 'supportive' else ''), ha='center', va='center', fontsize=6, color='#8a8984' if r_.status == 'unavailable' else 'black')
a.set_xlim(0, 4); a.set_ylim(len(claims), 0); a.xaxis.tick_top(); a.set_xticks(np.arange(4) + .5, [x.replace(' ', '\n') for x in D])
a.set_yticks(np.arange(len(claims)) + .5, ['\n'.join(textwrap.wrap(c, 62)) for c in claims], fontsize=5.8); a.tick_params(length=0)
for s in a.spines.values(): s.set_visible(False)
a.text(0, len(claims) + .35, '* narrow claim only: AML001 sign agreement is gene-level (not CAMERA-supported); the mouse FOXM1 result is not independent.', fontsize=5.8, va='top')
save(f, 'Fig6_evidence_matrix')
