# Usage (archive root): python3 scripts/01b_foxm1_human_targets.py   (needs signatures/foxm1_targets_mouse.csv from script 02)
import pandas as pd
o = pd.read_csv('inputs/Mouse_Human_Orthologs_Ensembl116.tsv', sep='\t', quoting=3); o.columns = ['mid','msym','hid','hsym','type','conf']
one = o[(o.type == 'ortholog_one2one') & (o.conf == 1)][['msym','hsym']].drop_duplicates().drop_duplicates('msym')
fx = pd.read_csv('signatures/foxm1_targets_mouse.csv')
h = fx.merge(one, left_on='target', right_on='msym', how='left').drop(columns='msym').rename(columns={'target': 'mouse_symbol', 'hsym': 'human_symbol'})
h.to_csv('signatures/foxm1_targets_human.csv', index=False)
