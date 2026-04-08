#!/usr/bin/env python3
"""
CAZyme Comparative Genomics Visualization Script
=================================================
Generate a heatmap and a fold-change horizontal bar chart from a CAZyme
gene-count TSV file, plus a figure_legends.txt with titles and explanations.

Usage:
    python plot_cazyme_param.py <input.tsv> [output_dir]

Arguments:
    input.tsv   : Path to the CAZyme gene-count TSV file.
                  Required format — first column: CAZyme family names (e.g. GH5,
                  AA9, CBM1 …); remaining columns: gene counts per organism.
    output_dir  : (Optional) Directory for output files. Default: ./output

Outputs:
    <output_dir>/heatmap.png          — Heatmap of key lignocellulose CAZyme families
    <output_dir>/foldchange_bar.png   — Fold-change horizontal bar chart
    <output_dir>/figure_legends.txt   — Figure titles and explanations

Example:
    python plot_cazyme_param.py Cazyme.tsv ./results
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch


# ─────────────────────────────────────────────────────────────
# 0. Parse arguments
# ─────────────────────────────────────────────────────────────
if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(1)

INPUT_TSV  = sys.argv[1]
OUTPUT_DIR = sys.argv[2] if len(sys.argv) >= 3 else './output'
os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"[INFO] Input TSV : {INPUT_TSV}")
print(f"[INFO] Output dir: {OUTPUT_DIR}")


# ─────────────────────────────────────────────────────────────
# 1. Font setup (try Noto Sans CJK, fallback to DejaVu Sans)
# ─────────────────────────────────────────────────────────────
_noto_regular = '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
_noto_bold    = '/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc'

if os.path.exists(_noto_regular):
    fp      = fm.FontProperties(fname=_noto_regular)
    fp_bold = fm.FontProperties(fname=_noto_bold if os.path.exists(_noto_bold) else _noto_regular)
    plt.rcParams['font.family'] = fp.get_name()
else:
    fp = fp_bold = fm.FontProperties()
    print("[WARN] Noto Sans CJK not found; falling back to default font.")

plt.rcParams['axes.unicode_minus'] = False


# ─────────────────────────────────────────────────────────────
# 2. Load data
# ─────────────────────────────────────────────────────────────
df = pd.read_csv(INPUT_TSV, sep='\t')
# Normalise the first column name
first_col = df.columns[0]
df = df.rename(columns={first_col: 'Family'})

all_orgs = [c for c in df.columns if c != 'Family']
print(f"[INFO] Detected {len(all_orgs)} organisms: {all_orgs}")
print(f"[INFO] Total CAZyme families: {len(df)}")

# ── Target / Others split ──
# Default: Cdom, Crad, Cmic as target group; everything else as others.
# Modify here if your species of interest differ.
TARGET_SPECIES = ['Cdom', 'Crad', 'Cmic']
target = [o for o in TARGET_SPECIES if o in all_orgs]
others = [o for o in all_orgs if o not in target]

if not target:
    print(f"[WARN] None of {TARGET_SPECIES} found in columns. "
          f"Using first 3 organisms as target group.")
    target = all_orgs[:3]
    others = all_orgs[3:]

print(f"[INFO] Target group : {target}")
print(f"[INFO] Others group : {others}")


# ─────────────────────────────────────────────────────────────
# 3. CHART 1 — HEATMAP
# ─────────────────────────────────────────────────────────────
hm_families = [
    # Cellulose GH
    'GH5','GH6','GH7','GH9','GH12','GH3','GH1',
    # Hemicellulose GH
    'GH10','GH11','GH43','GH30','GH16','GH74',
    'GH51','GH62','GH79','GH115','GH35',
    # Carbohydrate Esterases
    'CE3','CE4','CE1','CE15',
    # Auxiliary Activities (lignin + LPMO)
    'AA1','AA2','AA3','AA5','AA9','AA14',
    'AA7','AA8',
    # Carbohydrate-Binding Modules
    'CBM1','CBM2','CBM5',
]

# Keep only families present in the data
hm_present = [f for f in hm_families if f in df['Family'].values]
hm_df = df[df['Family'].isin(hm_present)].set_index('Family').loc[hm_present, all_orgs]
hm_z  = hm_df.apply(lambda r: (r - r.mean()) / (r.std() + 0.01), axis=1)

fig, ax = plt.subplots(figsize=(max(8, len(all_orgs) * 1.2), max(10, len(hm_present) * 0.42)))
cmap = LinearSegmentedColormap.from_list(
    'rbu', ['#2B6CB0','#EBF8FF','#FFFFF0','#FED7AA','#C53030'], N=256)
im = ax.imshow(hm_z.values, cmap=cmap, aspect='auto', vmin=-2, vmax=2)

ax.set_xticks(range(len(all_orgs)))
ax.set_xticklabels(all_orgs, fontsize=12, fontproperties=fp)
ax.set_yticks(range(len(hm_present)))
ax.set_yticklabels(hm_present, fontsize=11, fontproperties=fp)

# Highlight target species
for i, org in enumerate(all_orgs):
    if org in target:
        ax.get_xticklabels()[i].set_fontweight('bold')
        ax.get_xticklabels()[i].set_color('#C53030')

# Raw gene counts inside cells
for i in range(len(hm_present)):
    for j in range(len(all_orgs)):
        val = int(hm_df.iloc[i, j])
        color = 'white' if abs(hm_z.iloc[i, j]) > 1.3 else 'black'
        ax.text(j, i, str(val), ha='center', va='center',
                fontsize=9, color=color, fontweight='bold')

cbar = plt.colorbar(im, ax=ax, shrink=0.4, pad=0.02)
cbar.set_label('Row Z-score', fontsize=11, fontproperties=fp)

plt.tight_layout()
heatmap_path = os.path.join(OUTPUT_DIR, 'heatmap.png')
plt.savefig(heatmap_path, dpi=250, bbox_inches='tight', facecolor='white')
plt.close()
print(f"[OK]   Heatmap saved -> {heatmap_path}")


# ─────────────────────────────────────────────────────────────
# 4. CHART 2 — FOLD CHANGE HORIZONTAL BAR
# ─────────────────────────────────────────────────────────────
target_label = '+'.join(target)

df2 = df.copy()
df2 = df2[df2[all_orgs].sum(axis=1) > 0]
df2['t_mean'] = df2[target].mean(axis=1)
df2['o_mean'] = df2[others].mean(axis=1)
df2['fc']     = df2['t_mean'] / df2['o_mean'].replace(0, 0.1)

# Filter: target mean > others mean * 1.3  AND  target mean >= 2
df2 = df2[(df2['t_mean'] > df2['o_mean'] * 1.3) & (df2['t_mean'] >= 2)]
df2 = df2.sort_values('fc', ascending=True)

n_bars = len(df2)
fig, ax = plt.subplots(figsize=(11, max(7, n_bars * 0.38)))
bars = ax.barh(range(n_bars), df2['fc'].values, height=0.7, zorder=3)

for bar, fc in zip(bars, df2['fc'].values):
    if   fc >= 2.5: bar.set_color('#C53030'); bar.set_alpha(0.9)
    elif fc >= 2.0: bar.set_color('#DD6B20'); bar.set_alpha(0.85)
    elif fc >= 1.7: bar.set_color('#D69E2E'); bar.set_alpha(0.8)
    else:           bar.set_color('#38A169'); bar.set_alpha(0.75)

ax.set_yticks(range(n_bars))
ax.set_yticklabels(df2['Family'].values, fontsize=11, fontproperties=fp)

for i, (_, row) in enumerate(df2.iterrows()):
    fc = row['fc']
    ax.text(fc + 0.03, i, f'{fc:.1f}x',
            va='center', fontsize=10, fontweight='bold', color='#2D3748')
    ax.text(fc + 0.28, i, f'(avg {row["t_mean"]:.1f} vs {row["o_mean"]:.1f})',
            va='center', fontsize=8, color='#718096')

ax.axvline(x=1.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax.text(1.5, n_bars - 0.3, '1.5x', fontsize=9, color='gray', ha='center')

ax.set_xlabel(f'Fold Change ({target_label} mean / Others mean)',
              fontsize=12, fontproperties=fp, labelpad=6)
ax.set_xlim(1.0, max(df2['fc']) + 0.6)
ax.set_ylim(-0.6, n_bars - 0.4)
ax.grid(axis='x', alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

legend_elements = [
    Patch(facecolor='#C53030', alpha=0.9,  label='>=2.5x'),
    Patch(facecolor='#DD6B20', alpha=0.85, label='2.0-2.5x'),
    Patch(facecolor='#D69E2E', alpha=0.8,  label='1.7-2.0x'),
    Patch(facecolor='#38A169', alpha=0.75, label='1.3-1.7x'),
]
ax.legend(handles=legend_elements, loc='center right', fontsize=10,
          framealpha=0.9, prop=fp, bbox_to_anchor=(0.72, 0.5))

plt.tight_layout()
plt.subplots_adjust(bottom=0.05)
fc_path = os.path.join(OUTPUT_DIR, 'foldchange_bar.png')
plt.savefig(fc_path, dpi=250, bbox_inches='tight', facecolor='white')
plt.close()
print(f"[OK]   Fold-change bar saved -> {fc_path}")


# ─────────────────────────────────────────────────────────────
# 5. TXT — Figure legends
# ─────────────────────────────────────────────────────────────
fc_families = df2.sort_values('fc', ascending=False)
fc_table = ""
for _, r in fc_families.iterrows():
    fc_table += (f"  {r['Family']:<10s}  {r['fc']:.1f}x   "
                 f"({target_label} avg = {r['t_mean']:.1f}, "
                 f"Others avg = {r['o_mean']:.1f})\n")

txt = f"""\
================================================================================
Figure 1. Heatmap of Key Lignocellulose-Degrading CAZyme Families
================================================================================

Title (EN): Heatmap of key lignocellulose-degrading CAZyme gene counts across
            {len(all_orgs)} fungal species
Title (CN): {len(all_orgs)}种真菌木质纤维素降解关键CAZyme家族基因数热图

Description:
  The heatmap displays gene copy numbers for {len(hm_present)} CAZyme families involved
  in lignocellulose degradation across {len(all_orgs)} fungal species. Color intensity
  reflects row-wise Z-score normalization (red = above-average; blue =
  below-average within each family). Raw gene counts are shown in each cell.
  The target species ({target_label}) are highlighted in red on the x-axis.

  Rows are grouped into five functional categories (top to bottom):
    1. Cellulose GH   (GH5, GH6, GH7, GH9, GH12, GH3, GH1)
       - Endoglucanases, cellobiohydrolases, beta-glucosidases
    2. Hemicellulose GH   (GH10, GH11, GH43, GH30, GH16, GH74, GH51, GH62,
       GH79, GH115, GH35)
       - Endo-xylanases, xylosidases, arabinofuranosidases, xyloglucanases,
         glucuronidases, galactosidases
    3. Carbohydrate Esterases   (CE3, CE4, CE1, CE15)
       - Acetyl xylan esterases, feruloyl esterases,
         glucuronoyl methylesterases (LCC cleavage)
    4. Auxiliary Activities   (AA1, AA2, AA3, AA5, AA9, AA14, AA7, AA8)
       - Laccases, peroxidases (MnP/LiP/VP), GMC oxidoreductases,
         LPMOs (cellulose oxidative cleavage), auxiliary electron donors
    5. Carbohydrate-Binding Modules   (CBM1, CBM2, CBM5)
       - Cellulose-binding (fungal/bacterial type), chitin-binding


================================================================================
Figure 2. Fold Change of Dominant CAZyme Families in {target_label}
================================================================================

Title (EN): CAZyme families with elevated gene counts in {target_label}
            compared to other species (fold change ranking)
Title (CN): {target_label}优势CAZyme家族倍数变化排名

Description:
  Horizontal bar chart showing CAZyme families where the average gene count
  of {target_label} exceeds that of the remaining {len(others)} species ({', '.join(others)}).
  Fold change = mean({target_label}) / mean(Others).
  Bars are colored by fold change magnitude:
    - Red (>=2.5x): very strong enrichment
    - Orange (2.0-2.5x): strong enrichment
    - Yellow (1.7-2.0x): notable enrichment
    - Green (1.3-1.7x): moderate enrichment

  Each bar is annotated with the fold change value and the average gene
  counts (target avg vs others avg).
  Filtering criteria: target mean > others mean x 1.3, and target mean >= 2.

Ranked families (descending fold change):
{fc_table}"""

txt_path = os.path.join(OUTPUT_DIR, 'figure_legends.txt')
with open(txt_path, 'w', encoding='utf-8') as f:
    f.write(txt)
print(f"[OK]   Legends saved -> {txt_path}")
print(f"\n[DONE] All outputs in: {OUTPUT_DIR}")
