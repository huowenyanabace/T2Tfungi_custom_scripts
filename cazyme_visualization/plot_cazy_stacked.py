#!/usr/bin/env python3
"""
CAZyme Class Stacked Percentage Bar Chart
==========================================
Visualize the proportional composition of CAZyme classes across species.

Usage:
    python plot_cazy_stacked.py <input.tsv> [output_dir]

Arguments:
    input.tsv   : TSV with first column = CAZyme class (AA, CBM, CE, GH, GT, PL),
                  remaining columns = gene counts per organism.
    output_dir  : (Optional) Output directory. Default: ./output
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# ─────────────────────────────────────────────────────────────
# 0. Parse arguments
# ─────────────────────────────────────────────────────────────
if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(1)

INPUT_TSV  = sys.argv[1]
OUTPUT_DIR = sys.argv[2] if len(sys.argv) >= 3 else './output'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────
# 1. Font setup
# ─────────────────────────────────────────────────────────────
_noto = '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
if os.path.exists(_noto):
    fp = fm.FontProperties(fname=_noto)
    plt.rcParams['font.family'] = fp.get_name()
else:
    fp = fm.FontProperties()
plt.rcParams['axes.unicode_minus'] = False

# ─────────────────────────────────────────────────────────────
# 2. Load data
# ─────────────────────────────────────────────────────────────
df = pd.read_csv(INPUT_TSV, sep='\t')
first_col = df.columns[0]
df = df.rename(columns={first_col: 'Class'})
organisms = [c for c in df.columns if c != 'Class']
df = df.set_index('Class')

# Calculate percentages per organism
df_pct = df.div(df.sum(axis=0), axis=1) * 100

print(f"[INFO] Input: {INPUT_TSV}")
print(f"[INFO] Classes: {list(df.index)}")
print(f"[INFO] Organisms: {organisms}")

# ─────────────────────────────────────────────────────────────
# 3. Color palette for CAZyme classes
# ─────────────────────────────────────────────────────────────
class_colors = {
    'AA':  '#E53E3E',   # red
    'CBM': '#DD6B20',   # orange
    'CE':  '#D69E2E',   # yellow
    'GH':  '#38A169',   # green
    'GT':  '#3182CE',   # blue
    'PL':  '#805AD5',   # purple
}
# Fallback for unknown classes
fallback_colors = ['#FC8181','#F6AD55','#F6E05E','#68D391','#63B3ED','#B794F4',
                   '#FEB2B2','#FCCB8F']
colors = []
for cls in df.index:
    colors.append(class_colors.get(cls, fallback_colors[len(colors) % len(fallback_colors)]))

# ─────────────────────────────────────────────────────────────
# 4. Plot stacked percentage bar chart
# ─────────────────────────────────────────────────────────────
n_orgs = len(organisms)
x = np.arange(n_orgs)
bar_width = 0.6

fig, ax = plt.subplots(figsize=(max(8, n_orgs * 1.1), 6.5))

bottom = np.zeros(n_orgs)
bars_list = []

for i, cls in enumerate(df.index):
    vals = df_pct.loc[cls, organisms].values.astype(float)
    b = ax.bar(x, vals, bar_width, bottom=bottom, color=colors[i],
               edgecolor='white', linewidth=0.5, label=cls)
    bars_list.append((cls, vals, bottom.copy()))
    bottom += vals

# Add percentage labels inside bars (only if segment >= 4%)
for cls, vals, bot in bars_list:
    for j in range(n_orgs):
        if vals[j] >= 4:
            y_pos = bot[j] + vals[j] / 2
            ax.text(x[j], y_pos, f'{vals[j]:.1f}%',
                    ha='center', va='center', fontsize=8,
                    fontweight='bold', color='white')

# Add total gene count on top of each bar
for j, org in enumerate(organisms):
    total = int(df[org].sum())
    ax.text(x[j], 101, f'n={total}', ha='center', va='bottom',
            fontsize=9, color='#4A5568', fontproperties=fp)

ax.set_xticks(x)
ax.set_xticklabels(organisms, fontsize=12, fontproperties=fp)
ax.set_ylabel('Proportion (%)', fontsize=13, fontproperties=fp)
ax.set_ylim(0, 108)
ax.set_xlim(-0.5, n_orgs - 0.5)

# Legend
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),
          ncol=len(df.index), fontsize=10, frameon=False, prop=fp)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', alpha=0.3, linestyle='--')

plt.tight_layout()
out_path = os.path.join(OUTPUT_DIR, 'cazy_stacked_bar.png')
plt.savefig(out_path, dpi=250, bbox_inches='tight', facecolor='white')
plt.close()
print(f"[OK]   Stacked bar saved -> {out_path}")

# ─────────────────────────────────────────────────────────────
# 5. TXT legend
# ─────────────────────────────────────────────────────────────
pct_table = df_pct.round(1).to_string()
raw_table = df.to_string()

txt = f"""\
================================================================================
Figure. Stacked Percentage Bar Chart of CAZyme Class Composition
================================================================================

Title (EN): Proportional composition of CAZyme classes across {n_orgs} fungal species
Title (CN): {n_orgs}种真菌CAZyme类别组成百分比堆叠柱状图

Description:
  Stacked bar chart showing the relative proportion (%) of six CAZyme classes
  (AA, CBM, CE, GH, GT, PL) in each species. The total gene count (n=) is
  displayed above each bar. Percentage labels are shown inside segments >= 4%.

  CAZyme class definitions:
    AA  — Auxiliary Activities (redox enzymes, LPMOs)
    CBM — Carbohydrate-Binding Modules (non-catalytic adhesion)
    CE  — Carbohydrate Esterases (deacetylation of polysaccharides)
    GH  — Glycoside Hydrolases (hydrolysis of glycosidic bonds)
    GT  — Glycosyltransferases (formation of glycosidic bonds)
    PL  — Polysaccharide Lyases (non-hydrolytic cleavage)

Raw gene counts:
{raw_table}

Percentage composition:
{pct_table}
"""

txt_path = os.path.join(OUTPUT_DIR, 'figure_legend_stacked.txt')
with open(txt_path, 'w', encoding='utf-8') as f:
    f.write(txt)
print(f"[OK]   Legend saved -> {txt_path}")
print(f"\n[DONE] All outputs in: {OUTPUT_DIR}")
