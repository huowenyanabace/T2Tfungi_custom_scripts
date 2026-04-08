#!/usr/bin/env python3
"""
GO Enrichment Analysis Visualization for Coprinellus radians (辐毛小鬼伞)
Expanded gene families - GO functional enrichment results
Generates 3 publication-quality figures:
  1. Bubble Plot (Top 20 enriched GO terms)
  2. GO Hierarchy DAG
  3. Gene-Function Network (Cnetplot)

Usage:
    python plot_go_enrichment.py <enrichment_csv> <go_annotation_csv>
"""

import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
import matplotlib.patches as mpatches
import networkx as nx
from collections import defaultdict
import re
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Global SCI Publication Style Settings
# ============================================================
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 9,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 8.5,
    'ytick.labelsize': 8.5,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'pdf.fonttype': 42,       # Editable text in PDF
    'ps.fonttype': 42,
})

# Color palette
CMAP = plt.cm.RdYlBu_r  # Red-Yellow-Blue reversed, common in enrichment plots


def load_data(enrichment_file, annotation_file):
    """Load and merge the two CSV files."""
    enrich = pd.read_csv(enrichment_file)
    annot = pd.read_csv(annotation_file)

    # Build GO ID -> Name mapping from annotation file
    id2name = {}
    id2ns = {}
    for _, row in annot.iterrows():
        gid = row['GO_ID']
        name = row['Name']
        ns = row['Namespace'] if pd.notna(row.get('Namespace')) else ''
        if pd.notna(name) and name != 'NOT FOUND':
            id2name[gid] = name
            id2ns[gid] = ns

    # Add term names to enrichment data
    enrich['Term'] = enrich['ID'].map(id2name)
    enrich['Namespace'] = enrich['ID'].map(id2ns)

    # For terms without annotation, use ID as name
    enrich['Term'] = enrich['Term'].fillna(enrich['ID'])
    enrich['Namespace'] = enrich['Namespace'].fillna('unknown')

    return enrich, annot, id2name, id2ns


def wrap_label(text, max_width=45):
    """Wrap long labels for y-axis."""
    if len(text) <= max_width:
        return text
    words = text.split()
    lines = []
    current = ""
    for w in words:
        if len(current) + len(w) + 1 <= max_width:
            current = current + " " + w if current else w
        else:
            lines.append(current)
            current = w
    if current:
        lines.append(current)
    return "\n".join(lines)


def namespace_short(ns):
    """Convert namespace to short form."""
    mapping = {
        'biological_process': 'BP',
        'cellular_component': 'CC',
        'molecular_function': 'MF',
    }
    return mapping.get(ns, 'Other')


# ============================================================
# Figure 1: Bubble Plot (Top 20 GO terms)
# ============================================================
def plot_bubble(enrich, output_prefix):
    """Generate a publication-quality bubble plot."""
    # Select top 20 by p.adjust
    df = enrich.nsmallest(20, 'p.adjust').copy()
    df = df.sort_values('RichFactor', ascending=True)

    # Prepare labels
    df['Label'] = df['Term'].apply(lambda t: wrap_label(t))
    df['NS_short'] = df['Namespace'].map(namespace_short)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))

    # Normalize p.adjust for colormap
    padj_vals = df['p.adjust'].values
    norm = mcolors.LogNorm(vmin=max(padj_vals.min(), 1e-10), vmax=padj_vals.max())

    # Bubble sizes
    counts = df['Count'].values
    size_min, size_max = 40, 200
    if counts.max() == counts.min():
        sizes = np.full_like(counts, (size_min + size_max) / 2, dtype=float)
    else:
        sizes = size_min + (counts - counts.min()) / (counts.max() - counts.min()) * (size_max - size_min)

    # Namespace markers and colors
    ns_markers = {'BP': 'o', 'CC': 's', 'MF': 'D', 'Other': '^'}

    for ns in df['NS_short'].unique():
        mask = df['NS_short'] == ns
        subset = df[mask]
        idx = np.where(mask)[0]
        sc = ax.scatter(
            subset['RichFactor'],
            np.arange(len(df))[idx],
            s=sizes[idx],
            c=subset['p.adjust'].values,
            cmap=CMAP,
            norm=norm,
            marker=ns_markers.get(ns, 'o'),
            edgecolors='#333333',
            linewidths=0.5,
            alpha=0.9,
            zorder=3,
            label=ns,
        )

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['Label'].values)
    ax.set_xlabel('Rich Factor', fontweight='bold')
    # Title removed, saved to separate txt file

    # No gridlines for clean SCI look
    ax.yaxis.grid(False)
    ax.xaxis.grid(False)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_label('Adjusted p-value', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # Size legend
    unique_counts = sorted(df['Count'].unique())
    legend_sizes = []
    legend_labels = []
    for c in unique_counts:
        if counts.max() == counts.min():
            s = (size_min + size_max) / 2
        else:
            s = size_min + (c - counts.min()) / (counts.max() - counts.min()) * (size_max - size_min)
        legend_sizes.append(s)
        legend_labels.append(str(c))

    size_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
               markersize=np.sqrt(s) * 0.45, markeredgecolor='#333', markeredgewidth=0.5,
               label=f'Count = {lab}')
        for s, lab in zip(legend_sizes, legend_labels)
    ]

    # Namespace legend
    ns_handles = [
        Line2D([0], [0], marker=ns_markers[ns], color='w', markerfacecolor='grey',
               markersize=7, markeredgecolor='#333', markeredgewidth=0.5,
               label=ns)
        for ns in sorted(df['NS_short'].unique())
    ]

    all_handles = ns_handles + size_handles
    leg = ax.legend(
        handles=all_handles,
        loc='lower right',
        frameon=True,
        framealpha=0.9,
        edgecolor='#cccccc',
        fontsize=7.5,
        handletextpad=0.5,
        borderpad=0.6,
        title='Category / Gene Count',
        title_fontsize=8,
    )

    plt.tight_layout()
    outpath = f"{output_prefix}_bubble_plot.pdf"
    fig.savefig(outpath, format='pdf')
    outpath_png = f"{output_prefix}_bubble_plot.png"
    fig.savefig(outpath_png, format='png', dpi=300)
    plt.close(fig)
    print(f"  Bubble plot saved: {outpath}, {outpath_png}")
    return outpath_png


# ============================================================
# Figure 2: GO Hierarchy DAG
# ============================================================
def plot_dag(enrich, annot, id2name, output_prefix):
    """Generate a GO hierarchy DAG with enrichment info."""
    G = nx.DiGraph()

    # Build padj mapping
    padj_map = dict(zip(enrich['ID'], enrich['p.adjust']))
    count_map = dict(zip(enrich['ID'], enrich['Count']))
    gene_ratio_map = dict(zip(enrich['ID'], enrich['GeneRatio']))

    # Parse parent-child relationships from annotation file
    enriched_ids = set(enrich['ID'].values)
    nodes_in_dag = set()

    for _, row in annot.iterrows():
        gid = row['GO_ID']
        name = row['Name']
        parents_str = row['Parents'] if pd.notna(row.get('Parents')) else ''
        status = row['Status']

        if status == 'Not Found':
            continue

        node_label = f"{name}\n({gid})" if pd.notna(name) and name != 'NOT FOUND' else gid
        G.add_node(gid, label=node_label, name=name if pd.notna(name) else gid)
        nodes_in_dag.add(gid)

        if parents_str:
            parents = re.findall(r'(GO:\d+)', str(parents_str))
            for pid in parents:
                if pid not in G:
                    pname = id2name.get(pid, pid)
                    plabel = f"{pname}\n({pid})"
                    G.add_node(pid, label=plabel, name=pname)
                G.add_edge(pid, gid)  # parent -> child

    # Only keep enriched terms and their connecting paths
    enriched_in_graph = enriched_ids & set(G.nodes())

    # Find all nodes on paths between enriched terms
    relevant_nodes = set()
    enriched_list = list(enriched_in_graph)
    for node in enriched_list:
        relevant_nodes.add(node)
        # Add ancestors that connect to other enriched terms
        for ancestor in nx.ancestors(G, node):
            if ancestor in G:
                relevant_nodes.add(ancestor)

    # Prune to relevant subgraph but keep it manageable
    # Keep only nodes that are enriched or directly connect enriched nodes
    keep_nodes = set()
    for node in enriched_in_graph:
        keep_nodes.add(node)
    for node in relevant_nodes:
        # Check if this node connects enriched terms
        descendants = nx.descendants(G, node) & enriched_in_graph
        if len(descendants) >= 1:
            keep_nodes.add(node)

    # Further limit: keep top 25 enriched + their direct parents
    top_enriched = enrich.nsmallest(25, 'p.adjust')['ID'].values
    core_nodes = set()
    for gid in top_enriched:
        if gid in G:
            core_nodes.add(gid)
            for parent in G.predecessors(gid):
                core_nodes.add(parent)
                # Add grandparents too for context
                for gp in G.predecessors(parent):
                    core_nodes.add(gp)

    subG = G.subgraph(core_nodes & set(G.nodes())).copy()

    if len(subG.nodes()) < 2:
        print("  Warning: Too few nodes for DAG, skipping.")
        return None

    # Assign depth via topological sorting for hierarchical layout
    # Compute depth for each node (longest path from any root)
    roots = [n for n in subG.nodes() if subG.in_degree(n) == 0]
    depth_map = {}
    for node in nx.topological_sort(subG):
        preds = list(subG.predecessors(node))
        if not preds:
            depth_map[node] = 0
        else:
            depth_map[node] = max(depth_map.get(p, 0) for p in preds) + 1

    max_depth = max(depth_map.values()) if depth_map else 1

    # Group nodes by depth
    depth_groups = defaultdict(list)
    for node, d in depth_map.items():
        depth_groups[d].append(node)

    # Create hierarchical positions
    pos = {}
    for d, nodes_at_d in depth_groups.items():
        n_nodes = len(nodes_at_d)
        for i, node in enumerate(sorted(nodes_at_d)):
            x = (i - (n_nodes - 1) / 2) * 2.0  # spread horizontally
            y = -d * 2.5  # top to bottom
            pos[node] = (x, y)

    # Apply spring adjustment to reduce overlap while keeping hierarchy
    # (fix y, relax x)
    for _ in range(50):
        forces = {n: [0.0, 0.0] for n in subG.nodes()}
        nodes_list = list(subG.nodes())
        for i in range(len(nodes_list)):
            for j in range(i + 1, len(nodes_list)):
                n1, n2 = nodes_list[i], nodes_list[j]
                dx = pos[n1][0] - pos[n2][0]
                dy = pos[n1][1] - pos[n2][1]
                dist = max(np.sqrt(dx**2 + dy**2), 0.5)
                if dist < 3.0:
                    repulsion = 1.5 / (dist ** 2)
                    forces[n1][0] += repulsion * dx / dist
                    forces[n2][0] -= repulsion * dx / dist
        for node in subG.nodes():
            pos[node] = (pos[node][0] + forces[node][0] * 0.1,
                         pos[node][1])  # only adjust x

    fig, ax = plt.subplots(figsize=(14, 11))

    # Node colors based on p.adjust
    enriched_padj = [padj_map[n] for n in subG.nodes() if n in padj_map]
    if enriched_padj:
        norm = mcolors.LogNorm(vmin=max(min(enriched_padj), 1e-10), vmax=max(enriched_padj))
    else:
        norm = mcolors.LogNorm(vmin=1e-10, vmax=0.05)

    node_colors = []
    node_sizes = []
    node_edge_colors = []
    for node in subG.nodes():
        if node in padj_map:
            node_colors.append(CMAP(norm(padj_map[node])))
            node_sizes.append(1200)
            node_edge_colors.append('#222222')
        else:
            node_colors.append('#E0E0E0')
            node_sizes.append(650)
            node_edge_colors.append('#999999')

    # Draw edges
    nx.draw_networkx_edges(
        subG, pos, ax=ax,
        edge_color='#999999',
        arrows=True,
        arrowsize=12,
        arrowstyle='-|>',
        width=0.9,
        alpha=0.5,
        min_source_margin=20,
        min_target_margin=20,
    )

    # Draw nodes
    nx.draw_networkx_nodes(
        subG, pos, ax=ax,
        node_color=node_colors,
        node_size=node_sizes,
        edgecolors=node_edge_colors,
        linewidths=0.8,
        alpha=0.92,
    )

    # Assign numeric labels to all nodes
    node_list = list(subG.nodes())
    # Enriched nodes first (sorted by padj), then parental nodes
    enriched_nodes_sorted = sorted(
        [n for n in node_list if n in padj_map],
        key=lambda x: padj_map[x]
    )
    parental_nodes_sorted = sorted(
        [n for n in node_list if n not in padj_map]
    )
    ordered_nodes = enriched_nodes_sorted + parental_nodes_sorted

    node_num_map = {}  # node_id -> number
    for idx, node in enumerate(ordered_nodes, start=1):
        node_num_map[node] = idx

    # Save number-to-name legend file
    legend_path = f"{output_prefix}_DAG_node_legend.txt"
    with open(legend_path, 'w') as f:
        f.write("DAG Node Legend for C. radians GO Enrichment\n")
        f.write("=" * 60 + "\n")
        f.write(f"{'No.':<6}{'GO ID':<14}{'Term Name':<55}{'Type'}\n")
        f.write("-" * 60 + "\n")
        for node in ordered_nodes:
            num = node_num_map[node]
            name = subG.nodes[node].get('name', node)
            ntype = "Enriched" if node in padj_map else "Parental"
            f.write(f"{num:<6}{node:<14}{name:<55}{ntype}\n")
    print(f"  DAG node legend saved: {legend_path}")

    # Draw labels: numbers centered inside nodes
    for node in subG.nodes():
        num = node_num_map[node]
        x, y = pos[node]
        fontcolor = 'white' if node in padj_map else '#333'
        ax.text(
            x, y, str(num),
            fontsize=10,
            fontweight='bold',
            ha='center', va='center',
            color=fontcolor,
            fontfamily='DejaVu Sans',
            zorder=7,
        )

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.35, aspect=20, pad=0.01)
    cbar.set_label('Adjusted p-value', fontsize=9)
    cbar.ax.tick_params(labelsize=7)

    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=CMAP(0.2),
               markersize=10, markeredgecolor='#222', markeredgewidth=0.8,
               label='Enriched GO term'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#E0E0E0',
               markersize=8, markeredgecolor='#999', markeredgewidth=0.8,
               label='Parental GO term'),
    ]
    ax.legend(handles=legend_elements, loc='lower left', fontsize=8,
              frameon=True, framealpha=0.9, edgecolor='#ccc')

    # Title removed, saved to separate txt file
    ax.axis('off')

    plt.tight_layout()
    outpath = f"{output_prefix}_GO_DAG.pdf"
    fig.savefig(outpath, format='pdf')
    outpath_png = f"{output_prefix}_GO_DAG.png"
    fig.savefig(outpath_png, format='png', dpi=300)
    plt.close(fig)
    print(f"  GO DAG saved: {outpath}, {outpath_png}")
    return outpath_png


# ============================================================
# Figure 3: Gene-Function Network (Cnetplot)
# ============================================================
def plot_cnetplot(enrich, id2name, output_prefix):
    """Generate a gene-GO term network (Cnetplot style)."""
    # Strategy: include top 10 ESCRT-related terms + all singleton gene terms
    # to show both the dominant pattern and the diverse individual genes
    top = enrich.nsmallest(10, 'p.adjust').copy()

    # Also add terms associated with singleton genes (Count == 1)
    singletons = enrich[enrich['Count'] == 1].nsmallest(8, 'p.adjust')
    top = pd.concat([top, singletons]).drop_duplicates(subset='ID')

    G = nx.Graph()
    padj_map_local = dict(zip(top['ID'], top['p.adjust']))

    all_genes = set()
    go_terms = []

    for _, row in top.iterrows():
        go_id = row['ID']
        term = id2name.get(go_id, go_id)
        if len(term) > 38:
            term = term[:36] + '...'
        go_label = f"{term}\n({go_id})"
        go_terms.append(go_id)
        G.add_node(go_id, ntype='go', label=go_label, padj=row['p.adjust'],
                   count=row['Count'])

        genes = str(row['geneID']).split('/')
        for gene in genes:
            gene = gene.strip()
            if gene:
                short_gene = gene.replace('-RA', '')
                G.add_node(gene, ntype='gene', label=short_gene)
                G.add_edge(gene, go_id)
                all_genes.add(gene)

    # Separate gene and GO nodes
    gene_nodes = [n for n in G.nodes() if G.nodes[n]['ntype'] == 'gene']
    go_nodes = [n for n in G.nodes() if G.nodes[n]['ntype'] == 'go']

    # Custom layout: genes in center, GO terms around them
    # Use bipartite-like layout
    pos = {}
    n_go = len(go_nodes)

    # Place GO terms in outer ring
    radius_go = 4.0
    for i, node in enumerate(go_nodes):
        angle = 2 * np.pi * i / n_go - np.pi / 2
        pos[node] = (radius_go * np.cos(angle), radius_go * np.sin(angle))

    # Place genes: compute average position of connected GO terms
    for gene in gene_nodes:
        neighbors = list(G.neighbors(gene))
        if neighbors:
            avg_x = np.mean([pos[n][0] for n in neighbors if n in pos])
            avg_y = np.mean([pos[n][1] for n in neighbors if n in pos])
            # Pull toward center
            pos[gene] = (avg_x * 0.35, avg_y * 0.35)

    # Add jitter to separate overlapping genes
    gene_positions = [(g, pos[g]) for g in gene_nodes]
    for i in range(len(gene_positions)):
        for j in range(i + 1, len(gene_positions)):
            g1, p1 = gene_positions[i]
            g2, p2 = gene_positions[j]
            dx = p1[0] - p2[0]
            dy = p1[1] - p2[1]
            dist = np.sqrt(dx**2 + dy**2)
            if dist < 0.8:
                angle = np.arctan2(dy, dx) if dist > 0 else np.random.uniform(0, 2*np.pi)
                push = (0.8 - dist) / 2
                pos[g1] = (pos[g1][0] + push * np.cos(angle),
                           pos[g1][1] + push * np.sin(angle))
                pos[g2] = (pos[g2][0] - push * np.cos(angle),
                           pos[g2][1] - push * np.sin(angle))

    fig, ax = plt.subplots(figsize=(12, 10))

    # Color GO nodes by p.adjust
    go_padj = [G.nodes[n]['padj'] for n in go_nodes]
    if go_padj:
        norm = mcolors.LogNorm(vmin=max(min(go_padj), 1e-10), vmax=max(go_padj))
    else:
        norm = mcolors.LogNorm(vmin=1e-10, vmax=0.05)

    go_colors = [CMAP(norm(G.nodes[n]['padj'])) for n in go_nodes]

    # Gene degrees and colors
    gene_degrees = {g: G.degree(g) for g in gene_nodes}
    max_deg = max(gene_degrees.values()) if gene_degrees else 1
    gene_sizes = [250 + 350 * (gene_degrees[g] / max_deg) for g in gene_nodes]

    gene_cmap = plt.cm.YlOrRd
    gene_colors = [gene_cmap(0.3 + 0.6 * gene_degrees[g] / max_deg) for g in gene_nodes]

    # Draw edges with gradient color
    for u, v in G.edges():
        gene_n = u if G.nodes[u]['ntype'] == 'gene' else v
        go_n = v if G.nodes[v]['ntype'] == 'go' else u
        edge_color = CMAP(norm(G.nodes[go_n]['padj']))
        ax.plot(
            [pos[u][0], pos[v][0]],
            [pos[u][1], pos[v][1]],
            color=edge_color, alpha=0.30, linewidth=1.5, zorder=1
        )

    # Draw GO nodes (squares), size by count
    go_x = [pos[n][0] for n in go_nodes]
    go_y = [pos[n][1] for n in go_nodes]
    go_sizes = [300 + 200 * G.nodes[n].get('count', 1) for n in go_nodes]
    ax.scatter(go_x, go_y, s=go_sizes, c=go_colors, marker='s',
               edgecolors='#333', linewidths=0.8, zorder=3, alpha=0.92)

    # Draw Gene nodes (circles)
    gene_x = [pos[n][0] for n in gene_nodes]
    gene_y = [pos[n][1] for n in gene_nodes]
    ax.scatter(gene_x, gene_y, s=gene_sizes, c=gene_colors, marker='o',
               edgecolors='#444', linewidths=0.7, zorder=3, alpha=0.9)

    # Labels for GO terms - place outside the ring
    for n in go_nodes:
        label = G.nodes[n]['label']
        x, y = pos[n]
        # Push label further outward
        dist = np.sqrt(x**2 + y**2)
        if dist > 0:
            offset_x = 18 * x / dist
            offset_y = 18 * y / dist
        else:
            offset_x, offset_y = 15, 0

        ha = 'left' if x >= 0 else 'right'

        ax.annotate(
            label, (x, y),
            textcoords='offset points', xytext=(offset_x, offset_y),
            fontsize=6.0, color='#222',
            fontweight='normal',
            ha=ha, va='center',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                      edgecolor='#ddd', alpha=0.85),
            arrowprops=dict(arrowstyle='-', color='#bbb', lw=0.5),
            zorder=5,
        )

    # Labels for Gene nodes
    for n in gene_nodes:
        label = G.nodes[n]['label']
        x, y = pos[n]
        ax.annotate(
            label, (x, y),
            textcoords='offset points', xytext=(0, -13),
            fontsize=7, color='#333',
            fontweight='bold',
            ha='center', va='top',
            zorder=5,
        )

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.35, aspect=20, pad=0.01)
    cbar.set_label('Adjusted p-value', fontsize=9)
    cbar.ax.tick_params(labelsize=7)

    # Legend
    legend_elements = [
        Line2D([0], [0], marker='s', color='w', markerfacecolor=CMAP(0.3),
               markersize=11, markeredgecolor='#333', markeredgewidth=0.8,
               label='GO Term'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=gene_cmap(0.6),
               markersize=9, markeredgecolor='#444', markeredgewidth=0.6,
               label='Gene (hub)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=gene_cmap(0.3),
               markersize=7, markeredgecolor='#444', markeredgewidth=0.6,
               label='Gene (single)'),
    ]
    ax.legend(handles=legend_elements, loc='lower left', fontsize=8.5,
              frameon=True, framealpha=0.9, edgecolor='#ccc')

    ax.set_title('Gene–GO Term Network (Cnetplot)\nExpanded Gene Families in C. radians',
                 fontweight='bold', fontsize=11)
    ax.axis('off')

    plt.tight_layout()
    outpath = f"{output_prefix}_cnetplot.pdf"
    fig.savefig(outpath, format='pdf')
    outpath_png = f"{output_prefix}_cnetplot.png"
    fig.savefig(outpath_png, format='png', dpi=300)
    plt.close(fig)
    print(f"  Cnetplot saved: {outpath}, {outpath_png}")
    return outpath_png


# ============================================================
# Main
# ============================================================
def main():
    if len(sys.argv) < 3:
        print("Usage: python plot_go_enrichment.py <enrichment.csv> <go_annotation.csv> [output_prefix]")
        print("Example: python plot_go_enrichment.py Crad_GO_enrichment_results.csv Crad_results.csv Crad_GO")
        sys.exit(1)

    enrichment_file = sys.argv[1]
    annotation_file = sys.argv[2]
    output_prefix = sys.argv[3] if len(sys.argv) >= 4 else "Crad_GO"

    print("=" * 60)
    print("GO Enrichment Visualization for C. radians")
    print("=" * 60)

    # Load data
    print("\n[1/4] Loading data...")
    enrich, annot, id2name, id2ns = load_data(enrichment_file, annotation_file)
    print(f"  Enriched GO terms: {len(enrich)}")
    print(f"  Annotated GO terms: {len(annot)}")
    print(f"  Terms with names: {len(id2name)}")

    # Figure 1: Bubble Plot
    print("\n[2/4] Generating Bubble Plot...")
    bubble_png = plot_bubble(enrich, output_prefix)

    # Figure 2: GO DAG
    print("\n[3/4] Generating GO Hierarchy DAG...")
    dag_png = plot_dag(enrich, annot, id2name, output_prefix)

    # Save figure titles to txt file
    print("\n[4/4] Saving figure titles...")
    titles_path = f"{output_prefix}_figure_titles.txt"
    with open(titles_path, 'w') as f:
        f.write("Figure Titles for C. radians GO Enrichment Analysis\n")
        f.write("=" * 60 + "\n\n")
        f.write("Figure 1 (Bubble Plot):\n")
        f.write("Top 20 enriched GO terms of expanded gene families "
                "in Coprinellus radians. Bubble size represents gene count, "
                "color indicates adjusted p-value, and shape denotes GO "
                "category (BP: biological process; CC: cellular component; "
                "MF: molecular function).\n\n")
        f.write("Figure 2 (GO DAG):\n")
        f.write("Directed acyclic graph (DAG) showing the hierarchical "
                "relationships among enriched GO terms in expanded gene "
                "families of Coprinellus radians. Colored nodes represent "
                "significantly enriched GO terms (color mapped to adjusted "
                "p-value), and grey nodes represent parental GO terms "
                "providing hierarchical context.\n")
    print(f"  Titles saved: {titles_path}")

    print("\n" + "=" * 60)
    print("All figures generated successfully!")
    print("Output files:")
    print(f"  - {output_prefix}_bubble_plot.pdf / .png")
    print(f"  - {output_prefix}_GO_DAG.pdf / .png")
    print(f"  - {output_prefix}_DAG_node_legend.txt")
    print(f"  - {output_prefix}_figure_titles.txt")
    print("=" * 60)


if __name__ == '__main__':
    main()
