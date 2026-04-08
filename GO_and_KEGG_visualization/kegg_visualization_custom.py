#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
定制版KEGG富集分析结果可视化脚本
只生成bubble_plot和top_pathways两个图表，输出PDF格式
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import argparse
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# 设置图形参数
plt.style.use('default')
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.weight'] = 'normal'
plt.rcParams['axes.labelweight'] = 'normal'
plt.rcParams['axes.titleweight'] = 'normal'
plt.rcParams['figure.titleweight'] = 'normal'
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.figsize'] = (12, 8)

class KEGGVisualizer:
    """定制版KEGG富集分析可视化类"""
    
    def __init__(self, input_file, output_dir="./output", prefix="KEGG"):
        self.input_file = input_file
        self.output_dir = Path(output_dir)
        self.prefix = prefix
        
        # 创建输出目录
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 图表标题（将保存到txt文件）
        self.titles = {}
        
    def load_data(self):
        """加载和预处理数据"""
        print(f"Loading data from: {self.input_file}")
        
        # 读取数据
        df = pd.read_csv(self.input_file)
        
        # 检查必需的列
        required_columns = ['ID', 'Count', 'FoldEnrichment', 'p.adjust']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")
        
        # 如果没有Description列，使用ID作为描述
        if 'Description' not in df.columns:
            df['Description'] = df['ID']
        
        # 计算-log10(p.adjust)
        df['neg_log10_padj'] = -np.log10(df['p.adjust'])
        
        # 按p.adjust排序
        df = df.sort_values('p.adjust').reset_index(drop=True)
        
        self.data = df
        print(f"Loaded {len(df)} pathways")
        return df
    
    def generate_colors(self, n_colors):
        """生成n种不同的颜色"""
        if n_colors <= 10:
            # 使用tab10颜色映射
            colors = plt.cm.tab10(np.linspace(0, 1, min(n_colors, 10)))
        else:
            # 使用更多颜色的映射
            colors = plt.cm.tab20(np.linspace(0, 1, min(n_colors, 20)))
            if n_colors > 20:
                # 如果需要更多颜色，使用连续色谱
                additional_colors = plt.cm.Set3(np.linspace(0, 1, n_colors - 20))
                colors = np.vstack([colors, additional_colors])
        
        return colors[:n_colors]
    
    def create_bubble_plot(self, top_n=10):
        """创建富集分析气泡图"""
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # 选择前N个最显著的通路
        plot_data = self.data.head(top_n)
        
        # 生成颜色
        colors = self.generate_colors(len(plot_data))
        
        # 创建散点图
        scatter_plots = []
        labels = []
        
        for i, (idx, row) in enumerate(plot_data.iterrows()):
            scatter = ax.scatter(
                row['FoldEnrichment'],
                row['neg_log10_padj'],
                s=row['Count'] * 15,  # 基因数量决定气泡大小
                c=[colors[i]],
                alpha=0.7,
                edgecolors='black',
                linewidth=0.5
            )
            scatter_plots.append(scatter)
            
            # 简化图例标签：只保留ID
            labels.append(row['ID'])
        
        # 设置轴标签
        ax.set_xlabel('Fold Enrichment', fontsize=12)
        ax.set_ylabel('-log10(Adjusted P-value)', fontsize=12)
        
        # 添加图例 - 去除边框，增大行间距，稍微左移
        ax.legend(scatter_plots, labels, bbox_to_anchor=(1.02, 1), loc='upper left', 
                 frameon=False, fontsize=8, labelspacing=2.0)
        
        # 移除网格
        ax.grid(False)
        
        plt.tight_layout()
        
        # 保存标题
        self.titles['bubble_plot'] = f'KEGG Pathway Enrichment Analysis - Top {top_n} Most Significant Pathways'
        
        # 保存图片为PDF
        output_file = self.output_dir / f"{self.prefix}_bubble_plot.pdf"
        plt.savefig(output_file, format='pdf', bbox_inches='tight')
        print(f"Saved bubble plot: {output_file}")
        
        plt.show()
        plt.close()
    
    def create_top_pathways_plot(self, top_n=10):
        """创建Top富集通路条形图"""
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # 选择前N个最显著的通路
        plot_data = self.data.head(top_n).copy()
        
        # 截断过长的通路名称
        plot_data['Short_Name'] = plot_data['Description'].apply(
            lambda x: str(x)[:50] + '...' if len(str(x)) > 50 else str(x)
        )
        
        # 生成颜色，颜色按富集倍数映射
        norm = plt.Normalize(vmin=plot_data['FoldEnrichment'].min(), 
                            vmax=plot_data['FoldEnrichment'].max())
        colors = plt.cm.viridis(norm(plot_data['FoldEnrichment']))
        
        # 创建水平条形图
        bars = ax.barh(range(len(plot_data)), plot_data['FoldEnrichment'], 
                      color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
        
        # 设置Y轴标签
        ax.set_yticks(range(len(plot_data)))
        ax.set_yticklabels(plot_data['Short_Name'][::-1], fontsize=10)  # 反转顺序
        ax.invert_yaxis()  # 最显著的在上面
        
        # 设置标签
        ax.set_xlabel('Fold Enrichment', fontsize=12)
        
        # 移除网格
        ax.grid(False)
        
        plt.tight_layout()
        
        # 保存标题
        self.titles['top_pathways'] = f'Top {top_n} Most Significantly Enriched KEGG Pathways'
        
        # 保存图片为PDF
        output_file = self.output_dir / f"{self.prefix}_top_pathways.pdf"
        plt.savefig(output_file, format='pdf', bbox_inches='tight')
        print(f"Saved top pathways plot: {output_file}")
        
        plt.show()
        plt.close()
    
    def save_titles(self):
        """保存图表标题到txt文件"""
        output_file = self.output_dir / f"{self.prefix}_titles.txt"
        with open(output_file, 'w') as f:
            f.write("KEGG Enrichment Analysis - Figure Titles\n")
            f.write("=" * 50 + "\n\n")
            for key, title in self.titles.items():
                f.write(f"{key}:\n{title}\n\n")
        
        print(f"Saved titles: {output_file}")
    
    def run_visualization(self):
        """运行可视化"""
        print("=== KEGG Enrichment Analysis Visualization ===")
        
        # 加载数据
        self.load_data()
        
        print("\n1. Creating bubble plot...")
        self.create_bubble_plot()
        
        print("\n2. Creating top pathways plot...")
        self.create_top_pathways_plot()
        
        print("\n3. Saving titles...")
        self.save_titles()
        
        print(f"\n=== Visualization Complete ===")
        print(f"All plots saved to: {self.output_dir}")
        
        return self.data

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='Customized KEGG Enrichment Analysis Visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python kegg_visualization_custom.py --input enrichment_results.csv --output_dir ./plots
  python kegg_visualization_custom.py --input data.csv --output_dir ./results --prefix MyAnalysis
        """
    )
    
    parser.add_argument('--input', '-i', required=True,
                       help='Input CSV file with KEGG enrichment results')
    parser.add_argument('--output_dir', '-o', default='./output',
                       help='Output directory for plots (default: ./output)')
    parser.add_argument('--prefix', '-p', default='KEGG',
                       help='Prefix for output files (default: KEGG)')
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found!")
        return
    
    # 创建可视化器
    visualizer = KEGGVisualizer(
        input_file=args.input,
        output_dir=args.output_dir,
        prefix=args.prefix
    )
    
    # 运行可视化
    try:
        results = visualizer.run_visualization()
        print(f"\nAnalysis completed successfully!")
        print(f"Check the output directory: {args.output_dir}")
    except Exception as e:
        print(f"Error during visualization: {e}")
        print("Please check your input file format and try again.")

if __name__ == "__main__":
    main()
