#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
定制版GO富集分析结果可视化脚本
只生成bubble_plot和top_terms两个图表，输出PDF格式
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

class GOVisualizer:
    """定制版GO富集分析可视化类"""
    
    def __init__(self, enrichment_file, description_file, output_dir="./output", prefix="GO"):
        self.enrichment_file = enrichment_file
        self.description_file = description_file
        self.output_dir = Path(output_dir)
        self.prefix = prefix
        
        # 创建输出目录
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 标准GO命名空间颜色映射
        self.namespace_colors = {
            'biological_process': '#E31A1C',
            'molecular_function': '#1F78B4', 
            'cellular_component': '#33A02C',
            'unknown': '#999999'
        }
        
        # 图表标题（将保存到txt文件）
        self.titles = {}
    
    def load_data(self):
        """加载和合并数据"""
        print(f"Loading enrichment data from: {self.enrichment_file}")
        print(f"Loading description data from: {self.description_file}")
        
        # 读取富集数据
        enrichment_df = pd.read_csv(self.enrichment_file)
        
        # 检查必需的列
        required_enrich_columns = ['ID', 'Count', 'FoldEnrichment', 'p.adjust']
        missing_columns = [col for col in required_enrich_columns if col not in enrichment_df.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns in enrichment file: {missing_columns}")
        
        # 读取描述数据
        description_df = pd.read_csv(self.description_file)
        
        # 检查描述文件的必需列
        required_desc_columns = ['GO_ID']
        if not any(col in description_df.columns for col in required_desc_columns):
            raise ValueError("Description file must contain 'GO_ID' column")
        
        # 合并数据
        merged_df = enrichment_df.merge(
            description_df, 
            left_on='ID', right_on='GO_ID', how='left'
        )
        
        # 处理缺失的描述
        if 'Name' in description_df.columns:
            merged_df['Name'] = merged_df['Name'].fillna(merged_df['ID'])
        else:
            merged_df['Name'] = merged_df['ID']
            
        if 'Namespace' in description_df.columns:
            merged_df['Namespace'] = merged_df['Namespace'].fillna('unknown')
        else:
            merged_df['Namespace'] = 'unknown'
        
        # 计算-log10(p.adjust)
        merged_df['neg_log10_padj'] = -np.log10(merged_df['p.adjust'])
        
        # 按p.adjust排序
        merged_df = merged_df.sort_values('p.adjust').reset_index(drop=True)
        
        self.data = merged_df
        print(f"Loaded {len(merged_df)} GO terms")
        return merged_df
    
    def generate_colors_by_namespace(self, plot_data):
        """根据命名空间和条目生成颜色"""
        colors = []
        namespace_counts = {}
        
        for _, row in plot_data.iterrows():
            namespace = row['Namespace']
            base_color = self.namespace_colors.get(namespace, '#999999')
            
            # 为每个命名空间内的条目生成不同的色调
            if namespace not in namespace_counts:
                namespace_counts[namespace] = 0
            else:
                namespace_counts[namespace] += 1
            
            # 通过调整色调来生成同一命名空间内的不同颜色
            base_rgb = mcolors.to_rgb(base_color)
            # 添加一些变化来区分同命名空间内的不同条目
            variation = namespace_counts[namespace] * 0.1
            adjusted_rgb = tuple(max(0, min(1, c + variation * (0.5 - c))) for c in base_rgb)
            colors.append(adjusted_rgb)
        
        return colors
    
    def create_namespace_bubble_plot(self, top_n=10):
        """创建按命名空间分色的GO富集气泡图"""
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # 选择前N个最显著的GO条目
        plot_data = self.data.head(top_n)
        
        # 为每个GO条目生成不同的颜色
        colors = plt.cm.tab20(np.linspace(0, 1, len(plot_data)))
        
        # 创建散点图
        scatter_plots = []
        labels = []
        
        for i, (idx, row) in enumerate(plot_data.iterrows()):
            scatter = ax.scatter(
                row['FoldEnrichment'],
                row['neg_log10_padj'],
                s=row['Count'] * 20,
                c=[colors[i]],
                alpha=0.7,
                edgecolors='black',
                linewidth=0.5
            )
            scatter_plots.append(scatter)
            
            # 简化图例标签：只显示ID和命名空间分类
            namespace_abbr = {
                'biological_process': 'BP',
                'molecular_function': 'MF',
                'cellular_component': 'CC',
                'unknown': 'Unknown'
            }.get(row['Namespace'], 'Unknown')
            
            labels.append(f"{row['ID']} ({namespace_abbr})")
        
        # 设置轴标签
        ax.set_xlabel('Fold Enrichment', fontsize=12)
        ax.set_ylabel('-log10(Adjusted P-value)', fontsize=12)
        
        # 添加图例 - 去除边框，增大行间距
        ax.legend(scatter_plots, labels, bbox_to_anchor=(1.05, 1), loc='upper left', 
                 frameon=False, fontsize=8, labelspacing=2.0)
        
        # 移除网格
        ax.grid(False)
        
        plt.tight_layout()
        
        # 保存标题
        self.titles['bubble_plot'] = f'GO Enrichment Analysis by Namespace - Top {top_n} Most Significant Terms'
        
        # 保存图片为PDF
        output_file = self.output_dir / f"{self.prefix}_bubble_plot.pdf"
        plt.savefig(output_file, format='pdf', bbox_inches='tight')
        print(f"Saved namespace bubble plot: {output_file}")
        
        plt.show()
        plt.close()
    
    def create_top_terms_plot(self, top_n=10):
        """创建Top富集GO条目图"""
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # 选择前N个最显著的GO条目
        plot_data = self.data.head(top_n).copy()
        
        # 截断过长的GO名称
        plot_data['Short_Name'] = plot_data['Name'].apply(
            lambda x: (str(x)[:70] + '...') if len(str(x)) > 70 else str(x)
        )
        
        # 根据命名空间分配颜色
        colors = [self.namespace_colors.get(ns, '#999999') for ns in plot_data['Namespace']]
        
        # 创建水平条形图
        y_pos = np.arange(len(plot_data))
        bars = ax.barh(y_pos, plot_data['FoldEnrichment'], color=colors, 
                      alpha=0.8, edgecolor='black', linewidth=0.5)
        
        # 设置Y轴
        ax.set_yticks(y_pos)
        ax.set_yticklabels(plot_data['Short_Name'], fontsize=9)
        ax.invert_yaxis()
        
        # 设置标签
        ax.set_xlabel('Fold Enrichment', fontsize=12)
        
        # 创建图例（只有多个命名空间时）- 去除边框
        if len(plot_data['Namespace'].unique()) > 1:
            legend_elements = []
            for namespace, color in self.namespace_colors.items():
                if namespace in plot_data['Namespace'].values:
                    legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=color, 
                                                       label=namespace.replace('_', ' ').title()))
            
            if legend_elements:
                ax.legend(handles=legend_elements, title='GO Namespace',
                         bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
        
        # 移除网格
        ax.grid(False)
        
        plt.tight_layout()
        
        # 保存标题
        self.titles['top_terms'] = f'Top {top_n} Most Significantly Enriched GO Terms'
        
        # 保存图片为PDF
        output_file = self.output_dir / f"{self.prefix}_top_terms.pdf"
        plt.savefig(output_file, format='pdf', bbox_inches='tight')
        print(f"Saved top terms plot: {output_file}")
        
        plt.show()
        plt.close()
    
    def save_titles(self):
        """保存图表标题到txt文件"""
        output_file = self.output_dir / f"{self.prefix}_titles.txt"
        with open(output_file, 'w') as f:
            f.write("GO Enrichment Analysis - Figure Titles\n")
            f.write("=" * 50 + "\n\n")
            for key, title in self.titles.items():
                f.write(f"{key}:\n{title}\n\n")
        
        print(f"Saved titles: {output_file}")
    
    def run_visualization(self):
        """运行可视化"""
        print("=== GO Enrichment Analysis Visualization ===")
        
        # 加载数据
        self.load_data()
        
        print("\n1. Creating namespace bubble plot...")
        self.create_namespace_bubble_plot()
        
        print("\n2. Creating top terms plot...")
        self.create_top_terms_plot()
        
        print("\n3. Saving titles...")
        self.save_titles()
        
        print(f"\n=== Visualization Complete ===")
        print(f"All plots saved to: {self.output_dir}")
        
        return self.data

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='Customized GO Enrichment Analysis Visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python go_visualization_custom.py --enrichment enrichment.csv --description description.csv --output_dir ./plots
  python go_visualization_custom.py --enrichment results.csv --description go_desc.csv --output_dir ./results --prefix MyGO
        """
    )
    
    parser.add_argument('--enrichment', '-e', required=True,
                       help='Input CSV file with GO enrichment results')
    parser.add_argument('--description', '-d', required=True,
                       help='Input CSV file with GO descriptions')
    parser.add_argument('--output_dir', '-o', default='./output',
                       help='Output directory for plots (default: ./output)')
    parser.add_argument('--prefix', '-p', default='GO',
                       help='Prefix for output files (default: GO)')
    
    args = parser.parse_args()
    
    # 检查输入文件
    for file_path in [args.enrichment, args.description]:
        if not os.path.exists(file_path):
            print(f"Error: Input file '{file_path}' not found!")
            return
    
    # 创建可视化器
    visualizer = GOVisualizer(
        enrichment_file=args.enrichment,
        description_file=args.description,
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
