#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OrthoFinder结果统计分析脚本 (v2)

功能：
1. 统计每个物种在同源群中存在的基因数
2. 统计物种特异性的基因数（只存在于本物种的基因）
3. 统计每个物种的单拷贝基因数
4. 统计总共识别出的直系同源基因家族数
5. 统计核心基因家族数（存在于所有物种中）
6. 统计仅在指定的多个物种之间共享的同源群（支持2个、3个或更多物种）

作者：Claude
日期：2025
"""

import argparse
import pandas as pd
from collections import defaultdict
from datetime import datetime


def load_gene_count_file(filepath):
    """加载Orthogroups_GeneCount.tsv文件"""
    df = pd.read_csv(filepath, sep='\t', index_col=0)
    # 移除Total列（如果存在）
    if 'Total' in df.columns:
        df = df.drop(columns=['Total'])
    return df


def load_orthogroups_file(filepath):
    """加载Orthogroups.tsv文件"""
    df = pd.read_csv(filepath, sep='\t', index_col=0)
    return df


def count_genes_per_species(gene_count_df):
    """统计每个物种在同源群中存在的基因总数"""
    return gene_count_df.sum().to_dict()


def count_species_specific_genes(gene_count_df):
    """
    统计物种特异性基因数
    物种特异性基因：该同源群只在一个物种中存在基因
    """
    species_specific = defaultdict(int)
    
    for og_id, row in gene_count_df.iterrows():
        # 找出哪些物种在这个同源群中有基因
        species_with_genes = row[row > 0]
        
        # 如果只有一个物种有基因，则为物种特异性
        if len(species_with_genes) == 1:
            species_name = species_with_genes.index[0]
            species_specific[species_name] += species_with_genes.values[0]
    
    # 确保所有物种都在结果中
    for species in gene_count_df.columns:
        if species not in species_specific:
            species_specific[species] = 0
    
    return dict(species_specific)


def count_single_copy_genes(gene_count_df):
    """
    统计每个物种的单拷贝基因数
    单拷贝基因：在该物种中该同源群只有1个基因
    """
    single_copy = {}
    
    for species in gene_count_df.columns:
        single_copy[species] = (gene_count_df[species] == 1).sum()
    
    return single_copy


def count_total_orthogroups(gene_count_df):
    """统计总共识别出的直系同源基因家族数"""
    return len(gene_count_df)


def count_core_orthogroups(gene_count_df):
    """
    统计核心基因家族数
    核心基因：存在于所有物种中（每个物种至少有1个基因）
    """
    core_count = 0
    core_orthogroups = []
    
    for og_id, row in gene_count_df.iterrows():
        if (row > 0).all():
            core_count += 1
            core_orthogroups.append(og_id)
    
    return core_count, core_orthogroups


def find_shared_orthogroups(gene_count_df, species_list):
    """
    找出仅在指定的多个物种之间共享的同源群
    即：只有这些物种有基因，其他物种都没有
    
    参数：
        gene_count_df: 基因计数数据框
        species_list: 指定的物种列表（可以是2个、3个或更多）
    """
    shared_only = []
    target_species_set = set(species_list)
    
    for og_id, row in gene_count_df.iterrows():
        # 检查是否只有指定的物种有基因
        species_with_genes = set(row[row > 0].index.tolist())
        
        if species_with_genes == target_species_set:
            shared_only.append(og_id)
    
    return shared_only


def generate_report(gene_count_df, target_species_list, output_file):
    """生成统计报告"""
    
    # 获取所有物种名称
    species_list = gene_count_df.columns.tolist()
    
    # 进行各项统计
    total_genes = count_genes_per_species(gene_count_df)
    species_specific = count_species_specific_genes(gene_count_df)
    single_copy = count_single_copy_genes(gene_count_df)
    total_orthogroups = count_total_orthogroups(gene_count_df)
    core_count, core_orthogroups = count_core_orthogroups(gene_count_df)
    shared_only = find_shared_orthogroups(gene_count_df, target_species_list)
    
    # 生成报告
    report_lines = []
    report_lines.append("=" * 70)
    report_lines.append("OrthoFinder 结果统计报告")
    report_lines.append("=" * 70)
    report_lines.append(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append("")
    
    # 1. 基本信息
    report_lines.append("-" * 70)
    report_lines.append("1. 基本信息")
    report_lines.append("-" * 70)
    report_lines.append(f"分析的物种数量: {len(species_list)}")
    report_lines.append(f"物种列表: {', '.join(species_list)}")
    report_lines.append(f"直系同源基因家族总数: {total_orthogroups}")
    report_lines.append("")
    
    # 2. 每个物种的统计
    report_lines.append("-" * 70)
    report_lines.append("2. 各物种统计")
    report_lines.append("-" * 70)
    report_lines.append(f"{'物种':<20}{'基因数':<15}{'物种特异性基因数':<20}{'单拷贝基因数':<15}")
    report_lines.append("-" * 70)
    
    for species in species_list:
        report_lines.append(
            f"{species:<20}{total_genes[species]:<15}{species_specific[species]:<20}{single_copy[species]:<15}"
        )
    
    report_lines.append("")
    
    # 3. 核心基因统计
    report_lines.append("-" * 70)
    report_lines.append("3. 核心基因家族统计（存在于所有物种中）")
    report_lines.append("-" * 70)
    report_lines.append(f"核心基因家族数量: {core_count}")
    report_lines.append(f"占总家族比例: {core_count/total_orthogroups*100:.2f}%")
    report_lines.append("")
    
    # 如果核心基因家族数量不多，可以列出来
    if core_count <= 100:
        report_lines.append("核心基因家族列表:")
        for i in range(0, len(core_orthogroups), 10):
            report_lines.append("  " + ", ".join(core_orthogroups[i:i+10]))
    else:
        report_lines.append(f"（核心基因家族数量较多，不在此列出）")
    
    report_lines.append("")
    
    # 4. 指定物种间共享的同源群
    report_lines.append("-" * 70)
    species_names = " 和 ".join(target_species_list) if len(target_species_list) <= 3 else \
                    ", ".join(target_species_list[:-1]) + " 和 " + target_species_list[-1]
    report_lines.append(f"4. 仅在 {species_names} 之间共享的同源群")
    report_lines.append(f"   （指定物种数: {len(target_species_list)}）")
    report_lines.append("-" * 70)
    report_lines.append(f"共享同源群数量: {len(shared_only)}")
    report_lines.append("")
    
    if len(shared_only) > 0:
        report_lines.append("共享同源群列表:")
        for i in range(0, len(shared_only), 10):
            report_lines.append("  " + ", ".join(shared_only[i:i+10]))
    else:
        report_lines.append("（没有找到仅在这些物种间共享的同源群）")
    
    report_lines.append("")
    report_lines.append("=" * 70)
    report_lines.append("报告结束")
    report_lines.append("=" * 70)
    
    # 写入文件
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    print(f"统计报告已保存到: {output_file}")
    
    # 同时打印到控制台
    print('\n'.join(report_lines))
    
    return {
        'total_genes': total_genes,
        'species_specific': species_specific,
        'single_copy': single_copy,
        'total_orthogroups': total_orthogroups,
        'core_count': core_count,
        'core_orthogroups': core_orthogroups,
        'shared_only': shared_only
    }


def main():
    parser = argparse.ArgumentParser(
        description='OrthoFinder结果统计分析脚本',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 指定2个物种
  python orthogroup_stats_v2.py -i Orthogroups_GeneCount.tsv -o stats.txt -s Cdom.pep Crad.pep
  
  # 指定3个物种
  python orthogroup_stats_v2.py -i Orthogroups_GeneCount.tsv -o stats.txt -s Cdom.pep Crad.pep Ccin.pep
  
  # 指定4个或更多物种
  python orthogroup_stats_v2.py -i Orthogroups_GeneCount.tsv -o stats.txt -s Cdom.pep Crad.pep Ccin.pep Cmar.pep
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='输入文件路径 (Orthogroups_GeneCount.tsv)'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='orthogroup_statistics.txt',
        help='输出报告文件路径 (默认: orthogroup_statistics.txt)'
    )
    
    parser.add_argument(
        '-s', '--species',
        nargs='+',
        required=True,
        help='指定要查找共享同源群的物种列表（至少2个，用空格分隔）'
    )
    
    args = parser.parse_args()
    
    # 检查物种数量
    if len(args.species) < 2:
        print("错误: 至少需要指定2个物种")
        return
    
    # 加载数据
    print(f"正在加载文件: {args.input}")
    gene_count_df = load_gene_count_file(args.input)
    
    # 验证指定的物种是否存在
    all_species = gene_count_df.columns.tolist()
    
    invalid_species = []
    for sp in args.species:
        if sp not in all_species:
            invalid_species.append(sp)
    
    if invalid_species:
        print(f"错误: 以下物种不在数据中: {', '.join(invalid_species)}")
        print(f"可用的物种: {', '.join(all_species)}")
        return
    
    # 生成报告
    print(f"正在生成统计报告...")
    print(f"查找仅在以下物种间共享的同源群: {', '.join(args.species)}")
    generate_report(gene_count_df, args.species, args.output)


if __name__ == '__main__':
    main()
