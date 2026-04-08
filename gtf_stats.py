#!/usr/bin/env python3
"""
GTF文件统计脚本
统计内容：总基因数量、平均基因长度、每个基因平均外显子数、平均内含子长度
"""

import sys
from collections import defaultdict


def parse_attributes(attr_string):
    """解析GTF第9列的属性字段"""
    attributes = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if not item:
            continue
        # 处理 key "value" 格式
        parts = item.split(' ', 1)
        if len(parts) == 2:
            key = parts[0]
            value = parts[1].strip('"')
            attributes[key] = value
    return attributes


def parse_gtf(gtf_file):
    """解析GTF文件，返回基因信息字典"""
    genes = defaultdict(lambda: {
        'start': float('inf'),
        'end': 0,
        'strand': None,
        'exons': []  # 存储外显子坐标 [(start, end), ...]
    })
    
    with open(gtf_file, 'r') as f:
        for line in f:
            # 跳过注释行和空行
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom = fields[0]
            source = fields[1]
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = parse_attributes(fields[8])
            
            gene_id = attributes.get('gene_id')
            if not gene_id:
                continue
            
            # 更新基因的起始和结束位置
            genes[gene_id]['start'] = min(genes[gene_id]['start'], start)
            genes[gene_id]['end'] = max(genes[gene_id]['end'], end)
            genes[gene_id]['strand'] = strand
            
            # 记录外显子信息
            if feature == 'exon':
                genes[gene_id]['exons'].append((start, end))
    
    return genes


def calculate_introns(exons, strand):
    """根据外显子坐标计算内含子"""
    if len(exons) < 2:
        return []
    
    # 按起始位置排序外显子
    sorted_exons = sorted(exons, key=lambda x: x[0])
    
    introns = []
    for i in range(len(sorted_exons) - 1):
        intron_start = sorted_exons[i][1] + 1
        intron_end = sorted_exons[i + 1][0] - 1
        if intron_end >= intron_start:
            introns.append((intron_start, intron_end))
    
    return introns


def calculate_statistics(genes):
    """计算统计数据"""
    total_genes = len(genes)
    
    gene_lengths = []
    exon_counts = []
    intron_lengths = []
    
    for gene_id, gene_info in genes.items():
        # 基因长度
        gene_length = gene_info['end'] - gene_info['start'] + 1
        gene_lengths.append(gene_length)
        
        # 外显子数量
        exon_count = len(gene_info['exons'])
        exon_counts.append(exon_count)
        
        # 计算内含子
        introns = calculate_introns(gene_info['exons'], gene_info['strand'])
        for intron_start, intron_end in introns:
            intron_length = intron_end - intron_start + 1
            intron_lengths.append(intron_length)
    
    # 计算平均值
    avg_gene_length = sum(gene_lengths) / total_genes if total_genes > 0 else 0
    avg_exon_count = sum(exon_counts) / total_genes if total_genes > 0 else 0
    avg_intron_length = sum(intron_lengths) / len(intron_lengths) if intron_lengths else 0
    
    return {
        'total_genes': total_genes,
        'avg_gene_length': avg_gene_length,
        'avg_exon_count': avg_exon_count,
        'avg_intron_length': avg_intron_length,
        'total_exons': sum(exon_counts),
        'total_introns': len(intron_lengths)
    }


def main():
    if len(sys.argv) < 2:
        print("用法: python gtf_stats.py <GTF文件>")
        sys.exit(1)
    
    gtf_file = sys.argv[1]
    
    print(f"正在解析文件: {gtf_file}")
    print("-" * 50)
    
    # 解析GTF文件
    genes = parse_gtf(gtf_file)
    
    # 计算统计数据
    stats = calculate_statistics(genes)
    
    # 输出结果
    print(f"总基因数量:           {stats['total_genes']}")
    print(f"平均基因长度:         {stats['avg_gene_length']:.2f} bp")
    print(f"每个基因平均外显子数: {stats['avg_exon_count']:.2f}")
    print(f"平均内含子长度:       {stats['avg_intron_length']:.2f} bp")
    print("-" * 50)
    print(f"外显子总数:           {stats['total_exons']}")
    print(f"内含子总数:           {stats['total_introns']}")


if __name__ == "__main__":
    main()
