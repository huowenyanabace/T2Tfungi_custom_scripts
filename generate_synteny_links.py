#!/usr/bin/env python3
"""
生成共线性link文件

该脚本读取BED文件和anchors文件，输出共线性link文件。
共线性link文件格式：
第1-3列：第1列共线性区块所在的染色体名称，起始位置，终止位置
第4-6列：第2列共线性区块所在的染色体名称，起始位置，终止位置
"""

import sys
from collections import defaultdict


def parse_bed(bed_file):
    """
    解析BED文件，建立基因ID到位置信息的映射
    
    返回: 字典 {gene_id: (chrom, start, end)}
    """
    gene_info = {}
    with open(bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 4:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                gene_id = parts[3]
                gene_info[gene_id] = (chrom, start, end)
    return gene_info


def parse_anchors(anchors_file):
    """
    解析anchors文件，提取共线性区块信息
    
    返回: 共线性区块列表，每个区块是一个列表，包含(gene1, gene2)元组
    """
    blocks = []
    current_block = []
    
    with open(anchors_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('###'):
                # 遇到新的区块分隔符
                if current_block:
                    blocks.append(current_block)
                current_block = []
            else:
                parts = line.split('\t')
                if len(parts) >= 2:
                    gene1 = parts[0]
                    gene2 = parts[1]
                    current_block.append((gene1, gene2))
    
    # 添加最后一个区块
    if current_block:
        blocks.append(current_block)
    
    return blocks


def generate_synteny_links(bed_file, anchors_file, output_file):
    """
    生成共线性link文件
    
    对于每个共线性区块，计算：
    - 第1列基因的染色体、最小起始位置、最大终止位置
    - 第2列基因的染色体、最小起始位置、最大终止位置
    """
    # 解析BED文件
    gene_info = parse_bed(bed_file)
    print(f"从BED文件中读取了 {len(gene_info)} 个基因的位置信息")
    
    # 解析anchors文件
    blocks = parse_anchors(anchors_file)
    print(f"从anchors文件中读取了 {len(blocks)} 个共线性区块")
    
    # 生成link文件
    with open(output_file, 'w') as f:
        for i, block in enumerate(blocks):
            # 收集第1列和第2列基因的位置信息
            genes1_info = []
            genes2_info = []
            
            for gene1, gene2 in block:
                if gene1 in gene_info:
                    genes1_info.append(gene_info[gene1])
                else:
                    print(f"警告: 基因 {gene1} 未在BED文件中找到")
                
                if gene2 in gene_info:
                    genes2_info.append(gene_info[gene2])
                else:
                    print(f"警告: 基因 {gene2} 未在BED文件中找到")
            
            if not genes1_info or not genes2_info:
                print(f"警告: 区块 {i+1} 缺少有效的位置信息，跳过")
                continue
            
            # 计算第1列区块的范围
            chrom1 = genes1_info[0][0]  # 假设同一区块的基因在同一染色体上
            start1 = min(info[1] for info in genes1_info)
            end1 = max(info[2] for info in genes1_info)
            
            # 计算第2列区块的范围
            chrom2 = genes2_info[0][0]
            start2 = min(info[1] for info in genes2_info)
            end2 = max(info[2] for info in genes2_info)
            
            # 写入link文件
            f.write(f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\n")
    
    print(f"共线性link文件已生成: {output_file}")


if __name__ == "__main__":
    # 输入文件路径
    bed_file = "/mnt/user-data/uploads/Cdom.bed"
    anchors_file = "/mnt/user-data/uploads/Cdom_Cdom.anchors"
    output_file = "/mnt/user-data/outputs/synteny_links.txt"
    
    # 如果提供了命令行参数，使用命令行参数
    if len(sys.argv) >= 4:
        bed_file = sys.argv[1]
        anchors_file = sys.argv[2]
        output_file = sys.argv[3]
    
    generate_synteny_links(bed_file, anchors_file, output_file)
