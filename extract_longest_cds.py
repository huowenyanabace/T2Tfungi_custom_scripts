#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
从 *.cds.fa 文件中提取最长 CDS 序列的脚本

功能：
1. 从 *.cds.fa 中提取 CDS 的 ID 和基因 ID，以及序列的长度，并将它们的信息保存在一个 .txt 文件中
2. 如果一个基因 ID 对应 2 个或 2 个以上的 CDS ID，比较这几个 CDS 的长度，只保留长度最长的那个 ID
3. 最后根据 CDS ID 抽取最终的 cds.fa
结果输出到一个专门的文件夹，全程用中文
"""

import os
import re
from collections import defaultdict
import glob


def parse_fasta_file(fasta_file):
    """
    解析 FASTA 文件，提取 CDS ID、基因 ID 和序列
    
    Args:
        fasta_file (str): FASTA 文件路径
    
    Returns:
        list: 包含 (cds_id, gene_id, sequence) 的列表
    """
    cds_data = []
    
    with open(fasta_file, 'r', encoding='utf-8') as f:
        header = None
        sequence = ""
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # 如果已经有上一条记录，保存它
                if header is not None and sequence:
                    # 提取 CDS ID 和基因 ID
                    cds_id, gene_id = extract_ids_from_header(header)
                    cds_data.append((cds_id, gene_id, sequence))
                
                header = line[1:]  # 去掉开头的 '>'
                sequence = ""
            else:
                sequence += line
        
        # 处理最后一条记录
        if header is not None and sequence:
            cds_id, gene_id = extract_ids_from_header(header)
            cds_data.append((cds_id, gene_id, sequence))
    
    return cds_data


def extract_ids_from_header(header):
    """
    从 FASTA 头部提取 CDS ID 和基因 ID
    
    Args:
        header (str): FASTA 头部字符串
    
    Returns:
        tuple: (cds_id, gene_id)
    """
    # 示例头部: EKM83519 cds supercontig:Agabi_varbur_1:JH971385:518309:519915:-1 gene:AGABI1DRAFT_32645 gene_biotype:protein_coding transcript_biotype:protein_coding description:hypothetical protein
    
    # 提取 CDS ID (第一个空格前的部分)
    parts = header.split(' ')
    cds_id = parts[0]
    
    # 提取基因 ID (包含 gene: 的部分)
    gene_id = None
    for part in parts:
        if part.startswith('gene:'):
            gene_id = part.replace('gene:', '')
            break
    
    # 如果没有找到 gene: 格式的基因 ID，尝试其他格式
    if gene_id is None:
        # 尝试匹配 gene= 或其他格式
        gene_match = re.search(r'gene[:=]([A-Za-z0-9_]+)', header)
        if gene_match:
            gene_id = gene_match.group(1)
    
    return cds_id, gene_id


def write_info_to_txt(cds_data, output_txt_file):
    """
    将提取的信息保存到 .txt 文件中
    
    Args:
        cds_data (list): 包含 (cds_id, gene_id, sequence) 的列表
        output_txt_file (str): 输出的 .txt 文件路径
    """
    with open(output_txt_file, 'w', encoding='utf-8') as f:
        f.write("CDS ID\t基因 ID\t序列长度\n")
        for cds_id, gene_id, sequence in cds_data:
            seq_length = len(sequence)
            f.write(f"{cds_id}\t{gene_id}\t{seq_length}\n")


def find_longest_cds_per_gene(cds_data):
    """
    对于一个基因 ID 对应多个 CDS ID 的情况，比较长度并保留最长的 CDS ID
    
    Args:
        cds_data (list): 包含 (cds_id, gene_id, sequence) 的列表
    
    Returns:
        list: 包含 (cds_id, gene_id, sequence) 的列表，每个基因只保留最长的 CDS
    """
    # 按基因 ID 分组
    gene_to_cdss = defaultdict(list)
    for cds_id, gene_id, sequence in cds_data:
        if gene_id:  # 只处理有基因 ID 的记录
            gene_to_cdss[gene_id].append((cds_id, sequence))
    
    # 对每个基因，保留最长的 CDS
    longest_cdss = []
    for gene_id, cds_list in gene_to_cdss.items():
        if len(cds_list) == 1:
            # 如果只有一个 CDS，直接保留
            cds_id, sequence = cds_list[0]
            longest_cdss.append((cds_id, gene_id, sequence))
        else:
            # 如果有多个 CDS，保留最长的
            longest_cds = max(cds_list, key=lambda x: len(x[1]))
            longest_cdss.append((longest_cds[0], gene_id, longest_cds[1]))
    
    # 添加没有基因 ID 的 CDS（如果有的话）
    for cds_id, gene_id, sequence in cds_data:
        if not gene_id:
            longest_cdss.append((cds_id, gene_id, sequence))
    
    return longest_cdss


def write_filtered_fasta(longest_cdss, output_fasta_file):
    """
    将筛选后的 CDS 序列写入 FASTA 文件
    
    Args:
        longest_cdss (list): 包含 (cds_id, gene_id, sequence) 的列表
        output_fasta_file (str): 输出的 FASTA 文件路径
    """
    with open(output_fasta_file, 'w', encoding='utf-8') as f:
        for cds_id, gene_id, sequence in longest_cdss:
            if gene_id:
                f.write(f">{cds_id} gene:{gene_id}\n")
            else:
                f.write(f">{cds_id}\n")
            # 每行最多写 80 个字符
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")


def main():
    """
    主函数
    """
    print("开始处理 CDS 文件...")
    
    # 创建输出目录
    output_dir = "最长CDS结果"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 获取所有 .cds.fa 文件
    cds_files = glob.glob("*.cds.fa")
    
    if not cds_files:
        print("在当前目录中未找到 .cds.fa 文件")
        return
    
    print(f"找到 {len(cds_files)} 个 CDS 文件: {cds_files}")
    
    for cds_file in cds_files:
        print(f"正在处理 {cds_file}...")
        
        # 获取文件名（不含扩展名）
        base_name = os.path.splitext(cds_file)[0]
        
        # 解析 FASTA 文件
        cds_data = parse_fasta_file(cds_file)
        print(f"共提取到 {len(cds_data)} 条 CDS 记录")
        
        # 将信息保存到 .txt 文件
        txt_file = os.path.join(output_dir, f"{base_name}_cds_info.txt")
        write_info_to_txt(cds_data, txt_file)
        print(f"CDS 信息已保存到 {txt_file}")
        
        # 找到每个基因的最长 CDS
        longest_cdss = find_longest_cds_per_gene(cds_data)
        print(f"筛选后剩余 {len(longest_cdss)} 条最长 CDS 记录")
        
        # 将筛选后的 CDS 写入新的 FASTA 文件
        filtered_fasta = os.path.join(output_dir, f"{base_name}_longest_cds.fa")
        write_filtered_fasta(longest_cdss, filtered_fasta)
        print(f"最长 CDS 序列已保存到 {filtered_fasta}")
        
        print("-" * 50)
    
    print("所有文件处理完成！结果已保存到 '最长CDS结果' 目录中。")


if __name__ == "__main__":
    main()