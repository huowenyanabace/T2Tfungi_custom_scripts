#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
从*.pep.fa文件中提取最长蛋白序列的脚本
作者：Qwen
日期：2025-12-23
"""

import os
import re
from collections import defaultdict

def parse_fasta(fasta_file):
    """
    解析FASTA文件，返回序列字典
    键为蛋白质ID，值为(序列, 基因ID)的元组
    """
    sequences = {}
    current_seq_id = None
    current_gene_id = None
    current_seq = []

    with open(fasta_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # 如果当前有正在处理的序列，先保存
                if current_seq_id:
                    sequences[current_seq_id] = (''.join(current_seq), current_gene_id)

                # 解析新的序列头
                header = line[1:]  # 去掉'>'
                # 提取蛋白质ID（通常是第一个空格前的部分）
                parts = header.split()
                current_seq_id = parts[0] if parts else 'unknown_protein'

                # 提取基因ID（格式为 gene:GENEID）
                gene_match = re.search(r'gene:(\S+)', header)
                # 仅获取基因ID部分，如果匹配失败则使用'unknown_gene'
                current_gene_id = gene_match.group(1).split()[0] if gene_match else 'unknown_gene'

                current_seq = []
            else:
                current_seq.append(line)

    # 保存最后一个序列
    if current_seq_id:
        sequences[current_seq_id] = (''.join(current_seq), current_gene_id)

    return sequences

def get_longest_sequences(sequences):
    """
    根据基因ID获取最长的蛋白质序列
    返回包含最长序列信息的字典
    """
    gene_to_proteins = defaultdict(list)
    
    # 按基因ID分组蛋白质
    for protein_id, (sequence, gene_id) in sequences.items():
        gene_to_proteins[gene_id].append((protein_id, sequence, len(sequence)))
    
    longest_sequences = {}
    
    # 为每个基因选择最长的蛋白质序列
    for gene_id, protein_list in gene_to_proteins.items():
        # 按序列长度排序，取最长的
        longest_protein = max(protein_list, key=lambda x: x[2])
        protein_id, sequence, length = longest_protein
        longest_sequences[protein_id] = (sequence, gene_id, length)
    
    return longest_sequences

def write_results(longest_sequences, output_fasta, output_txt):
    """
    将结果写入FASTA文件和文本文件
    """
    # 写入FASTA文件
    with open(output_fasta, 'w', encoding='utf-8') as fa:
        for protein_id, (sequence, gene_id, length) in longest_sequences.items():
            fa.write(f'>{protein_id};gene={gene_id};length={length}\n')
            # 按每行60个字符格式化序列
            for i in range(0, len(sequence), 60):
                fa.write(sequence[i:i+60] + '\n')
    
    # 写入文本文件
    with open(output_txt, 'w', encoding='utf-8') as txt:
        txt.write('蛋白质ID\t基因ID\t序列长度\n')
        for protein_id, (sequence, gene_id, length) in longest_sequences.items():
            txt.write(f'{protein_id}\t{gene_id}\t{length}\n')

def main():
    """
    主函数
    """
    print('开始处理.pep.fa文件...')

    # 创建输出目录
    output_dir = 'longest_pep_results'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f'创建输出目录: {output_dir}')

    # 获取当前目录下所有的.pep.fa文件，但排除已生成的.longest.pep.fa文件
    pep_files = [f for f in os.listdir('.') if f.endswith('.pep.fa') and not f.endswith('.longest.pep.fa')]

    if not pep_files:
        print('当前目录下没有找到.pep.fa文件')
        return

    print(f'找到以下.pep.fa文件: {pep_files}')

    for pep_file in pep_files:
        print(f'\n正在处理文件: {pep_file}')

        # 解析FASTA文件
        sequences = parse_fasta(pep_file)
        print(f'共解析到 {len(sequences)} 个蛋白质序列')

        # 获取最长序列
        longest_sequences = get_longest_sequences(sequences)
        print(f'为 {len(longest_sequences)} 个基因选择了最长的蛋白质序列')

        # 生成输出文件名
        base_name = pep_file.replace('.pep.fa', '')
        output_fasta = os.path.join(output_dir, f'{base_name}.longest.pep.fa')
        output_txt = os.path.join(output_dir, f'{base_name}.longest.info.txt')

        # 写入结果
        write_results(longest_sequences, output_fasta, output_txt)

        print(f'结果已保存到:')
        print(f'  - {output_fasta} (最长蛋白质序列)')
        print(f'  - {output_txt} (蛋白质ID、基因ID和长度信息)')

    print('\n处理完成!')

if __name__ == '__main__':
    main()