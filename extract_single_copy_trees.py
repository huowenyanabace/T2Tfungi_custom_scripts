#!/usr/bin/env python3
"""
单拷贝直系同源基因树提取与处理脚本

功能：
1. 从基因树文件中提取单拷贝直系同源组的基因树
2. 将基因名称简化为物种名称
3. 过滤掉包含重复物种的基因树
4. 输出符合ASTRAL格式的基因树文件

作者: Claude
日期: 2026-01-06
"""

import re
import sys
from collections import Counter


def load_single_copy_orthogroups(filepath):
    """
    加载单拷贝直系同源组ID列表
    
    参数:
        filepath (str): 单拷贝直系同源组文件路径
    
    返回:
        set: 单拷贝直系同源组ID的集合
    """
    single_copy_ogs = set()
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line:  # 忽略空行
                    single_copy_ogs.add(line)
        print(f"✓ 成功加载 {len(single_copy_ogs)} 个单拷贝直系同源组")
        return single_copy_ogs
    except FileNotFoundError:
        print(f"✗ 错误: 找不到文件 '{filepath}'")
        sys.exit(1)
    except Exception as e:
        print(f"✗ 读取文件时出错: {e}")
        sys.exit(1)


def extract_species_name(gene_name):
    """
    从基因名称中提取物种名称
    
    参数:
        gene_name (str): 完整的基因名称，格式如 "Scere_pep_YNL317W"
    
    返回:
        str: 物种名称，如 "Scere"
    """
    # 提取第一个下划线之前的部分作为物种名称
    match = re.match(r'^([^_]+)_', gene_name)
    if match:
        return match.group(1)
    else:
        # 如果格式不符合预期，返回原始名称
        return gene_name


def simplify_tree_to_species(tree_str):
    """
    将基因树中的基因名称简化为物种名称
    
    参数:
        tree_str (str): 原始基因树字符串（Newick格式）
    
    返回:
        str: 简化后的基因树字符串
    """
    # 使用正则表达式查找所有基因名称（字母开头，包含下划线和字母数字字符）
    # 匹配模式：物种_pep_基因ID 或 物种_其他_其他
    pattern = r'([A-Za-z]+_[a-z]+_[A-Za-z0-9\-]+)'
    
    def replace_with_species(match):
        gene_name = match.group(1)
        species_name = extract_species_name(gene_name)
        return species_name
    
    simplified_tree = re.sub(pattern, replace_with_species, tree_str)
    return simplified_tree


def has_duplicate_species(tree_str):
    """
    检查基因树中是否有重复的物种
    
    参数:
        tree_str (str): 基因树字符串（Newick格式）
    
    返回:
        bool: 如果有重复物种返回True，否则返回False
    """
    # 提取所有物种名称（不包括节点标签如n1, n2等）
    # 物种名称模式：大写字母开头的单词（不包含冒号）
    pattern = r'\b([A-Z][a-z]+)(?=:)'
    species_list = re.findall(pattern, tree_str)
    
    # 使用Counter统计每个物种出现的次数
    species_count = Counter(species_list)
    
    # 检查是否有物种出现超过一次
    for species, count in species_count.items():
        if count > 1:
            return True
    
    return False


def process_gene_trees(input_file, single_copy_ogs, output_file):
    """
    处理基因树文件，提取并简化单拷贝直系同源组的基因树
    
    参数:
        input_file (str): 输入基因树文件路径
        single_copy_ogs (set): 单拷贝直系同源组ID集合
        output_file (str): 输出文件路径
    
    返回:
        dict: 处理统计信息
    """
    stats = {
        'total_trees': 0,
        'single_copy_trees': 0,
        'trees_with_duplicates': 0,
        'output_trees': 0
    }
    
    try:
        with open(input_file, 'r', encoding='utf-8') as infile, \
             open(output_file, 'w', encoding='utf-8') as outfile:
            
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                
                stats['total_trees'] += 1
                
                # 解析OG ID和基因树
                # 格式: OG0000000: (tree_structure);
                match = re.match(r'^(OG\d+):\s*(.+)$', line)
                if not match:
                    print(f"⚠ 警告: 无法解析行格式，跳过: {line[:50]}...")
                    continue
                
                og_id = match.group(1)
                tree_str = match.group(2)
                
                # 检查是否为单拷贝直系同源组
                if og_id not in single_copy_ogs:
                    continue
                
                stats['single_copy_trees'] += 1
                
                # 简化基因树为物种树
                simplified_tree = simplify_tree_to_species(tree_str)
                
                # 检查是否有重复物种
                if has_duplicate_species(simplified_tree):
                    stats['trees_with_duplicates'] += 1
                    continue
                
                # 写入输出文件（不包含OG ID，只保留树结构）
                outfile.write(simplified_tree + '\n')
                stats['output_trees'] += 1
        
        return stats
        
    except FileNotFoundError:
        print(f"✗ 错误: 找不到文件 '{input_file}'")
        sys.exit(1)
    except Exception as e:
        print(f"✗ 处理文件时出错: {e}")
        sys.exit(1)


def main():
    """主函数"""
    print("=" * 70)
    print("单拷贝直系同源基因树提取与处理工具")
    print("=" * 70)
    print()
    
    # 文件路径配置
    single_copy_file = "Orthogroups_SingleCopyOrthologues.txt"
    gene_trees_file = "Resolved_Gene_Trees.txt"
    output_file = "Astral_input.tree"
    
    # 检查命令行参数
    if len(sys.argv) > 1:
        single_copy_file = sys.argv[1]
    if len(sys.argv) > 2:
        gene_trees_file = sys.argv[2]
    if len(sys.argv) > 3:
        output_file = sys.argv[3]
    
    print(f"输入文件:")
    print(f"  - 单拷贝直系同源组: {single_copy_file}")
    print(f"  - 基因树文件: {gene_trees_file}")
    print(f"输出文件:")
    print(f"  - ASTRAL输入文件: {output_file}")
    print()
    
    # 步骤1: 加载单拷贝直系同源组列表
    print("[1/2] 加载单拷贝直系同源组...")
    single_copy_ogs = load_single_copy_orthogroups(single_copy_file)
    print()
    
    # 步骤2: 处理基因树文件
    print("[2/2] 处理基因树文件...")
    stats = process_gene_trees(gene_trees_file, single_copy_ogs, output_file)
    print(f"✓ 处理完成")
    print()
    
    # 显示统计信息
    print("=" * 70)
    print("处理统计信息")
    print("=" * 70)
    print(f"总基因树数量:              {stats['total_trees']:>8}")
    print(f"单拷贝直系同源组数量:      {stats['single_copy_trees']:>8}")
    print(f"包含重复物种的树 (已过滤): {stats['trees_with_duplicates']:>8}")
    print(f"最终输出的基因树数量:      {stats['output_trees']:>8}")
    print("=" * 70)
    print()
    
    if stats['output_trees'] > 0:
        print(f"✓ 成功! 结果已保存到: {output_file}")
        print(f"  该文件可直接用于ASTRAL物种树分析")
    else:
        print("⚠ 警告: 没有生成任何输出基因树")
        print("  请检查输入文件是否正确")
    print()


if __name__ == "__main__":
    main()
