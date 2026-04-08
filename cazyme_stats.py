#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CAZyme注释结果统计脚本

功能：
1. 遍历工作目录下所有*_result文件夹中的overview.txt文件
2. 统计各物种的CAZyme大类（GH、GT、PL、CE、AA、CBM）数量
3. 统计各小类（如GH1、GH2等）的数量
4. 输出统计结果到TSV文件

作者：Claude
日期：2025
"""

import argparse
import os
import re
from collections import defaultdict
from glob import glob


def parse_cazyme_families(cell_value):
    """
    从注释结果单元格中解析CAZyme家族
    例如: "GH18(76-461)" -> ["GH18"]
          "GT2_Chitin_synth_2(629-1151)" -> ["GT2"]
          "CBM48+GH13_8" -> ["CBM48", "GH13"]
          "GH13_40+3.2.1.20" -> ["GH13"]
    """
    if cell_value == '-' or not cell_value:
        return []
    
    families = []
    # 按+分割处理多个注释
    parts = cell_value.split('+')
    
    for part in parts:
        # 移除括号及其内容 (位置信息或数字)
        part = re.sub(r'\([^)]*\)', '', part)
        part = part.strip()
        
        # 跳过EC号（如3.2.1.20）
        if re.match(r'^\d+\.\d+\.\d+', part):
            continue
        
        # 提取CAZyme家族名（如GH13_40 -> GH13, GT2_Chitin_synth_2 -> GT2）
        # 匹配模式: AA/GH/GT/CE/PL/CBM + 数字
        match = re.match(r'^((?:AA|GH|GT|CE|PL|CBM)\d+)', part)
        if match:
            families.append(match.group(1))
    
    return families


def get_major_class(family):
    """
    获取CAZyme家族的大类
    例如: GH13 -> GH, GT2 -> GT
    """
    match = re.match(r'^(AA|GH|GT|CE|PL|CBM)', family)
    if match:
        return match.group(1)
    return None


def parse_overview_file(filepath, min_tools=2):
    """
    解析overview.txt文件，统计CAZyme家族
    
    参数:
        filepath: overview.txt文件路径
        min_tools: 最少需要多少个工具支持（默认2）
    
    返回:
        major_counts: 大类统计 {类别: 数量}
        minor_counts: 小类统计 {家族: 数量}
        gene_cazymes: 每个基因的CAZyme注释 {gene_id: [families]}
    """
    major_counts = defaultdict(int)
    minor_counts = defaultdict(int)
    gene_cazymes = {}
    
    with open(filepath, 'r') as f:
        header = f.readline()  # 跳过表头
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) < 5:
                continue
            
            gene_id = parts[0]
            hmmer = parts[1] if len(parts) > 1 else '-'
            hotpep = parts[2] if len(parts) > 2 else '-'
            diamond = parts[3] if len(parts) > 3 else '-'
            num_tools = int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 0
            
            # 根据工具数量过滤
            if num_tools < min_tools:
                continue
            
            # 从三个工具中收集所有CAZyme家族
            all_families = set()
            for cell in [hmmer, hotpep, diamond]:
                families = parse_cazyme_families(cell)
                all_families.update(families)
            
            if all_families:
                gene_cazymes[gene_id] = list(all_families)
                
                for family in all_families:
                    minor_counts[family] += 1
                    major_class = get_major_class(family)
                    if major_class:
                        major_counts[major_class] += 1
    
    return dict(major_counts), dict(minor_counts), gene_cazymes


def find_overview_files(work_dir, pattern='*_result'):
    """
    查找工作目录下所有符合条件的overview.txt文件
    
    返回:
        dict: {物种名: overview.txt路径}
    """
    overview_files = {}
    
    # 查找所有*_result文件夹
    result_dirs = glob(os.path.join(work_dir, pattern))
    
    for result_dir in result_dirs:
        if os.path.isdir(result_dir):
            overview_path = os.path.join(result_dir, 'overview.txt')
            if os.path.exists(overview_path):
                # 从文件夹名提取物种名（去掉_result后缀）
                dir_name = os.path.basename(result_dir)
                species_name = dir_name.replace('_result', '')
                overview_files[species_name] = overview_path
    
    return overview_files


def generate_stats_report(work_dir, output_prefix, min_tools=2, pattern='*_result'):
    """
    生成CAZyme统计报告
    
    参数:
        work_dir: 工作目录
        output_prefix: 输出文件前缀
        min_tools: 最少工具支持数
        pattern: 结果文件夹匹配模式
    """
    # 查找所有overview.txt文件
    overview_files = find_overview_files(work_dir, pattern)
    
    if not overview_files:
        print(f"Error: No result folders matching '{pattern}' found in {work_dir}")
        return
    
    print(f"Found {len(overview_files)} species overview.txt files:")
    for species, path in sorted(overview_files.items()):
        print(f"  - {species}: {path}")
    
    # 统计每个物种的CAZyme
    all_major_counts = {}  # {species: {class: count}}
    all_minor_counts = {}  # {species: {family: count}}
    all_minor_families = set()  # 所有出现过的小类
    
    for species, filepath in sorted(overview_files.items()):
        major_counts, minor_counts, _ = parse_overview_file(filepath, min_tools)
        all_major_counts[species] = major_counts
        all_minor_counts[species] = minor_counts
        all_minor_families.update(minor_counts.keys())
    
    species_list = sorted(all_major_counts.keys())
    
    # ========== 输出大类统计表格 ==========
    major_output = output_prefix + '_major_classes.tsv'
    
    major_classes = ['GH', 'GT', 'PL', 'CE', 'AA', 'CBM']
    major_names = {
        'GH': 'Glycoside Hydrolases (GH)',
        'GT': 'Glycosyl Transferases (GT)',
        'PL': 'Polysaccharide Lyases (PL)',
        'CE': 'Carbohydrate Esterases (CE)',
        'AA': 'Auxiliary Activities (AA)',
        'CBM': 'Carbohydrate-Binding Modules (CBM)'
    }
    
    with open(major_output, 'w') as f:
        # 表头
        header = ['CAZyme Class'] + species_list
        f.write('\t'.join(header) + '\n')
        
        # 各大类统计
        for major_class in major_classes:
            row = [major_names[major_class]]
            for species in species_list:
                count = all_major_counts[species].get(major_class, 0)
                row.append(str(count))
            f.write('\t'.join(row) + '\n')
        
        # 总计
        row = ['Total CAZymes']
        for species in species_list:
            total = sum(all_major_counts[species].values())
            row.append(str(total))
        f.write('\t'.join(row) + '\n')
    
    print(f"\nMajor class statistics saved to: {major_output}")
    
    # ========== 输出小类统计表格 ==========
    minor_output = output_prefix + '_subfamilies.tsv'
    
    # 按大类和数字排序小类
    def sort_key(family):
        match = re.match(r'^(AA|GH|GT|CE|PL|CBM)(\d+)', family)
        if match:
            return (match.group(1), int(match.group(2)))
        return (family, 0)
    
    sorted_families = sorted(all_minor_families, key=sort_key)
    
    with open(minor_output, 'w') as f:
        # 表头
        header = ['CAZyme Family'] + species_list
        f.write('\t'.join(header) + '\n')
        
        # 各小类统计
        for family in sorted_families:
            row = [family]
            for species in species_list:
                count = all_minor_counts[species].get(family, 0)
                row.append(str(count))
            f.write('\t'.join(row) + '\n')
    
    print(f"Subfamily statistics saved to: {minor_output}")
    
    # ========== 打印摘要 ==========
    print("\n" + "=" * 60)
    print("CAZyme Statistics Summary")
    print("=" * 60)
    print(f"{'Species':<15}{'GH':<8}{'GT':<8}{'PL':<8}{'CE':<8}{'AA':<8}{'CBM':<8}{'Total':<8}")
    print("-" * 60)
    for species in species_list:
        counts = all_major_counts[species]
        total = sum(counts.values())
        print(f"{species:<15}{counts.get('GH', 0):<8}{counts.get('GT', 0):<8}"
              f"{counts.get('PL', 0):<8}{counts.get('CE', 0):<8}"
              f"{counts.get('AA', 0):<8}{counts.get('CBM', 0):<8}{total:<8}")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description='CAZyme注释结果统计脚本',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 统计当前目录下所有*_result文件夹
  python cazyme_stats.py -d . -o cazyme_stats
  
  # 指定最少2个工具支持
  python cazyme_stats.py -d /path/to/work -o cazyme_stats -m 2
  
  # 自定义文件夹匹配模式
  python cazyme_stats.py -d . -o cazyme_stats -p "*_dbcan"
        """
    )
    
    parser.add_argument(
        '-d', '--dir',
        default='.',
        help='Working directory path (default: current directory)'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='cazyme_stats',
        help='Output file prefix (default: cazyme_stats)'
    )
    
    parser.add_argument(
        '-m', '--min-tools',
        type=int,
        default=2,
        help='Minimum number of tools required to support an annotation (default: 2)'
    )
    
    parser.add_argument(
        '-p', '--pattern',
        default='*_result',
        help='Result folder matching pattern (default: *_result)'
    )
    
    args = parser.parse_args()
    
    print(f"Working directory: {args.dir}")
    print(f"Folder matching pattern: {args.pattern}")
    print(f"Minimum tools required: {args.min_tools}")
    print()
    
    generate_stats_report(
        work_dir=args.dir,
        output_prefix=args.output,
        min_tools=args.min_tools,
        pattern=args.pattern
    )


if __name__ == '__main__':
    main()
