#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
从 OrthoFinder 结果中提取指定 Orthogroup 的指定物种基因

用法:
    python extract_genes_from_orthogroups.py -i ID.txt -o Orthogroups.tsv -c Cdom.pep -out genes.txt
"""

import argparse
import sys
import csv

def read_orthogroup_ids(id_file):
    """
    读取 Orthogroup ID 列表
    
    Args:
        id_file: ID文件路径
        
    Returns:
        set: Orthogroup ID 集合
    """
    og_ids = set()
    try:
        with open(id_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    og_ids.add(line)
        print(f"✓ 读取到 {len(og_ids)} 个 Orthogroup ID")
        print(f"  ID列表: {', '.join(sorted(og_ids))}")
        return og_ids
    except FileNotFoundError:
        print(f"✗ 错误: 找不到文件 {id_file}")
        sys.exit(1)
    except Exception as e:
        print(f"✗ 读取文件时出错: {e}")
        sys.exit(1)

def extract_genes_from_orthogroups(orthogroups_file, og_ids, target_column, output_file):
    """
    从 Orthogroups.tsv 中提取指定列的基因
    
    Args:
        orthogroups_file: Orthogroups.tsv 文件路径
        og_ids: 要提取的 Orthogroup ID 集合
        target_column: 目标列名（如 Cdom.pep）
        output_file: 输出文件路径
    """
    try:
        # 读取文件并查找目标列
        with open(orthogroups_file, 'r', newline='') as f:
            # 处理可能存在的 \r\n 换行符
            reader = csv.DictReader(f, delimiter='\t')
            
            # 检查目标列是否存在
            if target_column not in reader.fieldnames:
                print(f"\n✗ 错误: 列 '{target_column}' 不存在!")
                print(f"  可用的列: {', '.join(reader.fieldnames)}")
                sys.exit(1)
            
            print(f"\n✓ 找到目标列: {target_column}")
            
            # 提取基因
            all_genes = []
            found_ogs = set()
            og_gene_counts = {}
            
            for row in reader:
                og_id = row['Orthogroup'].strip()
                
                if og_id in og_ids:
                    found_ogs.add(og_id)
                    gene_field = row[target_column]
                    
                    # 解析基因列表（处理引号和逗号分隔）
                    genes = parse_gene_field(gene_field)
                    
                    if genes:
                        all_genes.extend(genes)
                        og_gene_counts[og_id] = len(genes)
                        print(f"  {og_id}: 找到 {len(genes)} 个基因")
                    else:
                        og_gene_counts[og_id] = 0
                        print(f"  {og_id}: 该列为空")
            
            # 检查是否所有ID都找到了
            not_found = og_ids - found_ogs
            if not_found:
                print(f"\n⚠️  警告: 以下 Orthogroup ID 未找到:")
                for og in sorted(not_found):
                    print(f"    - {og}")
            
            # 去重（保持顺序）
            unique_genes = []
            seen = set()
            for gene in all_genes:
                if gene not in seen:
                    unique_genes.append(gene)
                    seen.add(gene)
            
            # 统计信息
            print(f"\n提取统计:")
            print(f"  查询的 Orthogroup 数: {len(og_ids)}")
            print(f"  找到的 Orthogroup 数: {len(found_ogs)}")
            print(f"  提取的基因总数: {len(all_genes)}")
            print(f"  去重后基因数: {len(unique_genes)}")
            
            if len(all_genes) != len(unique_genes):
                print(f"  (移除了 {len(all_genes) - len(unique_genes)} 个重复基因)")
            
            # 输出到文件
            with open(output_file, 'w') as out:
                for gene in unique_genes:
                    out.write(f"{gene}\n")
            
            print(f"\n✓ 结果已保存到: {output_file}")
            
            # 显示每个Orthogroup的基因数
            if og_gene_counts:
                print(f"\n各 Orthogroup 的基因数:")
                for og in sorted(og_gene_counts.keys()):
                    print(f"  {og}: {og_gene_counts[og]} 个基因")
            
            # 显示前10个基因作为预览
            if unique_genes:
                print(f"\n基因列表预览（前10个）:")
                for gene in unique_genes[:10]:
                    print(f"  {gene}")
                if len(unique_genes) > 10:
                    print(f"  ...")
            
            return unique_genes
            
    except FileNotFoundError:
        print(f"✗ 错误: 找不到文件 {orthogroups_file}")
        sys.exit(1)
    except Exception as e:
        print(f"✗ 处理文件时出错: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

def parse_gene_field(field):
    """
    解析基因字段，处理引号和逗号分隔
    
    Args:
        field: 基因字段字符串
        
    Returns:
        list: 基因ID列表
    """
    if not field or field.strip() == '' or field.strip() == '-':
        return []
    
    # 去除外层引号
    field = field.strip()
    if field.startswith('"') and field.endswith('"'):
        field = field[1:-1]
    
    # 按逗号分隔并清理空格
    genes = [gene.strip() for gene in field.split(',')]
    
    # 过滤空字符串
    genes = [gene for gene in genes if gene]
    
    return genes

def main():
    """
    主函数
    """
    parser = argparse.ArgumentParser(
        description='从 OrthoFinder 结果中提取指定 Orthogroup 的指定物种基因',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 基本用法
  python extract_genes_from_orthogroups.py -i ID.txt -o Orthogroups.tsv -c Cdom.pep -out genes.txt
  
  # 提取不同物种的基因
  python extract_genes_from_orthogroups.py -i ID.txt -o Orthogroups.tsv -c Scere.pep -out yeast_genes.txt
  
  # 指定完整路径
  python extract_genes_from_orthogroups.py \\
    -i /path/to/ID.txt \\
    -o /path/to/Orthogroups.tsv \\
    -c Cdom.pep \\
    -out expanded_genes.txt

输出文件格式:
  每行一个基因ID，例如:
  COD001339-RA
  COD001732-RA
  COD002939-RA
  ...
        """
    )
    
    parser.add_argument('-i', '--id-file',
                        required=True,
                        help='Orthogroup ID 列表文件 (如 ID.txt)')
    
    parser.add_argument('-o', '--orthogroups',
                        required=True,
                        help='OrthoFinder 输出的 Orthogroups.tsv 文件')
    
    parser.add_argument('-c', '--column',
                        required=True,
                        help='目标物种列名 (如 Cdom.pep, Scere.pep)')
    
    parser.add_argument('-out', '--output',
                        default='extracted_genes.txt',
                        help='输出文件名 (默认: extracted_genes.txt)')
    
    parser.add_argument('--show-columns',
                        action='store_true',
                        help='显示 Orthogroups.tsv 中的所有列名然后退出')
    
    args = parser.parse_args()
    
    # 如果只是查看列名
    if args.show_columns:
        try:
            with open(args.orthogroups, 'r') as f:
                header = f.readline().strip().split('\t')
                print("Orthogroups.tsv 中的列:")
                for i, col in enumerate(header, 1):
                    print(f"  {i}. {col}")
            sys.exit(0)
        except FileNotFoundError:
            print(f"✗ 错误: 找不到文件 {args.orthogroups}")
            sys.exit(1)
    
    # 打印参数信息
    print("参数设置:")
    print(f"  Orthogroup ID 文件: {args.id_file}")
    print(f"  Orthogroups 文件: {args.orthogroups}")
    print(f"  目标列: {args.column}")
    print(f"  输出文件: {args.output}")
    print()
    
    # 执行提取
    og_ids = read_orthogroup_ids(args.id_file)
    genes = extract_genes_from_orthogroups(
        args.orthogroups, 
        og_ids, 
        args.column, 
        args.output
    )
    
    print("\n" + "="*50)
    print("✓ 提取完成！")
    print("="*50)

if __name__ == '__main__':
    main()
