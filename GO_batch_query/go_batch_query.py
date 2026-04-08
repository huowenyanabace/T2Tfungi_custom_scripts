#!/usr/bin/env python3
"""
GO Term 批量查询工具
支持从txt文件读取GO IDs，并提供文本和CSV两种输出格式
"""

from goatools import obo_parser
import csv
import sys
import os
import argparse

class GOBatchQuery:
    def __init__(self, obo_file="go-basic.obo"):
        """初始化GO数据库"""
        print(f"正在加载GO数据库: {obo_file}")
        if not os.path.exists(obo_file):
            print(f"错误: GO数据库文件 '{obo_file}' 不存在")
            print("请先下载: wget http://purl.obolibrary.org/obo/go/go-basic.obo")
            sys.exit(1)
        
        self.go = obo_parser.GODag(obo_file)
        print(f"加载完成！共包含 {len(self.go)} 个GO terms\n")
    
    def read_go_ids(self, filename):
        """从文件中读取GO ID列表"""
        if not os.path.exists(filename):
            print(f"错误: 输入文件 '{filename}' 不存在")
            sys.exit(1)
        
        with open(filename, 'r', encoding='utf-8') as f:
            go_ids = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        
        return go_ids
    
    def query_single_term(self, go_id):
        """查询单个GO term，返回字典"""
        go_id = go_id.strip()
        
        if go_id not in self.go:
            return {
                'id': go_id,
                'name': 'NOT FOUND',
                'namespace': '',
                'definition': '',
                'depth': '',
                'parents': '',
                'children_count': '',
                'found': False
            }
        
        term = self.go[go_id]
        parents = '; '.join([f"{p.id} ({p.name})" for p in term.parents]) if term.parents else ''
        
        return {
            'id': go_id,
            'name': term.name,
            'namespace': term.namespace,
            'definition': term.defn if hasattr(term, 'defn') else '',
            'depth': term.depth,
            'parents': parents,
            'children_count': len(term.children) if term.children else 0,
            'found': True
        }
    
    def output_text(self, go_ids, output_file=None):
        """以文本格式输出结果"""
        file_handle = None
        if output_file:
            file_handle = open(output_file, 'w', encoding='utf-8')
            print(f"结果将保存到: {output_file}\n")
        
        for i, go_id in enumerate(go_ids, 1):
            result = self.query_single_term(go_id)
            
            output = []
            output.append(f"\n{'='*70}")
            output.append(f"[{i}/{len(go_ids)}] GO ID: {result['id']}")
            output.append('='*70)
            
            if not result['found']:
                output.append("状态: 未找到此GO ID")
            else:
                output.append(f"名称: {result['name']}")
                output.append(f"命名空间: {result['namespace']}")
                output.append(f"定义: {result['definition']}")
                output.append(f"层级深度: {result['depth']}")
                output.append(f"父节点: {result['parents'] if result['parents'] else '无（根节点）'}")
                output.append(f"子节点数: {result['children_count']}")
            
            output.append('='*70)
            
            text = '\n'.join(output)
            print(text)
            
            if file_handle:
                file_handle.write(text + '\n')
        
        if file_handle:
            file_handle.close()
            print(f"\n结果已保存到: {output_file}")
        
        print(f"\n查询完成！共处理 {len(go_ids)} 个GO ID")
    
    def output_csv(self, go_ids, output_file):
        """以CSV格式输出结果"""
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['GO_ID', 'Name', 'Namespace', 'Definition', 'Depth', 'Parents', 'Children_Count', 'Status']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            
            for i, go_id in enumerate(go_ids, 1):
                result = self.query_single_term(go_id)
                
                writer.writerow({
                    'GO_ID': result['id'],
                    'Name': result['name'],
                    'Namespace': result['namespace'],
                    'Definition': result['definition'],
                    'Depth': result['depth'],
                    'Parents': result['parents'],
                    'Children_Count': result['children_count'],
                    'Status': 'Found' if result['found'] else 'Not Found'
                })
                
                print(f"[{i}/{len(go_ids)}] 已处理: {go_id}")
        
        print(f"\n结果已保存到: {output_file}")
        print(f"查询完成！共处理 {len(go_ids)} 个GO ID")

def main():
    parser = argparse.ArgumentParser(
        description='GO Term 批量查询工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 从go_ids.txt读取，结果输出到屏幕
  python3 go_batch_query.py go_ids.txt
  
  # 从go_ids.txt读取，保存为文本格式
  python3 go_batch_query.py go_ids.txt -o results.txt
  
  # 从go_ids.txt读取，保存为CSV格式
  python3 go_batch_query.py go_ids.txt -o results.csv -f csv
  
  # 指定自定义的OBO文件路径
  python3 go_batch_query.py go_ids.txt --obo /path/to/go.obo
        """
    )
    
    parser.add_argument('input_file', help='包含GO IDs的txt文件（每行一个GO ID）')
    parser.add_argument('-o', '--output', help='输出文件路径（不指定则只输出到屏幕）')
    parser.add_argument('-f', '--format', choices=['text', 'csv'], default='text',
                        help='输出格式: text 或 csv (默认: text)')
    parser.add_argument('--obo', default='go-basic.obo',
                        help='GO OBO数据库文件路径 (默认: go-basic.obo)')
    
    args = parser.parse_args()
    
    # 如果输出格式是CSV但没有指定输出文件，自动生成文件名
    if args.format == 'csv' and not args.output:
        args.output = args.input_file.rsplit('.', 1)[0] + '_results.csv'
        print(f"CSV格式需要输出文件，将使用: {args.output}\n")
    
    # 初始化查询工具
    query_tool = GOBatchQuery(args.obo)
    
    # 读取GO IDs
    go_ids = query_tool.read_go_ids(args.input_file)
    print(f"从 {args.input_file} 读取到 {len(go_ids)} 个GO ID\n")
    
    # 执行查询并输出
    if args.format == 'csv':
        query_tool.output_csv(go_ids, args.output)
    else:
        query_tool.output_text(go_ids, args.output)

if __name__ == "__main__":
    main()
