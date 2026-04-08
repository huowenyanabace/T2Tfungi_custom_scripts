# ============================================================
# eggNOG-mapper 注释结果富集分析（修复版）
# 修复了 "more columns than column names" 错误
# ============================================================

library(clusterProfiler)
library(tidyr)
library(dplyr)
library(ggplot2)
library(enrichplot)

cat("========================================\n")
cat("eggNOG-mapper 富集分析开始\n")
cat("========================================\n\n")

# ============================================================
# 第一步：读取 eggNOG-mapper 注释文件（改进版）
# ============================================================

# 修改为你的 eggNOG 文件路径
eggnog_file <- "20210058.emapper.annotations"  # 👈 修改这里

cat("正在读取文件:", eggnog_file, "\n")

# 改进的读取方法 - 分步处理
tryCatch({
  
  # 第1步：读取所有行
  all_lines <- readLines(eggnog_file)
  cat("  读取了", length(all_lines), "行\n")
  
  # 第2步：找到表头行（以 #query 或 query 开头）
  header_line <- NULL
  header_idx <- 0
  
  for(i in 1:min(50, length(all_lines))) {  # 只检查前50行
    line <- all_lines[i]
    if(grepl("^#query", line) || (grepl("^query", line) && grepl("\t", line))) {
      header_line <- line
      header_idx <- i
      break
    }
  }
  
  if(is.null(header_line)) {
    stop("未找到表头行！请检查文件格式。")
  }
  
  cat("  找到表头行（第", header_idx, "行）\n")
  
  # 第3步：解析表头
  header_line <- gsub("^#", "", header_line)  # 去掉开头的 #
  col_names <- unlist(strsplit(header_line, "\t"))
  cat("  列数:", length(col_names), "\n")
  cat("  列名:", paste(head(col_names, 5), collapse = ", "), "...\n")
  
  # 第4步：读取数据行（跳过注释和表头）
  data_lines <- all_lines[(header_idx + 1):length(all_lines)]
  
  # 移除空行和注释行
  data_lines <- data_lines[!grepl("^#", data_lines)]
  data_lines <- data_lines[data_lines != ""]
  
  cat("  数据行数:", length(data_lines), "\n")
  
  # 第5步：检查并修复列数不一致的问题
  cat("  检查数据格式...\n")
  
  # 统计每行的列数
  col_counts <- sapply(data_lines, function(x) length(strsplit(x, "\t")[[1]]))
  expected_cols <- length(col_names)
  
  if(any(col_counts != expected_cols)) {
    n_mismatch <- sum(col_counts != expected_cols)
    cat("  ⚠️  发现", n_mismatch, "行的列数与表头不一致\n")
    cat("  正在修复...\n")
    
    # 修复方法：对于列数不足的行，补充空值；对于列数过多的行，截断
    data_lines_fixed <- sapply(data_lines, function(line) {
      cols <- strsplit(line, "\t")[[1]]
      if(length(cols) < expected_cols) {
        # 补充空值
        cols <- c(cols, rep("", expected_cols - length(cols)))
      } else if(length(cols) > expected_cols) {
        # 截断多余的列
        cols <- cols[1:expected_cols]
      }
      paste(cols, collapse = "\t")
    })
    
    data_lines <- data_lines_fixed
    cat("  ✓ 修复完成\n")
  }
  
  # 第6步：组合成数据框
  cat("  构建数据框...\n")
  
  # 写入临时文件
  temp_file <- tempfile(fileext = ".tsv")
  writeLines(c(paste(col_names, collapse = "\t"), data_lines), temp_file)
  
  # 读取数据框
  eggnog <- read.table(temp_file, 
                      sep = "\t", 
                      header = TRUE, 
                      quote = "",
                      stringsAsFactors = FALSE,
                      fill = TRUE,
                      comment.char = "")
  
  # 删除临时文件
  unlink(temp_file)
  
  cat("✓ 文件读取成功！\n")
  cat("  总基因数：", nrow(eggnog), "\n")
  cat("  列数：", ncol(eggnog), "\n\n")
  
  # 显示前几行
  cat("数据预览（前3行，前6列）：\n")
  print(head(eggnog[, 1:min(6, ncol(eggnog))], 3))
  cat("\n")
  
}, error = function(e) {
  cat("✗ 读取文件失败！\n")
  cat("错误信息:", e$message, "\n\n")
  cat("请尝试以下解决方案：\n")
  cat("1. 检查文件路径是否正确\n")
  cat("2. 确认文件是 eggNOG-mapper 的标准输出格式\n")
  cat("3. 检查文件是否损坏\n")
  cat("4. 尝试用文本编辑器打开文件查看前几行\n\n")
  stop(e)
})

# ============================================================
# 第二步：提取 GO 注释
# ============================================================

cat("========================================\n")
cat("提取 GO 注释\n")
cat("========================================\n")

# 确保列名标准化
if(!"query" %in% colnames(eggnog)) {
  if(ncol(eggnog) > 0) {
    colnames(eggnog)[1] <- "query"
  } else {
    stop("数据框为空！")
  }
}

# 查找 GO 列
go_col <- NULL
for(col in c("GOs", "GO", "GO_terms", "go")) {
  if(col %in% colnames(eggnog)) {
    go_col <- col
    break
  }
}

if(is.null(go_col)) {
  cat("⚠️  未找到 GO 注释列\n")
  cat("可用的列名：", paste(colnames(eggnog), collapse = ", "), "\n\n")
  term2gene_go <- NULL
} else {
  cat("找到 GO 注释列：", go_col, "\n")
  
  # 提取有 GO 注释的基因
  go_data <- eggnog[eggnog[[go_col]] != "-" & 
                    eggnog[[go_col]] != "" & 
                    !is.na(eggnog[[go_col]]), 
                    c("query", go_col)]
  
  cat("有 GO 注释的基因数：", nrow(go_data), "\n")
  
  if(nrow(go_data) > 0) {
    # 拆分多个 GO terms
    go_expanded <- go_data %>%
      separate_rows(!!sym(go_col), sep = ",") %>%
      filter(!!sym(go_col) != "-" & !!sym(go_col) != "" & !is.na(!!sym(go_col)))
    
    colnames(go_expanded) <- c("gene", "go_term")
    
    # 构建 TERM2GENE
    term2gene_go <- data.frame(
      term = go_expanded$go_term,
      gene = go_expanded$gene,
      stringsAsFactors = FALSE
    )
    term2gene_go <- unique(term2gene_go)
    
    cat("GO 注释统计：\n")
    cat("  总注释数：", nrow(term2gene_go), "\n")
    cat("  涉及基因数：", length(unique(term2gene_go$gene)), "\n")
    cat("  涉及 GO terms 数：", length(unique(term2gene_go$term)), "\n\n")
    
  } else {
    term2gene_go <- NULL
    cat("⚠️  没有有效的 GO 注释\n\n")
  }
}

# ============================================================
# 第三步：提取 KEGG 通路注释
# ============================================================

cat("========================================\n")
cat("提取 KEGG 通路注释\n")
cat("========================================\n")

# 查找 KEGG 列
kegg_col <- NULL
for(col in c("KEGG_Pathway", "KEGG_pathway", "KEGG", "kegg")) {
  if(col %in% colnames(eggnog)) {
    kegg_col <- col
    break
  }
}

if(is.null(kegg_col)) {
  cat("⚠️  未找到 KEGG 注释列\n\n")
  term2gene_kegg <- NULL
} else {
  cat("找到 KEGG 注释列：", kegg_col, "\n")
  
  # 提取有 KEGG 注释的基因
  kegg_data <- eggnog[eggnog[[kegg_col]] != "-" & 
                      eggnog[[kegg_col]] != "" & 
                      !is.na(eggnog[[kegg_col]]), 
                      c("query", kegg_col)]
  
  cat("有 KEGG 注释的基因数：", nrow(kegg_data), "\n")
  
  if(nrow(kegg_data) > 0) {
    # 拆分多个 KEGG pathways
    kegg_expanded <- kegg_data %>%
      separate_rows(!!sym(kegg_col), sep = ",") %>%
      filter(!!sym(kegg_col) != "-" & !!sym(kegg_col) != "" & !is.na(!!sym(kegg_col)))
    
    colnames(kegg_expanded) <- c("gene", "kegg_pathway")
    
    # 去掉前缀
    kegg_expanded$pathway_id <- gsub("^(ko|map)", "", kegg_expanded$kegg_pathway)
    
    # 构建 TERM2GENE
    term2gene_kegg <- data.frame(
      term = kegg_expanded$pathway_id,
      gene = kegg_expanded$gene,
      stringsAsFactors = FALSE
    )
    term2gene_kegg <- unique(term2gene_kegg)
    
    cat("KEGG 注释统计：\n")
    cat("  总注释数：", nrow(term2gene_kegg), "\n")
    cat("  涉及基因数：", length(unique(term2gene_kegg$gene)), "\n")
    cat("  涉及 KEGG 通路数：", length(unique(term2gene_kegg$term)), "\n\n")
    
  } else {
    term2gene_kegg <- NULL
    cat("⚠️  没有有效的 KEGG 注释\n\n")
  }
}

# ============================================================
# 第四步：准备你的扩张基因家族列表
# ============================================================

cat("========================================\n")
cat("准备扩张基因家族列表\n")
cat("========================================\n")

# 方式1：从文件读取（推荐）
if(file.exists("expanded_genes.txt")) {
  expanded_genes <- read.table("expanded_genes.txt", 
                                header = FALSE, 
                                stringsAsFactors = FALSE)$V1
  cat("从 expanded_genes.txt 读取基因列表\n")
} else {
  # 方式2：直接在脚本中定义（如果没有文件）
  expanded_genes <- c(
    "COD011057-RA",
    "COD011067-RA",
    "COD011052-RA"
    # 👈 在这里添加你的基因
  )
  cat("使用脚本中定义的基因列表\n")
}

cat("扩张基因家族基因数：", length(expanded_genes), "\n")
cat("示例基因：", paste(head(expanded_genes, 5), collapse = ", "), "\n\n")

# ============================================================
# 第五步：检查基因覆盖度
# ============================================================

cat("========================================\n")
cat("检查基因在注释中的覆盖度\n")
cat("========================================\n")

# 检查基因是否在 eggNOG 文件中
all_genes_in_eggnog <- unique(eggnog$query)
genes_found_in_eggnog <- intersect(expanded_genes, all_genes_in_eggnog)
eggnog_coverage <- length(genes_found_in_eggnog) / length(expanded_genes) * 100

cat("eggNOG 文件中的匹配情况：\n")
cat("  输入基因数：", length(expanded_genes), "\n")
cat("  在 eggNOG 中找到的基因数：", length(genes_found_in_eggnog), "\n")
cat("  匹配率：", sprintf("%.1f%%", eggnog_coverage), "\n\n")

if(eggnog_coverage < 80) {
  cat("⚠️  警告：匹配率较低！\n")
  cat("eggNOG 中的基因ID示例：", paste(head(all_genes_in_eggnog, 3), collapse = ", "), "\n")
  cat("你的基因ID示例：", paste(head(expanded_genes, 3), collapse = ", "), "\n\n")
}

# GO 覆盖度
if(!is.null(term2gene_go)) {
  annotated_genes_go <- unique(term2gene_go$gene)
  genes_with_go <- intersect(expanded_genes, annotated_genes_go)
  go_coverage <- length(genes_with_go) / length(expanded_genes) * 100
  
  cat("GO 注释覆盖度：\n")
  cat("  有 GO 注释的基因数：", length(genes_with_go), "\n")
  cat("  覆盖率：", sprintf("%.1f%%", go_coverage), "\n\n")
} else {
  genes_with_go <- c()
  go_coverage <- 0
  cat("⚠️  无 GO 注释数据\n\n")
}

# KEGG 覆盖度
if(!is.null(term2gene_kegg)) {
  annotated_genes_kegg <- unique(term2gene_kegg$gene)
  genes_with_kegg <- intersect(expanded_genes, annotated_genes_kegg)
  kegg_coverage <- length(genes_with_kegg) / length(expanded_genes) * 100
  
  cat("KEGG 注释覆盖度：\n")
  cat("  有 KEGG 注释的基因数：", length(genes_with_kegg), "\n")
  cat("  覆盖率：", sprintf("%.1f%%", kegg_coverage), "\n\n")
} else {
  genes_with_kegg <- c()
  kegg_coverage <- 0
  cat("⚠️  无 KEGG 注释数据\n\n")
}

# ============================================================
# 第六步：GO 富集分析
# ============================================================

cat("========================================\n")
cat("进行 GO 富集分析\n")
cat("========================================\n")

if(!is.null(term2gene_go) && length(genes_with_go) >= 5) {
  
  cat("开始 GO 富集分析...\n")
  
  go_enrich <- enricher(
    gene = expanded_genes,
    TERM2GENE = term2gene_go,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 3,
    maxGSSize = 500,
    universe = unique(term2gene_go$gene)
  )
  
  if(!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
    
    go_result <- as.data.frame(go_enrich)
    cat("✓ 找到", nrow(go_result), "个显著富集的 GO terms\n\n")
    
    # 保存结果
    write.csv(go_result, "GO_enrichment_results.csv", row.names = FALSE)
    cat("结果已保存：GO_enrichment_results.csv\n\n")
    
    # 显示前10个结果
    cat("前10个显著富集的 GO terms：\n")
    print(go_result[1:min(10, nrow(go_result)), c("ID", "pvalue", "p.adjust", "Count", "geneID")])
    cat("\n")
    
    # 生成可视化图表
    cat("生成可视化图表...\n")
    
    pdf("GO_barplot.pdf", width = 10, height = 8)
    print(barplot(go_enrich, showCategory = 20, title = "GO Enrichment Analysis"))
    dev.off()
    
    pdf("GO_dotplot.pdf", width = 10, height = 8)
    print(dotplot(go_enrich, showCategory = 20, title = "GO Enrichment"))
    dev.off()
    
    if(nrow(go_result) >= 5) {
      pdf("GO_cnetplot.pdf", width = 12, height = 10)
      print(cnetplot(go_enrich, categorySize = "pvalue",
                    showCategory = min(5, nrow(go_result)),
                    colorEdge = TRUE, cex_label_gene = 0.7))
      dev.off()
      
      pdf("GO_heatplot.pdf", width = 12, height = 10)
      print(heatplot(go_enrich, showCategory = min(10, nrow(go_result))))
      dev.off()
      
      cat("✓ GO_barplot.pdf\n")
      cat("✓ GO_dotplot.pdf\n")
      cat("✓ GO_cnetplot.pdf\n")
      cat("✓ GO_heatplot.pdf\n\n")
    } else {
      cat("✓ GO_barplot.pdf\n")
      cat("✓ GO_dotplot.pdf\n\n")
    }
    
  } else {
    cat("⚠️  没有找到显著富集的 GO terms\n")
    cat("建议：尝试放宽阈值（pvalueCutoff = 0.1）\n\n")
    go_enrich <- NULL
  }
  
} else {
  cat("⚠️  跳过 GO 富集分析\n")
  if(is.null(term2gene_go)) {
    cat("原因：没有 GO 注释数据\n\n")
  } else {
    cat("原因：有注释的基因数太少（<5）\n\n")
  }
  go_enrich <- NULL
}

# ============================================================
# 第七步：KEGG 通路富集分析
# ============================================================

cat("========================================\n")
cat("进行 KEGG 通路富集分析\n")
cat("========================================\n")

if(!is.null(term2gene_kegg) && length(genes_with_kegg) >= 5) {
  
  cat("开始 KEGG 富集分析...\n")
  
  kegg_enrich <- enricher(
    gene = expanded_genes,
    TERM2GENE = term2gene_kegg,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 3,
    maxGSSize = 500,
    universe = unique(term2gene_kegg$gene)
  )
  
  if(!is.null(kegg_enrich) && nrow(as.data.frame(kegg_enrich)) > 0) {
    
    kegg_result <- as.data.frame(kegg_enrich)
    cat("✓ 找到", nrow(kegg_result), "个显著富集的 KEGG 通路\n\n")
    
    write.csv(kegg_result, "KEGG_enrichment_results.csv", row.names = FALSE)
    cat("结果已保存：KEGG_enrichment_results.csv\n\n")
    
    cat("前10个显著富集的 KEGG 通路：\n")
    print(kegg_result[1:min(10, nrow(kegg_result)), c("ID", "pvalue", "p.adjust", "Count", "geneID")])
    cat("\n")
    
    cat("生成可视化图表...\n")
    
    pdf("KEGG_barplot.pdf", width = 10, height = 8)
    print(barplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment"))
    dev.off()
    
    pdf("KEGG_dotplot.pdf", width = 10, height = 8)
    print(dotplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment"))
    dev.off()
    
    if(nrow(kegg_result) >= 5) {
      pdf("KEGG_cnetplot.pdf", width = 12, height = 10)
      print(cnetplot(kegg_enrich, categorySize = "pvalue",
                    showCategory = min(5, nrow(kegg_result)),
                    colorEdge = TRUE, cex_label_gene = 0.7))
      dev.off()
      
      cat("✓ KEGG_barplot.pdf\n")
      cat("✓ KEGG_dotplot.pdf\n")
      cat("✓ KEGG_cnetplot.pdf\n\n")
    } else {
      cat("✓ KEGG_barplot.pdf\n")
      cat("✓ KEGG_dotplot.pdf\n\n")
    }
    
  } else {
    cat("⚠️  没有找到显著富集的 KEGG 通路\n\n")
    kegg_enrich <- NULL
  }
  
} else {
  cat("⚠️  跳过 KEGG 富集分析\n")
  if(is.null(term2gene_kegg)) {
    cat("原因：没有 KEGG 注释数据\n\n")
  } else {
    cat("原因：有注释的基因数太少（<5）\n\n")
  }
  kegg_enrich <- NULL
}

# ============================================================
# 第八步：生成汇总报告
# ============================================================

cat("========================================\n")
cat("分析汇总报告\n")
cat("========================================\n")

summary_report <- data.frame(
  项目 = c(
    "eggNOG 注释总基因数",
    "扩张基因家族基因数",
    "eggNOG 匹配率",
    "GO 注释覆盖率",
    "KEGG 注释覆盖率",
    "显著富集的 GO terms 数",
    "显著富集的 KEGG 通路数"
  ),
  结果 = c(
    nrow(eggnog),
    length(expanded_genes),
    sprintf("%.1f%%", eggnog_coverage),
    sprintf("%.1f%%", go_coverage),
    sprintf("%.1f%%", kegg_coverage),
    if(!is.null(go_enrich)) nrow(as.data.frame(go_enrich)) else 0,
    if(!is.null(kegg_enrich)) nrow(as.data.frame(kegg_enrich)) else 0
  )
)

print(summary_report)
cat("\n")

write.csv(summary_report, "enrichment_analysis_summary.csv", row.names = FALSE)
cat("汇总报告已保存：enrichment_analysis_summary.csv\n\n")

# 保存注释数据
save(term2gene_go, term2gene_kegg, eggnog, 
     file = "eggnog_annotation_data.RData")
cat("注释数据已保存：eggnog_annotation_data.RData\n\n")

cat("========================================\n")
cat("✅ 分析完成！\n")
cat("========================================\n")

cat("\n生成的文件：\n")
if(!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
  cat("  ✓ GO_enrichment_results.csv\n")
  cat("  ✓ GO_barplot.pdf\n")
  cat("  ✓ GO_dotplot.pdf\n")
  if(nrow(as.data.frame(go_enrich)) >= 5) {
    cat("  ✓ GO_cnetplot.pdf\n")
    cat("  ✓ GO_heatplot.pdf\n")
  }
}
if(!is.null(kegg_enrich) && nrow(as.data.frame(kegg_enrich)) > 0) {
  cat("  ✓ KEGG_enrichment_results.csv\n")
  cat("  ✓ KEGG_barplot.pdf\n")
  cat("  ✓ KEGG_dotplot.pdf\n")
  if(nrow(as.data.frame(kegg_enrich)) >= 5) {
    cat("  ✓ KEGG_cnetplot.pdf\n")
  }
}
cat("  ✓ enrichment_analysis_summary.csv\n")
cat("  ✓ eggnog_annotation_data.RData\n\n")

cat("========================================\n")
cat("脚本执行完毕！\n")
cat("========================================\n")
