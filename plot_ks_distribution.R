# Ks 分布曲线图绘制脚本
# 用于绘制 Cdom、Crad、Cdom_Crad 三个物种/比较类型的 Ks 分布

# 加载必要的包
library(ggplot2)

# ============== 参数设置 ==============
input_file <- "all.ks"  # 输入文件路径，请根据实际情况修改
output_file <- "Ks_distribution.pdf"  # 输出文件名

# Ks 值范围和区间设置
ks_min <- 0
ks_max <- 5
bin_width <- 0.1

# ============== 读取数据 ==============
data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 检查列名并标准化
colnames(data) <- c("pair", "Ks", "Species")

# 过滤 Ks > 5 的数据
data <- data[data$Ks <= ks_max, ]

# 查看数据概览
cat("数据概览：\n")
cat("总基因对数量：", nrow(data), "\n")
cat("各类型基因对数量：\n")
print(table(data$Species))

# ============== 计算 Ks 分布 ==============
# 创建 Ks 区间
breaks <- seq(ks_min, ks_max, by = bin_width)

# 为每个物种类型计算频率分布
calc_distribution <- function(ks_values, breaks) {
  # 使用 hist 计算频率，但不绘图
  h <- hist(ks_values, breaks = breaks, plot = FALSE)
  # 返回区间中点和计数
  data.frame(
    Ks = h$mids,
    Count = h$counts
  )
}

# 分别计算三个类型的分布
species_list <- unique(data$Species)
dist_list <- list()

for (sp in species_list) {
  sp_data <- data[data$Species == sp, ]
  dist <- calc_distribution(sp_data$Ks, breaks)
  dist$Species <- sp
  dist_list[[sp]] <- dist
}

# 合并所有分布数据
plot_data <- do.call(rbind, dist_list)

# ============== 绘制曲线图 ==============
# 设置颜色
colors <- c("Cdis" = "#E41A1C", "Cdis_Cmic" = "#377EB8", "Cdis_Cdom" = "#4DAF4A", "Cdis_Crad" = "#FF7F00")

# 绘图
p <- ggplot(plot_data, aes(x = Ks, y = Count, color = Species)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.0, alpha = 0.5) +
  scale_color_manual(values = colors,
                     labels = c("Cdis" = "Cdis",
                                "Cdis_Cmic" = "Cdis_Cmic",
                                "Cdis_Cdom" = "Cdis_Cdom",
                                "Cdis_Crad" = "Cdis_Crad")) +
  scale_x_continuous(breaks = seq(0, 5, by = 0.5),
                     limits = c(0, 5)) +
  labs(x = "Ks",
       y = "Duplications") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = c(0.75, 0.75),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

# 保存图片
ggsave(output_file, p, width = 10, height = 6, dpi = 600)
cat("\n图片已保存至：", output_file, "\n")

# 同时保存 PNG 格式
ggsave(gsub(".pdf", ".png", output_file), p, width = 10, height = 6, dpi = 600)
cat("PNG 格式已保存至：", gsub(".pdf", ".png", output_file), "\n")

# ============== 输出统计信息 ==============
cat("\n============== Ks 统计信息 ==============\n")
for (sp in species_list) {
  sp_data <- data[data$Species == sp, ]
  cat("\n", sp, ":\n", sep = "")
  cat("  基因对数量：", nrow(sp_data), "\n")
  cat("  Ks 均值：", round(mean(sp_data$Ks), 4), "\n")
  cat("  Ks 中位数：", round(median(sp_data$Ks), 4), "\n")
  cat("  Ks 范围：", round(min(sp_data$Ks), 4), " - ", round(max(sp_data$Ks), 4), "\n")
}
