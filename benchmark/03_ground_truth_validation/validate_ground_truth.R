#!/usr/bin/env Rscript
# ============================================================================
# Ground Truth 验证脚本
# ============================================================================
# 
# 目的: 验证scMetaLink是否能检测到已知的代谢物-细胞通讯关系
# 
# 科学原理:
# 1. 使用文献中已验证的代谢物通讯作为"Ground Truth"
# 2. 检查scMetaLink是否能检测到这些已知通讯
# 3. 计算检测率、排名等指标
# 
# Ground Truth来源:
# - 乳酸(Lactate) → 巨噬细胞M2极化 (PMID: 30431439, 29468929)
# - 乳酸 → CD8+ T细胞抑制 (PMID: 33603751)
# - 腺苷(Adenosine) → T细胞抑制 (PMID: 27066002)
# - 等
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# 设置工作目录
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(sub("--file=", "", file_arg)))
  }
  return(".")
}
script_dir <- get_script_dir()
setwd(script_dir)

output_dir <- "results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=" , paste(rep("=", 60), collapse = ""), "\n")
cat("Ground Truth 验证\n")
cat("=" , paste(rep("=", 60), collapse = ""), "\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================================
# 1. 加载Ground Truth数据
# ============================================================================
cat("\n=== 1. 加载Ground Truth ===\n")

ground_truth <- read.csv("known_interactions.csv", stringsAsFactors = FALSE)
cat("Ground Truth 通讯数量:", nrow(ground_truth), "\n")

# 显示Ground Truth
cat("\nGround Truth 代谢物-细胞通讯:\n")
print(ground_truth[, c("metabolite", "sender_cell_type", "receiver_cell_type", "expected_effect")])

# ============================================================================
# 2. 加载scMetaLink结果
# ============================================================================
cat("\n=== 2. 加载scMetaLink结果 ===\n")

scMetaLink_dir <- "../02_algorithm_comparison/results/scMetaLink"

if (!file.exists(file.path(scMetaLink_dir, "significant_interactions.rds"))) {
  stop("请先运行 02_algorithm_comparison/01_run_scMetaLink.R")
}

sig_interactions <- readRDS(file.path(scMetaLink_dir, "significant_interactions.rds"))
all_comm_scores <- readRDS(file.path(scMetaLink_dir, "communication_scores.rds"))
all_pvalues <- readRDS(file.path(scMetaLink_dir, "communication_pvalues.rds"))

cat("显著通讯总数:", nrow(sig_interactions), "\n")
cat("可用细胞类型:", paste(dimnames(all_comm_scores)[[1]], collapse = ", "), "\n")

# ============================================================================
# 3. 细胞类型映射
# ============================================================================
cat("\n=== 3. 细胞类型映射 ===\n")

# 数据集中的细胞类型
available_types <- dimnames(all_comm_scores)[[1]]
cat("数据集细胞类型:", paste(available_types, collapse = ", "), "\n")

# 创建映射表 (Ground Truth细胞类型 -> 数据集细胞类型)
# 注意: 这需要根据实际数据集调整
cell_type_mapping <- list(
  "Tumor" = c("Tumor", "Cancer", "Malignant", "Epithelial"),
  "Macrophage" = c("Macrophage", "Mono/Macro", "TAM", "Myeloid"),
  "CD8 T cell" = c("CD8", "CD8 T", "CD8+ T", "T cells CD8"),
  "T cell" = c("T cell", "T cells", "CD4", "CD8", "Treg"),
  "Treg" = c("Treg", "T reg", "Regulatory T"),
  "Dendritic cell" = c("DC", "Dendritic", "cDC", "pDC"),
  "NK cell" = c("NK", "NK cell"),
  "Endothelial" = c("Endothelial", "Endo"),
  "Fibroblast" = c("Fibroblast", "CAF", "Stromal"),
  "B cell" = c("B cell", "B cells", "Plasma")
)

# 函数: 找到匹配的细胞类型
find_matching_celltype <- function(gt_type, available_types, mapping) {
  # 直接匹配
  if (gt_type %in% available_types) return(gt_type)
  
  # 通过映射查找
  if (gt_type %in% names(mapping)) {
    for (alias in mapping[[gt_type]]) {
      matches <- grep(alias, available_types, ignore.case = TRUE, value = TRUE)
      if (length(matches) > 0) return(matches[1])
    }
  }
  
  # 部分匹配
  matches <- grep(gt_type, available_types, ignore.case = TRUE, value = TRUE)
  if (length(matches) > 0) return(matches[1])
  
  return(NA)
}

# ============================================================================
# 4. 验证每个Ground Truth通讯
# ============================================================================
cat("\n=== 4. 验证Ground Truth通讯 ===\n")

validation_results <- data.frame()

for (i in seq_len(nrow(ground_truth))) {
  gt <- ground_truth[i, ]
  
  # 映射细胞类型
  sender_mapped <- find_matching_celltype(gt$sender_cell_type, available_types, cell_type_mapping)
  receiver_mapped <- find_matching_celltype(gt$receiver_cell_type, available_types, cell_type_mapping)
  
  # 检查代谢物是否存在
  metabolite_found <- gt$metabolite %in% dimnames(all_comm_scores)[[3]]
  
  # 初始化结果
  result <- data.frame(
    metabolite = gt$metabolite,
    gt_sender = gt$sender_cell_type,
    gt_receiver = gt$receiver_cell_type,
    mapped_sender = sender_mapped,
    mapped_receiver = receiver_mapped,
    expected_effect = gt$expected_effect,
    evidence_strength = gt$evidence_strength,
    metabolite_found = metabolite_found,
    celltype_mapped = !is.na(sender_mapped) & !is.na(receiver_mapped),
    communication_score = NA,
    pvalue = NA,
    pvalue_adjusted = NA,
    is_significant = FALSE,
    rank_in_metabolite = NA,
    rank_overall = NA,
    validation_status = "Not testable",
    stringsAsFactors = FALSE
  )
  
  # 如果可以测试
  if (metabolite_found && !is.na(sender_mapped) && !is.na(receiver_mapped)) {
    # 提取通讯分数
    score <- all_comm_scores[sender_mapped, receiver_mapped, gt$metabolite]
    pval <- all_pvalues[sender_mapped, receiver_mapped, gt$metabolite]
    
    result$communication_score <- score
    result$pvalue <- pval
    
    # 检查是否在显著结果中
    sig_match <- sig_interactions %>%
      filter(metabolite_name == gt$metabolite,
             sender == sender_mapped,
             receiver == receiver_mapped)
    
    if (nrow(sig_match) > 0) {
      result$is_significant <- TRUE
      result$pvalue_adjusted <- sig_match$pvalue_adjusted[1]
      result$validation_status <- "Detected (Significant)"
    } else if (score > 0) {
      result$validation_status <- "Detected (Not Significant)"
    } else {
      result$validation_status <- "Not Detected"
    }
    
    # 计算排名
    # 在该代谢物内的排名
    met_scores <- all_comm_scores[, , gt$metabolite]
    met_scores_vec <- as.vector(met_scores)
    result$rank_in_metabolite <- sum(met_scores_vec >= score, na.rm = TRUE)
    
    # 总体排名
    all_scores_vec <- as.vector(all_comm_scores)
    result$rank_overall <- sum(all_scores_vec >= score, na.rm = TRUE)
    
  } else if (!metabolite_found) {
    result$validation_status <- "Metabolite not in database"
  } else {
    result$validation_status <- "Cell type not in data"
  }
  
  validation_results <- rbind(validation_results, result)
}

# ============================================================================
# 5. 计算验证统计
# ============================================================================
cat("\n=== 5. 验证统计 ===\n")

# 可测试的数量
testable <- validation_results %>% filter(metabolite_found & celltype_mapped)
n_testable <- nrow(testable)

# 检测到的数量
detected <- testable %>% filter(validation_status %in% c("Detected (Significant)", "Detected (Not Significant)"))
n_detected <- nrow(detected)

# 显著的数量
significant <- testable %>% filter(is_significant)
n_significant <- nrow(significant)

# 计算检测率
detection_rate <- ifelse(n_testable > 0, n_detected / n_testable, 0)
significance_rate <- ifelse(n_testable > 0, n_significant / n_testable, 0)

cat("\nGround Truth 验证结果:\n")
cat("-----------------------------------\n")
cat("总Ground Truth数量:", nrow(ground_truth), "\n")
cat("可测试数量:", n_testable, "\n")
cat("检测到数量:", n_detected, "\n")
cat("显著数量:", n_significant, "\n")
cat("检测率:", sprintf("%.1f%%", detection_rate * 100), "\n")
cat("显著率:", sprintf("%.1f%%", significance_rate * 100), "\n")

# 按证据强度分组
cat("\n按证据强度分组:\n")
by_strength <- validation_results %>%
  filter(metabolite_found & celltype_mapped) %>%
  group_by(evidence_strength) %>%
  summarise(
    n_total = n(),
    n_detected = sum(validation_status %in% c("Detected (Significant)", "Detected (Not Significant)")),
    n_significant = sum(is_significant),
    detection_rate = n_detected / n_total,
    significance_rate = n_significant / n_total,
    .groups = "drop"
  )
print(by_strength)

# ============================================================================
# 6. 详细结果输出
# ============================================================================
cat("\n=== 6. 详细结果 ===\n")

# 显示验证结果
cat("\n验证详情:\n")
print(validation_results[, c("metabolite", "gt_sender", "gt_receiver", 
                             "validation_status", "communication_score", 
                             "pvalue", "is_significant")])

# 保存结果
write.csv(validation_results, file.path(output_dir, "validation_results.csv"), row.names = FALSE)
saveRDS(validation_results, file.path(output_dir, "validation_results.rds"))

# ============================================================================
# 7. 生成验证报告图
# ============================================================================
cat("\n=== 7. 生成验证报告图 ===\n")

# 图1: 验证状态分布
p1 <- ggplot(validation_results, aes(x = validation_status)) +
  geom_bar(aes(fill = validation_status), color = "black") +
  scale_fill_manual(values = c(
    "Detected (Significant)" = "#2E7D32",
    "Detected (Not Significant)" = "#FFA726",
    "Not Detected" = "#D32F2F",
    "Metabolite not in database" = "#9E9E9E",
    "Cell type not in data" = "#BDBDBD",
    "Not testable" = "#E0E0E0"
  )) +
  labs(
    title = "Ground Truth 验证结果",
    subtitle = sprintf("检测率: %.1f%%, 显著率: %.1f%%", detection_rate * 100, significance_rate * 100),
    x = "验证状态",
    y = "数量"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(file.path(output_dir, "fig1_validation_status.pdf"), p1, width = 8, height = 6)
ggsave(file.path(output_dir, "fig1_validation_status.png"), p1, width = 8, height = 6, dpi = 150)

# 图2: 通讯分数分布 (仅可测试的)
testable_results <- validation_results %>% 
  filter(metabolite_found & celltype_mapped & !is.na(communication_score))

if (nrow(testable_results) > 0) {
  p2 <- ggplot(testable_results, aes(x = reorder(paste(metabolite, gt_sender, "->", gt_receiver), -communication_score),
                                      y = communication_score)) +
    geom_bar(stat = "identity", aes(fill = is_significant), color = "black") +
    scale_fill_manual(values = c("TRUE" = "#2E7D32", "FALSE" = "#FFA726"),
                      labels = c("TRUE" = "显著", "FALSE" = "不显著")) +
    labs(
      title = "Ground Truth 通讯分数",
      x = "代谢物-细胞对",
      y = "通讯分数",
      fill = "统计显著性"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  
  ggsave(file.path(output_dir, "fig2_communication_scores.pdf"), p2, width = 12, height = 6)
  ggsave(file.path(output_dir, "fig2_communication_scores.png"), p2, width = 12, height = 6, dpi = 150)
}

# 图3: 按代谢物分组的验证结果
met_summary <- validation_results %>%
  filter(metabolite_found & celltype_mapped) %>%
  group_by(metabolite) %>%
  summarise(
    n_total = n(),
    n_significant = sum(is_significant),
    mean_score = mean(communication_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_significant))

if (nrow(met_summary) > 0) {
  p3 <- ggplot(met_summary, aes(x = reorder(metabolite, n_significant), y = n_significant)) +
    geom_bar(stat = "identity", fill = "#1976D2", color = "black") +
    geom_text(aes(label = paste0(n_significant, "/", n_total)), hjust = -0.2) +
    coord_flip() +
    labs(
      title = "按代谢物分组的验证结果",
      x = "代谢物",
      y = "显著通讯数 / 总数"
    ) +
    theme_minimal()
  
  ggsave(file.path(output_dir, "fig3_by_metabolite.pdf"), p3, width = 10, height = 6)
  ggsave(file.path(output_dir, "fig3_by_metabolite.png"), p3, width = 10, height = 6, dpi = 150)
}

# ============================================================================
# 8. 生成验证摘要
# ============================================================================
cat("\n=== 8. 生成验证摘要 ===\n")

validation_summary <- list(
  total_ground_truth = nrow(ground_truth),
  n_testable = n_testable,
  n_detected = n_detected,
  n_significant = n_significant,
  detection_rate = detection_rate,
  significance_rate = significance_rate,
  by_strength = by_strength,
  by_metabolite = met_summary,
  timestamp = Sys.time()
)

saveRDS(validation_summary, file.path(output_dir, "validation_summary.rds"))

# ============================================================================
# 完成
# ============================================================================
cat("\n")
cat("=" , paste(rep("=", 60), collapse = ""), "\n")
cat("Ground Truth 验证完成!\n")
cat("=" , paste(rep("=", 60), collapse = ""), "\n")
cat("结果保存在:", normalizePath(output_dir), "\n")
cat("\n主要发现:\n")
cat("- 检测率:", sprintf("%.1f%%", detection_rate * 100), "\n")
cat("- 显著率:", sprintf("%.1f%%", significance_rate * 100), "\n")
cat("\n")
