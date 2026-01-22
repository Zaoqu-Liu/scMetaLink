#!/usr/bin/env Rscript
# ============================================================================
# scMetaLink Benchmark: 运行 scMetaLink 分析
# ============================================================================
# 
# 目的: 在示例数据集上运行scMetaLink，生成结果用于与MEBOCOST对比
# 数据: crc_example (scMetaLink内置的结直肠癌数据)
# 
# 科学标准:
# 1. 使用相同的输入数据
# 2. 使用合理的默认参数
# 3. 记录所有参数设置
# 4. 保存完整的中间结果
# ============================================================================

suppressPackageStartupMessages({
  library(scMetaLink)
  library(Matrix)
  library(dplyr)
})

# 设置工作目录和输出目录
# 获取脚本目录 (支持命令行和RStudio)
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

output_dir <- "results/scMetaLink"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=" , paste(rep("=", 60), collapse = ""), "\n")
cat("scMetaLink Benchmark Analysis\n")
cat("=" , paste(rep("=", 60), collapse = ""), "\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("输出目录:", normalizePath(output_dir), "\n\n")

# ============================================================================
# 1. 加载数据
# ============================================================================
cat("\n=== 1. 加载数据 ===\n")

data(crc_example, package = "scMetaLink")

cat("表达矩阵维度:", dim(crc_expr)[1], "genes x", dim(crc_expr)[2], "cells\n")
cat("细胞类型分布:\n")
print(table(crc_meta$cell_type))

# 保存数据信息
data_info <- list(
  n_genes = nrow(crc_expr),
  n_cells = ncol(crc_expr),
  cell_types = unique(crc_meta$cell_type),
  n_cell_types = length(unique(crc_meta$cell_type)),
  sparsity = 1 - sum(crc_expr > 0) / prod(dim(crc_expr))
)
saveRDS(data_info, file.path(output_dir, "data_info.rds"))

# ============================================================================
# 2. 创建 scMetaLink 对象
# ============================================================================
cat("\n=== 2. 创建 scMetaLink 对象 ===\n")

start_time <- Sys.time()

obj <- createScMetaLink(
  expression_data = crc_expr,
  cell_meta = crc_meta,
  cell_type_column = "cell_type",
  min_cells = 10  # 至少10个细胞的细胞类型才保留
)

creation_time <- difftime(Sys.time(), start_time, units = "secs")
cat("对象创建时间:", round(creation_time, 2), "秒\n")

# ============================================================================
# 3. 推断代谢物生产 (Metabolite Production)
# ============================================================================
cat("\n=== 3. 推断代谢物生产 ===\n")

# 参数说明:
# - method: "combined" = mean_expr * pct_expr (综合表达量和表达比例)
# - consider_degradation: TRUE = 考虑降解酶的影响
# - consider_secretion: TRUE = 考虑分泌潜力
# - normalize: TRUE = Z-score标准化后缩放到[0,1]

production_params <- list(
  method = "combined",
  consider_degradation = TRUE,
  consider_secretion = TRUE,
  normalize = TRUE,
  min_expression = 0
)

start_time <- Sys.time()

obj <- inferProduction(
  obj,
  method = production_params$method,
  consider_degradation = production_params$consider_degradation,
  consider_secretion = production_params$consider_secretion,
  normalize = production_params$normalize,
  min_expression = production_params$min_expression,
  verbose = TRUE
)

production_time <- difftime(Sys.time(), start_time, units = "secs")
cat("生产推断时间:", round(production_time, 2), "秒\n")

# 保存生产分数
saveRDS(obj@production_scores, file.path(output_dir, "production_scores.rds"))

# ============================================================================
# 4. 推断代谢物感知 (Metabolite Sensing)
# ============================================================================
cat("\n=== 4. 推断代谢物感知 ===\n")

# 参数说明:
# - weight_by_affinity: TRUE = 根据受体亲和力加权
# - include_transporters: TRUE = 包含转运体
# - use_hill: FALSE = 不使用Hill函数转换 (保持线性)
# - normalize: TRUE = Z-score标准化后缩放到[0,1]

sensing_params <- list(
  method = "combined",
  weight_by_affinity = TRUE,
  include_transporters = TRUE,
  use_hill = FALSE,
  normalize = TRUE
)

start_time <- Sys.time()

obj <- inferSensing(
  obj,
  method = sensing_params$method,
  weight_by_affinity = sensing_params$weight_by_affinity,
  include_transporters = sensing_params$include_transporters,
  use_hill = sensing_params$use_hill,
  normalize = sensing_params$normalize,
  verbose = TRUE
)

sensing_time <- difftime(Sys.time(), start_time, units = "secs")
cat("感知推断时间:", round(sensing_time, 2), "秒\n")

# 保存感知分数
saveRDS(obj@sensing_scores, file.path(output_dir, "sensing_scores.rds"))

# ============================================================================
# 5. 计算通讯分数 (Communication Scores)
# ============================================================================
cat("\n=== 5. 计算通讯分数 ===\n")

# 参数说明:
# - method: "geometric" = sqrt(production * sensing)
# - n_permutations: 1000 = 置换检验次数
# - min_production/min_sensing: 0.1 = 最低阈值

communication_params <- list(
  method = "geometric",
  n_permutations = 1000,
  min_production = 0.1,
  min_sensing = 0.1,
  n_cores = 1,
  seed = 42
)

start_time <- Sys.time()

obj <- computeCommunication(
  obj,
  method = communication_params$method,
  n_permutations = communication_params$n_permutations,
  min_production = communication_params$min_production,
  min_sensing = communication_params$min_sensing,
  n_cores = communication_params$n_cores,
  seed = communication_params$seed,
  verbose = TRUE
)

communication_time <- difftime(Sys.time(), start_time, units = "secs")
cat("通讯计算时间:", round(communication_time, 2), "秒\n")

# 保存通讯分数
saveRDS(obj@communication_scores, file.path(output_dir, "communication_scores.rds"))
saveRDS(obj@communication_pvalues, file.path(output_dir, "communication_pvalues.rds"))

# ============================================================================
# 6. 筛选显著通讯
# ============================================================================
cat("\n=== 6. 筛选显著通讯 ===\n")

# 使用代谢物分层FDR校正 (scMetaLink独特功能)
filter_params <- list(
  pvalue_threshold = 0.05,
  adjust_method = "metabolite_stratified",  # scMetaLink独特的校正方法
  min_score = 0
)

obj <- filterSignificantInteractions(
  obj,
  pvalue_threshold = filter_params$pvalue_threshold,
  adjust_method = filter_params$adjust_method,
  min_score = filter_params$min_score
)

sig_interactions <- obj@significant_interactions
cat("显著通讯数量:", nrow(sig_interactions), "\n")

# 保存显著通讯
saveRDS(sig_interactions, file.path(output_dir, "significant_interactions.rds"))
write.csv(sig_interactions, file.path(output_dir, "significant_interactions.csv"), row.names = FALSE)

# ============================================================================
# 7. 生成摘要统计
# ============================================================================
cat("\n=== 7. 生成摘要统计 ===\n")

# 总体统计
summary_stats <- list(
  # 数据信息
  n_genes = data_info$n_genes,
  n_cells = data_info$n_cells,
  n_cell_types = data_info$n_cell_types,
  
  # 结果统计
  n_metabolites_production = nrow(obj@production_scores),
  n_metabolites_sensing = nrow(obj@sensing_scores),
  n_significant_interactions = nrow(sig_interactions),
  
  # 通讯统计
  n_unique_metabolites = length(unique(sig_interactions$metabolite_name)),
  n_unique_senders = length(unique(sig_interactions$sender)),
  n_unique_receivers = length(unique(sig_interactions$receiver)),
  
  # 计算时间
  time_creation = as.numeric(creation_time),
  time_production = as.numeric(production_time),
  time_sensing = as.numeric(sensing_time),
  time_communication = as.numeric(communication_time),
  time_total = as.numeric(creation_time + production_time + sensing_time + communication_time),
  
  # 参数
  params_production = production_params,
  params_sensing = sensing_params,
  params_communication = communication_params,
  params_filter = filter_params
)

saveRDS(summary_stats, file.path(output_dir, "summary_stats.rds"))

# 打印摘要
cat("\n--- 结果摘要 ---\n")
cat("代谢物数量 (生产):", summary_stats$n_metabolites_production, "\n")
cat("代谢物数量 (感知):", summary_stats$n_metabolites_sensing, "\n")
cat("显著通讯数量:", summary_stats$n_significant_interactions, "\n")
cat("涉及代谢物数:", summary_stats$n_unique_metabolites, "\n")
cat("涉及发送细胞类型:", summary_stats$n_unique_senders, "\n")
cat("涉及接收细胞类型:", summary_stats$n_unique_receivers, "\n")
cat("总运行时间:", round(summary_stats$time_total, 2), "秒\n")

# ============================================================================
# 8. 提取关键代谢物的通讯 (用于Ground Truth验证)
# ============================================================================
cat("\n=== 8. 提取关键代谢物通讯 ===\n")

# 关注的代谢物 (Ground Truth相关)
key_metabolites <- c(
  "L-Lactic acid",      # 乳酸 - 肿瘤代谢
  "Adenosine",          # 腺苷 - 免疫抑制
  "L-Glutamine",        # 谷氨酰胺
  "L-Glutamic acid",    # 谷氨酸
  "Glucose",            # 葡萄糖
  "ATP"                 # ATP - DAMP信号
)

key_results <- sig_interactions %>%
  filter(metabolite_name %in% key_metabolites) %>%
  arrange(metabolite_name, pvalue_adjusted)

cat("关键代谢物通讯数量:", nrow(key_results), "\n")
if (nrow(key_results) > 0) {
  cat("\n关键代谢物通讯预览:\n")
  print(head(key_results[, c("sender", "receiver", "metabolite_name", 
                              "communication_score", "pvalue_adjusted")], 20))
}

saveRDS(key_results, file.path(output_dir, "key_metabolite_interactions.rds"))
write.csv(key_results, file.path(output_dir, "key_metabolite_interactions.csv"), row.names = FALSE)

# ============================================================================
# 9. 保存完整对象
# ============================================================================
cat("\n=== 9. 保存完整对象 ===\n")

saveRDS(obj, file.path(output_dir, "scMetaLink_object.rds"))
cat("完整对象已保存:", file.path(output_dir, "scMetaLink_object.rds"), "\n")

# ============================================================================
# 完成
# ============================================================================
cat("\n")
cat("=" , paste(rep("=", 60), collapse = ""), "\n")
cat("scMetaLink 分析完成!\n")
cat("=" , paste(rep("=", 60), collapse = ""), "\n")
cat("结果保存在:", normalizePath(output_dir), "\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
