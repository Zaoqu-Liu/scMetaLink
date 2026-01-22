#!/usr/bin/env Rscript
# ============================================================================
# scMetaLink 空间分析优势展示
# ============================================================================
# 
# 目的: 展示scMetaLink独特的空间转录组分析能力
# (MEBOCOST不具备原生空间分析功能)
# 
# 科学意义:
# 1. 代谢物是扩散性信号分子，空间距离直接影响通讯强度
# 2. 空间分析可以验证预测的通讯是否发生在空间邻近的细胞之间
# 3. 可以识别通讯热点区域
# 
# 方法:
# 1. 距离加权通讯: C_spatial = sqrt(P × S) × w(d)
# 2. 空间验证指标: Distance Enrichment Score (DES)
# ============================================================================

suppressPackageStartupMessages({
  library(scMetaLink)
  library(SpatialExperiment)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
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
cat("scMetaLink 空间分析优势展示\n")
cat("=" , paste(rep("=", 60), collapse = ""), "\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================================
# 1. 检查空间数据
# ============================================================================
cat("\n=== 1. 检查空间数据 ===\n")

# 尝试加载TENxVisiumData
spatial_data_path <- "../data/spatial/breast_cancer_visium.rds"
use_builtin_data <- FALSE

if (file.exists(spatial_data_path)) {
  cat("加载已下载的空间数据...\n")
  spe <- readRDS(spatial_data_path)
  cat("数据来源: TENxVisiumData (HumanBreastCancerIDC)\n")
} else {
  # 尝试直接从TENxVisiumData加载
  cat("尝试从TENxVisiumData加载数据...\n")
  tryCatch({
    library(TENxVisiumData)
    spe <- HumanBreastCancerIDC()
    saveRDS(spe, spatial_data_path)
    cat("数据已下载并保存\n")
  }, error = function(e) {
    cat("空间数据不可用，将使用scMetaLink内置的st_colon数据\n")
    use_builtin_data <<- TRUE
  })
}

# 如果使用内置数据
if (use_builtin_data) {
  tryCatch({
    data(st_colon, package = "scMetaLink")
    cat("使用scMetaLink内置空间数据 (st_colon)\n")
    
    # 创建SpatialExperiment对象或直接使用提供的数据
    expr_data <- st_expr
    spatial_coords <- st_meta[, c("x", "y")]
    cell_meta <- st_meta
    scale_factors <- st_scalefactors
    
    # 创建scMetaLink空间对象
    obj_spatial <- createScMetaLinkFromSpatial(
      expression_data = expr_data,
      spatial_coords = spatial_coords,
      cell_meta = cell_meta,
      cell_type_column = "cell_type",
      scale_factors = scale_factors,
      min_cells = 5
    )
    
  }, error = function(e) {
    cat("内置空间数据也不可用:", e$message, "\n")
    cat("\n请先下载空间数据:\n")
    cat("  setwd('benchmark/data')\n")
    cat("  source('download_spatial_data.R')\n")
    stop("空间数据不可用")
  })
} else {
  # 从SpatialExperiment创建scMetaLink对象
  cat("\n处理SpatialExperiment数据...\n")
  
  # 提取数据
  expr_data <- assay(spe, "counts")
  spatial_coords <- spatialCoords(spe)
  
  # 创建简单的细胞类型注释 (基于marker基因)
  # 注意: 实际应用中应使用更好的注释方法
  cat("创建细胞类型注释...\n")
  
  # 归一化
  col_sums <- Matrix::colSums(expr_data)
  expr_norm <- sweep(expr_data, 2, col_sums, "/") * 10000
  expr_norm <- log1p(expr_norm)
  
  # 基于marker基因的简单分类
  # 这里使用简化的方法，实际应用建议使用专业的细胞类型注释工具
  markers <- list(
    Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19"),
    Tumor = c("MKI67", "TOP2A", "PCNA"),
    T_cell = c("CD3D", "CD3E", "CD8A", "CD4"),
    Macrophage = c("CD68", "CD163", "CSF1R", "MARCO"),
    B_cell = c("CD79A", "CD79B", "MS4A1"),
    Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM"),
    Endothelial = c("PECAM1", "VWF", "CDH5")
  )
  
  # 计算每个spot的细胞类型得分
  cell_types <- character(ncol(expr_data))
  
  for (i in seq_len(ncol(expr_data))) {
    spot_expr <- expr_norm[, i]
    scores <- sapply(markers, function(m) {
      available_markers <- m[m %in% rownames(expr_data)]
      if (length(available_markers) == 0) return(0)
      mean(spot_expr[available_markers], na.rm = TRUE)
    })
    cell_types[i] <- names(which.max(scores))
  }
  
  # 创建metadata
  cell_meta <- data.frame(
    barcode = colnames(expr_data),
    cell_type = cell_types,
    row.names = colnames(expr_data)
  )
  
  cat("细胞类型分布:\n")
  print(table(cell_types))
  
  # 创建scMetaLink空间对象
  obj_spatial <- createScMetaLinkFromSpatial(
    expression_data = expr_data,
    spatial_coords = as.data.frame(spatial_coords),
    cell_meta = cell_meta,
    cell_type_column = "cell_type",
    scale_factors = list(pixels_per_um = 1),  # 假设坐标已经是微米
    min_cells = 5
  )
}

cat("\n空间对象创建完成\n")
print(obj_spatial)

# ============================================================================
# 2. 推断代谢物生产和感知
# ============================================================================
cat("\n=== 2. 推断代谢物生产和感知 ===\n")

obj_spatial <- inferProduction(
  obj_spatial,
  method = "combined",
  consider_degradation = TRUE,
  consider_secretion = TRUE,
  verbose = TRUE
)

obj_spatial <- inferSensing(
  obj_spatial,
  method = "combined",
  weight_by_affinity = TRUE,
  include_transporters = TRUE,
  verbose = TRUE
)

# ============================================================================
# 3. 空间通讯分析 (scMetaLink独特功能)
# ============================================================================
cat("\n=== 3. 空间通讯分析 (scMetaLink独特功能) ===\n")

# 3.1 KNN方法 (推荐用于Visium数据)
cat("\n3.1 KNN空间加权方法...\n")
obj_knn <- computeSpatialCommunication(
  obj_spatial,
  method = "knn",
  k_neighbors = 6,          # Visium六边形网格
  symmetric = TRUE,
  comm_method = "geometric",
  min_production = 0.1,
  min_sensing = 0.1,
  n_permutations = 100,     # 演示用，正式分析用1000+
  n_cores = 1,
  seed = 42,
  verbose = TRUE
)

# 3.2 Gaussian方法
cat("\n3.2 Gaussian空间加权方法...\n")
obj_gaussian <- computeSpatialCommunication(
  obj_spatial,
  method = "gaussian",
  sigma = 50,               # 50微米特征衰减距离
  distance_threshold = 150, # 最大150微米
  comm_method = "geometric",
  min_production = 0.1,
  min_sensing = 0.1,
  n_permutations = 100,
  n_cores = 1,
  seed = 42,
  verbose = TRUE
)

# ============================================================================
# 4. 筛选显著空间通讯
# ============================================================================
cat("\n=== 4. 筛选显著空间通讯 ===\n")

obj_knn <- filterSignificantInteractions(obj_knn, pvalue_threshold = 0.05)
obj_gaussian <- filterSignificantInteractions(obj_gaussian, pvalue_threshold = 0.05)

cat("KNN方法显著通讯:", nrow(obj_knn@significant_interactions), "\n")
cat("Gaussian方法显著通讯:", nrow(obj_gaussian@significant_interactions), "\n")

# 保存结果
saveRDS(obj_knn, file.path(output_dir, "spatial_obj_knn.rds"))
saveRDS(obj_gaussian, file.path(output_dir, "spatial_obj_gaussian.rds"))

# ============================================================================
# 5. 空间可视化
# ============================================================================
cat("\n=== 5. 空间可视化 ===\n")

# 5.1 空间细胞类型分布
cat("生成空间细胞类型分布图...\n")
tryCatch({
  p1 <- plotSpatialCellTypes(obj_knn, point_size = 1.5, alpha = 0.8)
  ggsave(file.path(output_dir, "fig1_spatial_celltypes.pdf"), p1, width = 10, height = 8)
  ggsave(file.path(output_dir, "fig1_spatial_celltypes.png"), p1, width = 10, height = 8, dpi = 150)
}, error = function(e) cat("  跳过:", e$message, "\n"))

# 5.2 乳酸生产的空间分布
cat("生成乳酸生产空间分布图...\n")
tryCatch({
  p2 <- plotSpatialFeature(obj_knn, metabolite = "L-Lactic acid", type = "production", point_size = 1.5)
  ggsave(file.path(output_dir, "fig2_lactate_production_spatial.pdf"), p2, width = 10, height = 8)
  ggsave(file.path(output_dir, "fig2_lactate_production_spatial.png"), p2, width = 10, height = 8, dpi = 150)
}, error = function(e) cat("  跳过:", e$message, "\n"))

# 5.3 空间通讯网络
cat("生成空间通讯网络图...\n")
tryCatch({
  p3 <- plotSpatialCommunicationNetwork(obj_knn, metabolite = NULL, top_n = 15, arrow_scale = 1.5)
  ggsave(file.path(output_dir, "fig3_spatial_network.pdf"), p3, width = 12, height = 10)
  ggsave(file.path(output_dir, "fig3_spatial_network.png"), p3, width = 12, height = 10, dpi = 150)
}, error = function(e) cat("  跳过:", e$message, "\n"))

# ============================================================================
# 6. 空间方法对比
# ============================================================================
cat("\n=== 6. 空间方法对比 ===\n")

# 对比KNN和Gaussian方法的结果
sig_knn <- obj_knn@significant_interactions
sig_gaussian <- obj_gaussian@significant_interactions

# 创建唯一标识符
if (nrow(sig_knn) > 0 && nrow(sig_gaussian) > 0) {
  sig_knn$interaction_id <- paste(sig_knn$sender, sig_knn$receiver, sig_knn$metabolite_name, sep = "_")
  sig_gaussian$interaction_id <- paste(sig_gaussian$sender, sig_gaussian$receiver, sig_gaussian$metabolite_name, sep = "_")
  
  # 计算重叠
  overlap <- length(intersect(sig_knn$interaction_id, sig_gaussian$interaction_id))
  total_unique <- length(union(sig_knn$interaction_id, sig_gaussian$interaction_id))
  jaccard <- overlap / total_unique
  
  cat("\n空间方法一致性:\n")
  cat("KNN显著通讯:", nrow(sig_knn), "\n")
  cat("Gaussian显著通讯:", nrow(sig_gaussian), "\n")
  cat("重叠数:", overlap, "\n")
  cat("Jaccard指数:", round(jaccard, 3), "\n")
}

# ============================================================================
# 7. 空间距离统计
# ============================================================================
cat("\n=== 7. 空间距离统计 ===\n")

dist_stats <- getSpatialDistanceStats(obj_knn)

cat("\n空间距离统计:\n")
cat("Spot数量:", dist_stats$n_spots, "\n")
cat("最小距离:", round(dist_stats$min_distance_um, 1), "um\n")
cat("中位距离:", round(dist_stats$median_distance_um, 1), "um\n")
cat("最大距离:", round(dist_stats$max_distance_um, 1), "um\n")

# ============================================================================
# 8. 生成空间分析报告
# ============================================================================
cat("\n=== 8. 生成空间分析报告 ===\n")

spatial_report <- list(
  # 数据信息
  n_spots = dist_stats$n_spots,
  n_cell_types = length(unique(obj_knn@cell_meta[[obj_knn@cell_type_column]])),
  
  # 空间统计
  distance_stats = dist_stats,
  
  # 通讯结果
  n_sig_knn = nrow(obj_knn@significant_interactions),
  n_sig_gaussian = nrow(obj_gaussian@significant_interactions),
  
  # 方法一致性
  method_consistency = if (exists("jaccard")) jaccard else NA,
  
  # scMetaLink空间功能
  spatial_methods = c("knn", "gaussian", "exponential", "linear", "threshold"),
  scMetaLink_advantage = "原生空间转录组支持，MEBOCOST不具备",
  
  timestamp = Sys.time()
)

saveRDS(spatial_report, file.path(output_dir, "spatial_analysis_report.rds"))

# 创建Markdown报告
report_md <- sprintf(
"# scMetaLink 空间分析报告

## 数据概览
- Spot数量: %d
- 细胞类型数: %d
- 最小距离: %.1f um
- 中位距离: %.1f um

## 空间通讯结果
- KNN方法显著通讯: %d
- Gaussian方法显著通讯: %d

## scMetaLink空间分析优势

### MEBOCOST不具备的功能:
1. **空间距离加权通讯**: C_spatial = √(P × S) × w(d)
2. **多种空间加权方法**: KNN, Gaussian, Exponential, Linear, Threshold
3. **空间可视化**: 组织坐标上的通讯网络
4. **热点识别**: 识别高通讯活性区域

### 生物学意义:
- 代谢物是扩散性信号分子
- 空间距离直接影响通讯效率
- 空间分析更符合生物学现实

## 结论
scMetaLink的空间分析模块是与MEBOCOST的关键差异化功能。

---
生成时间: %s
",
  dist_stats$n_spots,
  length(unique(obj_knn@cell_meta[[obj_knn@cell_type_column]])),
  dist_stats$min_distance_um,
  dist_stats$median_distance_um,
  nrow(obj_knn@significant_interactions),
  nrow(obj_gaussian@significant_interactions),
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

writeLines(report_md, file.path(output_dir, "spatial_analysis_report.md"))

# ============================================================================
# 完成
# ============================================================================
cat("\n")
cat("=" , paste(rep("=", 60), collapse = ""), "\n")
cat("空间分析演示完成!\n")
cat("=" , paste(rep("=", 60), collapse = ""), "\n")
cat("结果保存在:", normalizePath(output_dir), "\n")
cat("\n核心发现:\n")
cat("- scMetaLink提供原生空间转录组支持\n")
cat("- 支持多种空间加权方法\n")
cat("- MEBOCOST不具备空间分析功能\n")
cat("\n")
