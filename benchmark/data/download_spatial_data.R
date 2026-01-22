#!/usr/bin/env Rscript
# ============================================================================
# 下载空间转录组数据
# 在RStudio中运行此脚本
# ============================================================================

library(TENxVisiumData)
library(SpatialExperiment)

# 设置保存目录
save_dir <- "spatial"
if (!dir.exists(save_dir)) dir.create(save_dir)

# ============================================================================
# 下载主要的肿瘤空间数据
# ============================================================================

cat("\n=== 下载空间转录组数据 ===\n")

# 1. 乳腺癌 (浸润性导管癌) - 推荐用于benchmark
cat("\n1. 下载乳腺癌数据 (HumanBreastCancerIDC)...\n")
tryCatch({
  spe_breast <- HumanBreastCancerIDC()
  saveRDS(spe_breast, file.path(save_dir, "breast_cancer_visium.rds"))
  cat("   已保存:", file.path(save_dir, "breast_cancer_visium.rds"), "\n")
  cat("   Spots:", ncol(spe_breast), "Genes:", nrow(spe_breast), "\n")
}, error = function(e) cat("   下载失败:", e$message, "\n"))

# 2. 结直肠癌 - 用于验证肿瘤免疫通讯
cat("\n2. 下载结直肠癌数据 (HumanColorectalCancer)...\n")
tryCatch({
  spe_crc <- HumanColorectalCancer()
  saveRDS(spe_crc, file.path(save_dir, "colorectal_cancer_visium.rds"))
  cat("   已保存:", file.path(save_dir, "colorectal_cancer_visium.rds"), "\n")
  cat("   Spots:", ncol(spe_crc), "Genes:", nrow(spe_crc), "\n")
}, error = function(e) cat("   下载失败:", e$message, "\n"))

# 3. 胶质母细胞瘤 - 用于验证肿瘤代谢通讯
cat("\n3. 下载胶质母细胞瘤数据 (HumanGlioblastoma)...\n")
tryCatch({
  spe_gbm <- HumanGlioblastoma()
  saveRDS(spe_gbm, file.path(save_dir, "glioblastoma_visium.rds"))
  cat("   已保存:", file.path(save_dir, "glioblastoma_visium.rds"), "\n")
  cat("   Spots:", ncol(spe_gbm), "Genes:", nrow(spe_gbm), "\n")
}, error = function(e) cat("   下载失败:", e$message, "\n"))

# 4. 卵巢癌 - 有免疫panel
cat("\n4. 下载卵巢癌数据 (HumanOvarianCancer)...\n")
tryCatch({
  spe_ov <- HumanOvarianCancer()
  saveRDS(spe_ov, file.path(save_dir, "ovarian_cancer_visium.rds"))
  cat("   已保存:", file.path(save_dir, "ovarian_cancer_visium.rds"), "\n")
  cat("   Spots:", ncol(spe_ov), "Genes:", nrow(spe_ov), "\n")
}, error = function(e) cat("   下载失败:", e$message, "\n"))

cat("\n=== 下载完成 ===\n")
cat("数据保存在:", normalizePath(save_dir), "\n")

# 列出下载的文件
cat("\n已下载的文件:\n")
files <- list.files(save_dir, pattern = "\\.rds$", full.names = TRUE)
for (f in files) {
  info <- file.info(f)
  cat(sprintf("  %s (%.1f MB)\n", basename(f), info$size / 1e6))
}
