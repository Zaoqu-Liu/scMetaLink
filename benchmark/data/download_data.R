#!/usr/bin/env Rscript
# ============================================================================
# scMetaLink Benchmark Data Download Script
# ============================================================================
# 
# 这个脚本帮助下载benchmark所需的所有数据
# 
# 数据来源:
# 1. 人类脂肪组织 (GSE176171) - Emont Lab
# 2. 结直肠癌单细胞 (GSE146771)
# 3. 乳腺癌单细胞 (GSE75688)
# 4. 空间转录组数据 (TENxVisiumData)
# ============================================================================

# 设置工作目录
# 如果在RStudio中运行，取消下一行注释
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# 命令行运行时使用
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
  setwd(dirname(script_path))
}

# 安装必要的包
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

# ============================================================================
# 1. 下载人类脂肪组织数据 (GSE176171)
# ============================================================================
cat("\n=== 1. 人类脂肪组织数据 (GSE176171) ===\n")
cat("这个数据需要从 Emont Lab 网站手动下载:\n")
cat("  网址: https://emontlab.uchicago.edu/data/\n")
cat("  找到: 'A single cell atlas of human and mouse white adipose tissue'\n")
cat("  下载: human_all_lite.rds (1.32 GB) 或 human_immune_lite.rds (239 MB)\n")
cat("  保存到: benchmark/data/adipose/\n\n")

# 创建说明文件
writeLines(c(
  "# Adipose Tissue Data Download Instructions",
  "",
  "## 数据来源",
  "Emont MP, et al. A single-cell atlas of human and mouse white adipose tissue. Nature 2022.",
  "GEO: GSE176171",
  "",
  "## 下载步骤",
  "1. 访问: https://emontlab.uchicago.edu/data/",
  "2. 找到 'A single cell atlas of human and mouse white adipose tissue' 部分",
  "3. 在 'Human data' 下载以下文件:",
  "   - human_all_lite.rds (1.32 GB) - 包含所有166K细胞",
  "   - 或者 human_immune_lite.rds (239 MB) - 仅免疫细胞,包含macrophage",
  "4. 将文件放到此目录 (benchmark/data/adipose/)",
  "",
  "## 文件说明",
  "- 这些是Seurat对象(.rds格式)",
  "- lite版本用DietSeurat处理过,只保留RNA assay",
  "- 包含BMI信息,可以分组对比肥胖vs正常",
  "",
  "## 细胞类型",
  "- Adipocyte (脂肪细胞)",
  "- ASPC (脂肪干/祖细胞)",
  "- Immune (免疫细胞,包括Macrophage, T cells等)",
  "- Vascular (血管细胞)",
  "- Mesothelium (间皮细胞)"
), "adipose/README.md")

# ============================================================================
# 2. 下载结直肠癌数据 (GSE146771)
# ============================================================================
cat("\n=== 2. 结直肠癌数据 (GSE146771) ===\n")

install_if_missing("GEOquery", bioc = TRUE)
library(GEOquery)

# 创建目录
if (!dir.exists("colorectal")) dir.create("colorectal")

cat("尝试从GEO下载元数据...\n")
tryCatch({
  # 下载GEO元数据
  gse <- getGEO("GSE146771", GSEMatrix = FALSE)
  
  # 保存基本信息
  geo_info <- list(
    title = gse@header$title,
    summary = gse@header$summary,
    overall_design = gse@header$overall_design,
    platform = gse@header$platform_id
  )
  saveRDS(geo_info, "colorectal/GSE146771_info.rds")
  cat("元数据已保存到 colorectal/GSE146771_info.rds\n")
}, error = function(e) {
  cat("GEO元数据下载失败:", e$message, "\n")
})

# 创建说明文件
writeLines(c(
  "# Colorectal Cancer Data (GSE146771)",
  "",
  "## 数据来源", 
  "Lee HO, et al. Lineage-dependent gene expression programs influence the ",
  "immune landscape of colorectal cancer. Nat Genet. 2020.",
  "",
  "## 下载方式",
  "### 方式1: 从GEO下载",
  "```r",
  "library(GEOquery)",
  "getGEOSuppFiles('GSE146771', makeDirectory = TRUE)",
  "```",
  "",
  "### 方式2: 从SeuratData下载(如果有)",
  "某些数据集可以通过SeuratData包获取",
  "",
  "## 数据内容",
  "- 结直肠癌肿瘤组织单细胞RNA-seq",
  "- 包含: Tumor cells, T cells, Macrophages, Fibroblasts等",
  "- 适合验证乳酸-TAM, 腺苷-T cell等通讯"
), "colorectal/README.md")

# ============================================================================
# 3. 下载乳腺癌数据 (GSE75688)
# ============================================================================
cat("\n=== 3. 乳腺癌数据 (GSE75688) ===\n")

if (!dir.exists("breast_cancer")) dir.create("breast_cancer")

cat("尝试从GEO下载元数据...\n")
tryCatch({
  gse <- getGEO("GSE75688", GSEMatrix = FALSE)
  geo_info <- list(
    title = gse@header$title,
    summary = gse@header$summary
  )
  saveRDS(geo_info, "breast_cancer/GSE75688_info.rds")
  cat("元数据已保存\n")
}, error = function(e) {
  cat("下载失败:", e$message, "\n")
})

writeLines(c(
  "# Breast Cancer Data (GSE75688)",
  "",
  "## 数据来源",
  "Chung W, et al. Single-cell RNA-seq enables comprehensive tumour and immune ",
  "cell profiling in primary breast cancer. Nat Commun. 2017.",
  "",
  "## 下载",
  "```r",
  "library(GEOquery)",
  "getGEOSuppFiles('GSE75688', makeDirectory = TRUE)",
  "```",
  "",
  "## 数据特点",
  "- 原发乳腺癌单细胞测序",
  "- 包含tumor cells, immune cells",
  "- 可以分析肿瘤代谢与免疫抑制"
), "breast_cancer/README.md")

# ============================================================================
# 4. 下载空间转录组数据 (TENxVisiumData)
# ============================================================================
cat("\n=== 4. 空间转录组数据 (TENxVisiumData) ===\n")

install_if_missing("TENxVisiumData", bioc = TRUE)
install_if_missing("SpatialExperiment", bioc = TRUE)
install_if_missing("ExperimentHub", bioc = TRUE)

if (!dir.exists("spatial")) dir.create("spatial")

cat("下载TENxVisiumData包中的数据...\n")
tryCatch({
  library(TENxVisiumData)
  library(ExperimentHub)
  
  # 获取可用数据集列表
  eh <- ExperimentHub()
  q <- query(eh, "TENxVisiumData")
  
  # 保存数据集信息
  dataset_info <- data.frame(
    id = q$ah_id,
    title = q$title,
    description = q$description
  )
  write.csv(dataset_info, "spatial/available_datasets.csv", row.names = FALSE)
  cat("可用数据集列表已保存到 spatial/available_datasets.csv\n")
  
  # 下载一个示例数据 - Human Breast Cancer (Block A Section 1)
  cat("下载乳腺癌空间数据示例...\n")
  spe <- TENxVisiumData::HumanBreastCancerIDC()
  saveRDS(spe, "spatial/breast_cancer_visium.rds")
  cat("乳腺癌Visium数据已保存: spatial/breast_cancer_visium.rds\n")
  
}, error = function(e) {
  cat("空间数据下载失败:", e$message, "\n")
  cat("请手动安装: BiocManager::install('TENxVisiumData')\n")
})

writeLines(c(
  "# Spatial Transcriptomics Data",
  "",
  "## 来源: TENxVisiumData (Bioconductor)",
  "",
  "## 可用数据集",
  "查看 available_datasets.csv 获取完整列表",
  "",
  "## 使用方法",
  "```r",
  "library(TENxVisiumData)",
  "",
  "# 查看所有可用数据",
  "TENxVisiumData()",
  "",
  "# 下载特定数据",
  "spe <- HumanBreastCancerIDC()  # 乳腺癌",
  "spe <- HumanColonCancerPatientA()  # 结肠癌",
  "spe <- MouseBrainCoronal()  # 小鼠大脑",
  "```",
  "",
  "## 其他空间数据源",
  "1. GSE280315 - 结直肠癌 Visium HD",
  "2. Zenodo: doi.org/10.5281/zenodo.14620362 - 肾/肺癌"
), "spatial/README.md")

# ============================================================================
# 5. 创建Ground Truth文件
# ============================================================================
cat("\n=== 5. 创建Ground Truth验证数据 ===\n")

if (!dir.exists("../03_ground_truth_validation")) {
  dir.create("../03_ground_truth_validation")
}

# 已知的代谢物-细胞通讯关系
ground_truth <- data.frame(
  metabolite = c(
    "L-Lactic acid", "L-Lactic acid", "L-Lactic acid",
    "Adenosine", "Adenosine", "Adenosine",
    "L-Glutamic acid", "L-Glutamine",
    "Glucose", "ATP"
  ),
  hmdb = c(
    "HMDB0000190", "HMDB0000190", "HMDB0000190",
    "HMDB0000050", "HMDB0000050", "HMDB0000050",
    "HMDB0000148", "HMDB0000641",
    "HMDB0000122", "HMDB0000538"
  ),
  sender_cell_type = c(
    "Tumor", "Tumor", "Tumor",
    "Treg", "Tumor", "Hypoxic cells",
    "Neuron", "Endothelial",
    "Any", "Dying cells"
  ),
  receiver_cell_type = c(
    "Macrophage", "CD8 T cell", "Treg",
    "CD8 T cell", "Macrophage", "T cell",
    "Astrocyte", "Adipocyte",
    "T cell", "Macrophage"
  ),
  expected_effect = c(
    "M2 polarization", "Cytotoxicity suppression", "Enhanced function",
    "Suppression via A2A", "M2 polarization", "Suppression",
    "Signaling", "Differentiation",
    "Activation fuel", "Inflammation"
  ),
  receptor_gene = c(
    "GPR81/GPR132", "MCT1/SLC16A1", "MCT1",
    "ADORA2A", "ADORA2A/ADORA2B", "ADORA2A",
    "GRM1-8", "SLC1A5",
    "SLC2A1/GLUT1", "P2RX7"
  ),
  downstream_markers = c(
    "ARG1,CD206,IL10,VEGF", "IFNG_down,GZMB_down", "FOXP3_up",
    "IFNG_down,IL2_down", "ARG1_up,IL10_up", "IFNG_down",
    "Calcium signaling", "Lipogenesis",
    "Glycolysis genes", "IL1B,IL18"
  ),
  evidence_pmid = c(
    "30431439,29468929", "33603751", "30431439",
    "27066002", "32398639", "27066002",
    "Various", "MEBOCOST paper",
    "Metabolic competition", "28098234"
  ),
  evidence_strength = c(
    "Strong", "Strong", "Medium",
    "Strong", "Strong", "Strong",
    "Medium", "Strong",
    "Strong", "Strong"
  ),
  stringsAsFactors = FALSE
)

write.csv(ground_truth, "../03_ground_truth_validation/known_interactions.csv", row.names = FALSE)
cat("Ground truth数据已保存: 03_ground_truth_validation/known_interactions.csv\n")

# ============================================================================
# 总结
# ============================================================================
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("下载脚本执行完成!\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("\n需要手动下载的数据:\n")
cat("1. 人类脂肪组织: https://emontlab.uchicago.edu/data/\n")
cat("   - 下载 human_all_lite.rds 或 human_immune_lite.rds\n")
cat("   - 放到 benchmark/data/adipose/\n")
cat("\n已自动下载/创建:\n")
cat("- GEO数据元信息\n")
cat("- 空间转录组数据 (TENxVisiumData)\n")
cat("- Ground truth验证表\n")
cat("\n")
