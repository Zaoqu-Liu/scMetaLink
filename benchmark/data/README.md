# Benchmark 数据目录

## 📁 目录结构

```
benchmark/data/
├── adipose/          # 人类脂肪组织 (GSE176171) - 需手动下载
├── breast_cancer/    # 乳腺癌单细胞
├── colorectal/       # 结直肠癌单细胞
└── spatial/          # 空间转录组数据
```

## 🔴 需要手动下载的数据

### 1. 人类脂肪组织数据 (MEBOCOST论文复现)

**来源**: GSE176171 - Emont et al. Nature 2022

**下载步骤**:
1. 访问: https://emontlab.uchicago.edu/data/
2. 找到 "A single cell atlas of human and mouse white adipose tissue"
3. 在 "Human data" 部分下载:
   - `human_all_lite.rds` (1.32 GB) - 包含所有166K细胞
   - 或 `human_immune_lite.rds` (239 MB) - 仅免疫细胞
4. 将文件放到 `benchmark/data/adipose/` 目录

**用途**: 复现MEBOCOST论文的分析，对比肥胖vs正常的代谢物通讯

---

## 🟢 自动下载的数据

### 2. GEO元数据
- `colorectal/GSE146771_info.rds` - 结直肠癌数据信息
- `breast_cancer/GSE75688_info.rds` - 乳腺癌数据信息

### 3. 空间转录组数据

在RStudio中运行 `download_spatial_data.R` 下载:

```r
setwd("benchmark/data")
source("download_spatial_data.R")
```

将下载以下数据:
- `spatial/breast_cancer_visium.rds` - 乳腺癌Visium
- `spatial/colorectal_cancer_visium.rds` - 结直肠癌Visium  
- `spatial/glioblastoma_visium.rds` - 胶质母细胞瘤Visium
- `spatial/ovarian_cancer_visium.rds` - 卵巢癌Visium

---

## 📊 数据用途

| 数据集 | 用途 | 验证目标 |
|--------|------|----------|
| 脂肪组织 | 复现MEBOCOST | 与MEBOCOST结果对比 |
| 结直肠癌 | Ground truth | 乳酸→TAM, 腺苷→T cell |
| 乳腺癌 | Ground truth | 肿瘤代谢免疫抑制 |
| 空间数据 | 空间验证 | DES (Distance Enrichment Score) |

---

## 🧬 Ground Truth

已知的代谢物-细胞通讯关系保存在:
`benchmark/03_ground_truth_validation/known_interactions.csv`

包含:
- L-Lactic acid → Macrophage (M2极化)
- L-Lactic acid → CD8 T cell (抑制)
- Adenosine → T cell (抑制)
- 等20+已验证的通讯关系

---

## 🚀 下一步

1. 下载脂肪组织数据 (手动)
2. 运行 `download_spatial_data.R` 获取空间数据
3. 开始benchmark: `benchmark/02_algorithm_comparison/`
