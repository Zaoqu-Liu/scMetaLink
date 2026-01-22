# Breast Cancer Data (GSE75688)

## 数据来源
Chung W, et al. Single-cell RNA-seq enables comprehensive tumour and immune 
cell profiling in primary breast cancer. Nat Commun. 2017.

## 下载
```r
library(GEOquery)
getGEOSuppFiles('GSE75688', makeDirectory = TRUE)
```

## 数据特点
- 原发乳腺癌单细胞测序
- 包含tumor cells, immune cells
- 可以分析肿瘤代谢与免疫抑制
