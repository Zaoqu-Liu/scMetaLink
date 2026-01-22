# Colorectal Cancer Data (GSE146771)

## 数据来源
Lee HO, et al. Lineage-dependent gene expression programs influence the 
immune landscape of colorectal cancer. Nat Genet. 2020.

## 下载方式
### 方式1: 从GEO下载
```r
library(GEOquery)
getGEOSuppFiles('GSE146771', makeDirectory = TRUE)
```

### 方式2: 从SeuratData下载(如果有)
某些数据集可以通过SeuratData包获取

## 数据内容
- 结直肠癌肿瘤组织单细胞RNA-seq
- 包含: Tumor cells, T cells, Macrophages, Fibroblasts等
- 适合验证乳酸-TAM, 腺苷-T cell等通讯
