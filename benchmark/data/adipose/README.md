# Adipose Tissue Data Download Instructions

## 数据来源
Emont MP, et al. A single-cell atlas of human and mouse white adipose tissue. Nature 2022.
GEO: GSE176171

## 下载步骤
1. 访问: https://emontlab.uchicago.edu/data/
2. 找到 'A single cell atlas of human and mouse white adipose tissue' 部分
3. 在 'Human data' 下载以下文件:
   - human_all_lite.rds (1.32 GB) - 包含所有166K细胞
   - 或者 human_immune_lite.rds (239 MB) - 仅免疫细胞,包含macrophage
4. 将文件放到此目录 (benchmark/data/adipose/)

## 文件说明
- 这些是Seurat对象(.rds格式)
- lite版本用DietSeurat处理过,只保留RNA assay
- 包含BMI信息,可以分组对比肥胖vs正常

## 细胞类型
- Adipocyte (脂肪细胞)
- ASPC (脂肪干/祖细胞)
- Immune (免疫细胞,包括Macrophage, T cells等)
- Vascular (血管细胞)
- Mesothelium (间皮细胞)
