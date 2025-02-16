setwd("F:/R/GEO1") 
# 加载必要的R包
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)

# 读取上传的CSV文件
data <- read.csv("GSE111016.csv")

# 假设ENTREZID列包含以ENSG00000000xxx为形式的基因ID
gene_ids <- data$ENTREZID

# 进行ID转换：从Ensembl ID（ENSEMBL）转换为基因符号（SYMBOL）
gene.df <- bitr(gene_ids, 
                fromType = "ENSEMBL", 
                toType = "SYMBOL", 
                OrgDb = org.Hs.eg.db)

# 查看gene.df的列名并确保正确
colnames(gene.df)

# 修改列名：将ENSEMBL列改为ENTREZID，以便与原始数据合并
colnames(gene.df)[colnames(gene.df) == "ENSEMBL"] <- "ENTREZID"

# 合并数据：将原始数据与转换后的基因符号数据合并
data_with_symbols <- merge(data, gene.df, by = "ENTREZID", all.x = TRUE)

# 保存转换后的数据
write.csv(data_with_symbols, "converted_gene_symbols.csv", row.names = FALSE)
