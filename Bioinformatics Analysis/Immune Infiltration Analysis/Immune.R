# 设置工作目录 (可选，如果数据文件在当前工作目录下则不需要)
setwd("F:/R")          # 设置工作目录路径，用于指定R的工作文件夹 

# 加载必要的R包
library(ggpubr)    # 提供基于ggplot2的高级绘图函数，适合出版级别的图形绘制
library(ggplot2)   # 强大的数据可视化包，用于创建复杂的统计图形
library(reshape2)  # 提供数据重塑功能，如将宽格式数据转换为长格式（melt和dcast函数）
library(dplyr)     # 数据操作的核心包，提供高效的数据筛选、排序、分组和汇总功能
library(tibble)    # 提供增强版的数据框（tibble），比传统data.frame更友好和高效
library(ggsci)     # 提供科学期刊风格的颜色调色板，适用于ggplot2等绘图包
library(pheatmap)  # 绘制热图的专用包，支持聚类分析和自定义颜色
library(tidyHeatmap)   # 基于tidyverse的热图绘制工具，与tidy数据格式无缝集成
library(tidyverse)     # 包含一系列数据科学相关包（如dplyr、ggplot2等）的集合，提供统一的工作流
library(RColorBrewer)  # 提供多种预定义的颜色调色板，适用于数据可视化
library(forcats)  # 专注于因子变量的操作，提供重新排序、合并和修改因子的工具
library(corrplot)      # 用于绘制相关性矩阵的可视化工具，支持多种图形样式
library(tidyr)    # 提供数据整理功能，如处理缺失值、分离或合并列（pivot_longer/pivot_wider）
library(linkET)   # 用于基因组学数据的可视化和分析，特别是eQTL和GWAS结果
library(IOBR)     # 专注于免疫肿瘤学数据分析的工具包，提供生物信息学相关的功能

# &**&一般进行细胞免疫浸润分析对预处理后的全部基因表达数据进行,不是筛选后的数据,如果基因数量较少（例如少于 100 个）,可能会影响分析的准确性
# 1. 数据读取和预处理
# 读取想要进行细胞免疫浸润分析的基因数据文件,假设文件名为 "gene.csv",将第一行作为列名，将第一列作为行名。
#==========================================================================================================================================================
#header = TRUE (列名/表头):
#如果 CSV 文件第一行是列名,务必设置 header = TRUE,否则 R 会将列名当作数据。
#row.names = 1 (行名):
#row.names 可以是数字,表示第几列是行名；也可以是字符串,表示列名是什么。例如,如果行名在名为 "ID" 的列中,可以使用 row.names = "ID"。
#重要提示：当使用 row.names = 1 时,R 会将第一列作为行名,并且这一列将不再作为数据的一部分出现在数据框中。这意味着数据框的列数会比 CSV 文件中少一列。
#==========================================================================================================================================================
data <- read.csv("gene.csv", header = TRUE, row.names = 1)

# 2. CIBERSORT细胞免疫浸润分析分析
cibersort <- deconvo_tme(eset = data,
                         method = "cibersort",
                         arrays = FALSE, # 是否使用微阵列数据
                         perm = 500)     # 排列次数

# 3. 将分析结果保存到当前路径
cibersort_result <- as.data.frame(cibersort)
write.csv(cibersort_result,file = "cibersort.csv",row.names = T)

# 4. 绘制CIBERSORT细胞免疫浸润分析结果图
#--------------------------------------------------------------------------------------------------------------

# @绘制样本的免疫细胞比例
# 将格式变换长数据
cibersort_long<-cibersort%>%
  select(`P-value_CIBERSORT`,Correlation_CIBERSORT,RMSE_CIBERSORT,ID,everything())%>%
  pivot_longer(-c(1:4),names_to="cell_type",values_to="fraction")%>%     #变换为长数据
  dplyr::mutate(cell_type=gsub("_CIBERSORT","",cell_type),    #去除细胞名称的尾巴
                cell_type=gsub("_"," ",cell_type))            #“_”改空格，名字好看一些

# 绘制部分样本的堆叠条形图,展示出每个样本的免疫细胞比例。每一行表示一个样本,横坐标轴代表了免疫细胞的种类以及在样本中的丰度
# 示例代码是绘制前10的样本为例，如果要绘制前20,则修改：cibersort[1:20,][-c(24:26)]
png(filename = "cibersort_cell_part.png", width = 6500, height = 3000, res = 300)
cell_bar_plot(cibersort[1:10,][-c(24:26)],     # **以前10个样本为例
              pattern="CIBERSORT",
              title="CIBERSORT Cell Fraction",
              palette = "palette1")
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 绘制全部样本的堆叠条形图,展示出每个样本的免疫细胞比例。每一列表示一个样本，纵坐标轴代表了免疫细胞的种类以及在样本中的丰度
png(filename = "cibersort_cell_all.png", width = 6500, height = 3000, res = 300)
cibersort_long%>%
  ggplot(aes(ID,fraction))+
  geom_bar(stat="identity",position="stack",aes(fill=cell_type))+
  labs(x=NULL)+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=palette4,name=NULL)+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom")
dev.off()

#--------------------------------------------------------------------------------------------------------------

# @绘制分组箱线图,展示免疫细胞之间含量的比较
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 在图片上标注的 "ns"、"*"、"**" 等符号表示统计学显著性水平：
# "ns"：表示无显著性差异（p-value ≥ 0.05）。
# "*"：表示 p-value < 0.05，有显著性差异。
# "**"：表示 p-value < 0.01，有非常显著的差异。
# "***"：表示 p-value < 0.001，有极显著的差异。
# 这些标注是由 stat_compare_means 函数根据Wilcoxon秩和检验的结果自动生成的。
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cibersort_long <- cibersort[,-c(1,24:26)]      # 从数据框cibersort中移除第1列以及第24到26列,将剩余的列保存到一个新的数据框cibersort_long中
cibersort_long$ID <- rownames(cibersort_long)  # 添加样本ID列,将行名作为新的一列 ID 添加到数据框中。

# **根据样本顺序创建分组：示例为前5个为对照组，后5个为IDD组
# **重要**: 这个需要根据每组的样本数目进行修改。可以根据这一部分修改分组图例的名称，示例为control和IDD
cibersort_long$Group <- c(rep("control", 5), rep("IDD", 5))

# 转换为长格式：细胞类型和其对应的比例
cibersort_long <- cibersort_long %>%
  gather(key = "Cell_Type", value = "Proportion", -Group, -ID)

# 去掉 Cell_Type 中的 "_CIBERSORT" 后缀
cibersort_long$Cell_Type <- gsub("_CIBERSORT", "", cibersort_long$Cell_Type)

# 过滤掉比例全为0的免疫细胞列
cibersort_filtered <- cibersort_long %>%
  group_by(Cell_Type) %>%
  filter(any(Proportion > 0)) %>%  # 只保留至少有一个非零值的免疫细胞
  ungroup()

# 绘制图形并保存为PNG文件
png(filename = "boxplot.png", width = 3500, height = 1500, res = 150)
ggplot(cibersort_filtered, aes(Cell_Type, Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  scale_fill_manual(values = c("control" = "skyblue", "IDD" = "salmon")) +     # 要和上面的c(rep("control", 15), rep("IDD", 8))里面的分组名称一样,示例为control和IDD
  theme_bw() +
  labs(x = NULL, 
       y = element_text("Estimated Proportion", size = 18)) +  # 设置Y轴标题字体大小
  theme(legend.position = "top",
        legend.text = element_text(size = 20),                # 图例文字大小
        legend.title = element_text(size = 20),               # 图例标题大小
        axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # X轴标签字体大小
        axis.text.y = element_text(size = 18),                # Y轴标签字体大小
        axis.title.y = element_text(size = 18),               # Y轴标题字体大小
        panel.grid.major = element_line(color = "gray80", linetype = "dashed", size = 0.5),  # 主网格线
        panel.grid.minor = element_line(color = "gray90", linetype = "dotted", size = 0.5),  # 次网格线
        panel.border = element_rect(color = "black", fill = NA, size = 1.5)) +  # 边框粗细
  stat_compare_means(aes(group = Group, label = ..p.signif..),
                     method = "wilcox.test", 
                     label.y = 0.5, 
                     size = 8)  # 调整显著性标注的字体大小
dev.off()

#--------------------------------------------------------------------------------------------------------------

# @绘制免疫细胞含量之间的相关性图
# 数据预处理
data_for_cor <- cibersort[,-c(1,24:26)]

# 检查并移除标准差为0的列
sd_zero <- sapply(data_for_cor, sd, na.rm = TRUE) == 0
if(any(sd_zero)) {
  data_for_cor <- data_for_cor[, !sd_zero]
  message("移除了标准差为0的列：", paste(names(data_for_cor)[sd_zero], collapse = ", "))
}

# 移除含有NA的行
data_for_cor <- na.omit(data_for_cor)

# 检查是否有无穷值并替换
data_for_cor[is.infinite(as.matrix(data_for_cor))] <- NA
data_for_cor <- na.omit(data_for_cor)

# 确保数据有效性
if(ncol(data_for_cor) < 2 || nrow(data_for_cor) < 3) {
  stop("数据预处理后样本量不足，无法进行相关性分析")
}

# 计算相关性
cor_cibersort <- cor(data_for_cor, use = "pairwise.complete.obs")
colnames(cor_cibersort) <- gsub("_CIBERSORT", "", colnames(cor_cibersort))
rownames(cor_cibersort) <- gsub("_CIBERSORT", "", rownames(cor_cibersort))

# 检查相关性矩阵是否有效
if(any(is.na(cor_cibersort)) || any(is.infinite(as.matrix(cor_cibersort)))) {
  stop("相关性矩阵包含NA或无穷值，请检查数据")
}

# 计算p值矩阵
pmat <- tryCatch({
  cor.mtest(data_for_cor, conf.level = 0.9)
}, error = function(e) {
  message("计算p值时出错：", e$message)
  return(NULL)
})

# 绘制相关性热图
png(filename = "cell_correlation.png", width = 2200, height = 1800, res = 300)
tryCatch({
  corrplot(cor_cibersort,
           method = "circle",
           type = "lower",
           tl.col = "black",
           tl.cex = 0.73,
           tl.srt = 20,
           col = COL2('RdYlBu', 9), # 使用 Set1 调色板的9种颜色
           order = "hclust",  # 改用层次聚类排序，避免特征向量计算的问题
           mar = c(0,0,3,0),
           p.mat = if(!is.null(pmat)) pmat$p else NULL,
           sig.level = 0.1,
           insig = "blank",
           title = " ")
}, error = function(e) {
  message("绘制相关性热图时出错：", e$message)
})
dev.off()

#--------------------------------------------------------------------------------------------------------------

# @绘制热图,展示待选基因和所有细胞的关系（多对多）
# 确保genes包含在数据行名中,这里示例4个基因
genes <- c("FOXO1","ESR1","MYD88","MAPK1")
genes_exp <- as.data.frame(t(data[rownames(data) %in% genes, ]))    # 提取基因的表达数据
genes_exp <- genes_exp[match(cibersort$ID, rownames(genes_exp)), ]  # 确保与 cibersort 数据的样本ID一致

cibersort_filtered <- cibersort[, -c(1, 24:26)]  # 去除第一列（ID）和最后几列（不相关的列）

# 过滤掉所有比例全为零的免疫细胞列
cibersort_filtered <- cibersort_filtered[, colSums(cibersort_filtered != 0) > 0]

# 执行相关性分析
cor_res <- correlate(genes_exp, cibersort_filtered, method = "spearman")  #计算相关性系数

# 继续进行后续的数据整合和绘图部分
df_r <- cor_res$r %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(-1, names_to = "cell_type", values_to = "correlation")

df_p <- cor_res$p %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(-1, names_to = "cell_type", values_to = "pvalue")

df_cor <- df_r %>%
  left_join(df_p) %>%
  mutate(stars = cut(pvalue, breaks = c(-Inf, 0.05, 0.01, 0.001, Inf), right = F, labels = c("***", "**", "*", " ")))

df_cor$cell_type <- gsub('_CIBERSORT', '', df_cor$cell_type) # 移除细胞类型名称的后缀

# 绘制热图
png(filename = "correlation_part.png", width = 3000, height = 1200, res = 300)
ggplot(df_cor, aes(cell_type, gene)) +
  geom_tile(aes(fill = correlation)) +  # 绘制热图
  geom_text(aes(label = stars), color = "black", size = 8) +     # 图内的字体设置
  scale_fill_gradient2(low = '#318CE7', high = 'red', mid = 'white',
                       limit = c(-1, 1), name = paste0("*p<0.05", "\n\n", "**p<0.01", "\n\n", "***p<0.001", "\n\n", "Correlation")) +  # 添加显著性检验结果
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.y = element_blank(),
        panel.background = element_blank())
dev.off()

#--------------------------------------------------------------------------------------------------------------

# @绘制免疫细胞分组丰度热图
cibersort_long <- cibersort[,-c(1,24:26)]      # 从数据框cibersort中移除第1列以及第24到26列,将剩余的列保存到一个新的数据框cibersort_long中
cibersort_long$ID <- rownames(cibersort_long)  # 添加样本ID列,将行名作为新的一列 ID 添加到数据框中。

# **根据样本顺序创建分组：示例为前5个为对照组，后5个为IDD组
# **重要**: 这个需要根据每组的样本数目进行修改。可以根据这一部分修改分组图例的名称，示例为control和IDD
cibersort_long$Group <- c(rep("control", 5), rep("IDD", 5))

# 转换为长格式：细胞类型和其对应的比例
cibersort_long <- cibersort_long %>%
  gather(key = "Cell_Type", value = "Proportion", -Group, -ID)

# 去掉 Cell_Type 中的 "_CIBERSORT" 后缀
cibersort_long$Cell_Type <- gsub("_CIBERSORT", "", cibersort_long$Cell_Type)

# 过滤掉比例全为0的免疫细胞列
cibersort_filtered <- cibersort_long %>%
  group_by(Cell_Type) %>%
  filter(any(Proportion > 0)) %>%  # 只保留至少有一个非零值的免疫细胞
  ungroup()

png(filename = "heatmap.png", width = 1800, height = 1200, res = 300)
cibersort_filtered %>%
  group_by(Group) %>%
  heatmap(
    .row = Cell_Type,               # 行为细胞类型
    .column = ID,                   # 列为样本ID
    .value = Proportion,            # 值为比例
    scale = "column",               # 对列进行缩放（标准化）
    
    # 设置热图颜色调色板
    palette_value = circlize::colorRamp2(
      seq(-2, 2, length.out = 11),  # 颜色范围从-2到2，分为11个等级
      RColorBrewer::brewer.pal(11, "Spectral")  # 使用Spectral调色板
    ),
    
    # 设置分组标签的颜色
    palette_grouping = list(c("#1F78B4", "#E31A1C")),  # "normal" 为蓝色，"tumor" 为红色
    
    show_column_names = FALSE,  # 不显示列名
    row_names_gp = gpar(fontsize = 10),  # 行名字体大小
    column_title_gp = gpar(fontsize = 7),  # 列标题字体大小
    row_title_gp = gpar(fontsize = 7)      # 行标题字体大小
  )
dev.off()

#--------------------------------------------------------------------------------------------------------------

# @绘制箱线图，不分组，其中，每一列表示一个免疫细胞种类，纵坐标轴表示该细胞在样本中的丰度
# 将格式变换长数据
cibersort_long<-cibersort%>%
  select(`P-value_CIBERSORT`,Correlation_CIBERSORT,RMSE_CIBERSORT,ID,everything())%>%
  pivot_longer(-c(1:4),names_to="cell_type",values_to="fraction")%>%     #变换为长数据
  dplyr::mutate(cell_type=gsub("_CIBERSORT","",cell_type),    #去除细胞名称的尾巴
                cell_type=gsub("_"," ",cell_type))            #“_”改空格，名字好看一些

png(filename = "boxplot_nogroup.png", width = 3000, height = 1800, res = 300)
ggplot(cibersort_long,aes(fct_reorder(cell_type,fraction),fraction,fill=cell_type))+
  geom_boxplot()+
  theme_bw()+
  labs(x="Cell Type",y="Estimated Proportion")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom")+
  scale_fill_manual(values=palette4)
dev.off()

