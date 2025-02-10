# 设置工作目录 (可选，如果数据文件在当前工作目录下则不需要)
setwd("F:/R") 

# 加载必要的R包
library(pheatmap)

# 1. 数据读取和预处理
# 读取想要绘制热图的基因数据文件，假设文件名为 "DEGs_intersect.csv",将第一行作为列名，将第一列作为行名。
#==========================================================================================================================================================
#header = TRUE (列名/表头):
#如果 CSV 文件第一行是列名,务必设置 header = TRUE,否则 R 会将列名当作数据。
#row.names = 1 (行名):
#row.names 可以是数字,表示第几列是行名；也可以是字符串,表示列名是什么。例如,如果行名在名为 "ID" 的列中,可以使用 row.names = "ID"。
#重要提示：当使用 row.names = 1 时,R 会将第一列作为行名,并且这一列将不再作为数据的一部分出现在数据框中。这意味着数据框的列数会比 CSV 文件中少一列。
#==========================================================================================================================================================
data_matrix <- read.csv("DEGs_intersect.csv", header = TRUE, row.names = 1)

# 转置数据矩阵
data_matrix <- t(data_matrix)

# 创建样本分组信息,假设前 5 个样本是 Group1，后 5 个样本是 Group2，可以根据这一部分修改分组图例的名称，示例为Group 1和Group 2
#sample_groups <- c(rep("Group1", 5), rep("Group2", 5))
# 定义样本分组信息（按列名匹配的为一组，剩下的为一组）
group1_samples <- c("GSM1354764","GSM1354765","GSM1354766","GSM1354767","GSM1354768")
sample_groups <- factor(ifelse(colnames(expression_data) %in% group1_samples, "Group1", "Group2"))

# 2. 绘制热力图
# 创建样本分组信息的 data.frame，并确保行名匹配
annotation_row <- data.frame(Group = sample_groups, row.names = rownames(data_matrix))
# 指定行注释信息的颜色和名称，要和上面的c(rep("Group1", 5), rep("Group2", 5))里面的分组名称一样,示例为Group 1和Group 2
annotation_colors <- list(Group = c(Group1 = "#9D1871", Group2 = "#5572CC"))

# 绘制热图并保存为图片文件（例如，PNG格式）
png("Gene_Cross_heatmap.png", width = 25, height = 9.5, units = "in", res = 500)

pheatmap(data_matrix,
         color = colorRampPalette(c("#3B9AAF", "white", "#E0531E"))(100),
         cluster_rows = TRUE,   #基于基因表达数据对样本进行聚类。
         cluster_cols = TRUE,   #基于基因表达数据对基因进行聚类。
         show_rownames = FALSE, #控制是否显示行名（样本名）。
         show_colnames = TRUE,  #控制是否显示列名（基因名）。
         annotation_row = annotation_row, #为热图的每一行添加注释信息,注释信息是样本分组。
         annotation_colors = annotation_colors, #指定行注释信息的颜色。
         
         # 调整聚类树线条粗细
         treeheight_row = 34,    #  调整行聚类树的高度
         treeheight_col = 34,    #  调整列聚类树的高度
         
         # 修改聚类树线条粗细
         gaps_row = NULL, # 设置聚类树的间隔
         gaps_col = NULL, # 设置聚类树的间隔
         
         # 增加聚类树美观度
         clustering_distance_rows = "euclidean", # 设置行聚类的距离方法
         clustering_distance_cols = "euclidean", # 设置列聚类的距离方法
         clustering_method = "ward.D2", # 设置聚类方法（例如，ward.D2）
         
         
         #调整边距和标签
         fontsize_row = 13.5,
         fontsize_col = 13.5,
         
         # 避免行名或列名太长导致重叠
         # cutree_rows = 5 #将行聚类树分成5组，根据需要调整
         # cutree_cols = 3 #将列聚类树分成3组，根据需要调整
         # 或是增加画布大小，如上述png设置
         
)

dev.off()
