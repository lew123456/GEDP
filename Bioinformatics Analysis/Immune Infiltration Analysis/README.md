# 免疫细胞浸润分析流程说明

## 1. 数据准备

### 数据源：

- 准备好要进行免疫细胞浸润分析的基因。（基因数目一般要足够多）

## 2. R环境设置

### 加载必要的R包：

包括`ggpubr`、`ggplot2`、`reshape2`、`dplyr`、`tibble`、`ggsci`、`pheatmap`、`tidyHeatmap`、`tidyverse`、`RColorBrewer`、`forcats`、`corrplot`、`tidyr`、`linkET`、`IOBR`、`AnnotationDbi`、`org.Hs.eg.db`、`clusterProfiler`等。

### 设置工作目录：

- 设置工作目录路径，例如`F:/R`。

## 3. 数据读取与预处理

### 读取基因数据文件：

- 读取要进行分析的基因数据文件`gene.csv`。
- 设置`header = TRUE`以将第一行作为列名。
- 设置`row.names = 1`以将第一列作为行名。

## 4. 免疫细胞浸润分析

### CIBERSORT分析：

- 使用`deconvo_tme`函数进行CIBERSORT细胞免疫浸润分析。
  - 参数包括`eset`（基因表达数据）、`method`（CIBERSORT方法）、`arrays`（是否使用微阵列数据）和`perm`（排列次数）。

### 保存分析结果：

- 将分析结果保存为`cibersort.csv`文件。

## 5. 免疫细胞浸润分析结果可视化

### 免疫细胞比例图绘制

1. **部分样本的免疫细胞比例图**
   
   - 将分析结果转换为长格式数据。
   - 使用`cell_bar_plot`函数绘制堆叠条形图，展示每个样本的免疫细胞比例。
   - 保存图像为`cibersort_cell_part.png`。
2. **全部样本的免疫细胞比例图**
   
   - 使用`ggplot2`绘制堆叠条形图，展示每个样本的免疫细胞比例。
   - 保存图像为`cibersort_cell_all.png`。

### 绘制免疫细胞比较箱线图：

- 根据样本创建分组（如对照组和IDD组）。
- 使用`ggplot2`绘制箱线图，展示免疫细胞之间含量的比较。
- 添加显著性标注。
- 保存图像为`boxplot.png`。

### 绘制免疫细胞含量相关性图：

- 计算相关性矩阵。
- 使用`corrplot`绘制相关性热图。
- 保存图像为`cell_correlation.png`。

### 绘制基因与免疫细胞关系热图：

- 提取特定基因的表达数据。
- 进行相关性分析。
- 使用`ggplot2`绘制热图，展示特定基因与免疫细胞的关系。
- 保存图像为`correlation_part.png`。

### 绘制免疫细胞分组丰度热图：

- 根据样本创建分组（如对照组和IDD组）。
- 使用`heatmap`绘制热图，展示细胞在样本组织中的丰度。
- 保存图像为`heatmap.png`。

### 绘制免疫细胞丰度箱线图（不分组）：

- 将分析结果转换为长格式数据。
- 使用`ggplot2`绘制箱线图，展示细胞在样本中的丰度。
- 保存图像为`boxplot_nogroup.png`。


### 备注：

- 所有输出文件默认保存在当前工作目录。
- 如果需要调整图表样式或参数，请参考代码中的注释说明。
