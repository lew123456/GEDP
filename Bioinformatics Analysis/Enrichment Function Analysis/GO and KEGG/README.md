# 基因 GO & KEGG 富集分析流程说明

## 1. 数据准备

### 数据源：

- 准备好要进行富集分析的基因。
- 数据文件需包含基因符号列（SYMBOL），文件名为 `gene.csv`。



## 2. R 环境设置

- **加载必要的 R 包**：
  
  - 包括 `AnnotationDbi`、`org.Hs.eg.db`、`clusterProfiler`、`dplyr`、`ggplot2` 等。
- **设置工作目录**：
  
  - 设置工作目录路径，例如 `F:/R`。

## 3. 数据读取与基因 ID 转换

- **读取基因数据文件**：
  
  - 读取要进行富集分析的基因数据文件 `gene.csv`。
- **基因 ID 转换**：
  
  - 使用 `bitr` 函数将基因符号（SYMBOL）转换为 Entrez ID。

## 4. GO 富集分析

- **定义 GO 富集分析函数**：
  
  - 定义一个函数 `perform_go_enrichment`，用于进行 GO 富集分析。
  - 参数包括基因列表和 GO 本体类型（BP/CC/MF）。
- **进行三个本体的 GO 富集分析**：
  
  - 分别进行生物过程（BP）、细胞组分（CC）和分子功能（MF）的 GO 富集分析。
- **保存结果**：
  
  - BP的分析结果保存至  `ego_result_BP.csv` 文件。
  - CC的分析结果保存至  `ego_result_CC.csv` 文件。
  - MF的分析结果保存至  `ego_result_MF.csv` 文件。
  - 三者分析结果的合集保存至  `ego.csv` 文件。

## 5. GO 富集分析结果可视化

- 绘制并保存<font color="red">气泡图</font>：
  
  - 为每个本体绘制气泡图，展示显著的 GO 项。
  - 生成图片 `GO_BP_bubbleplot.png `、 ` GO_CC_bubbleplot.png`、`GO_MF_bubbleplot.png ` 并保存至当前文件夹。
- 绘制并保存<font color="red">柱状图</font>：
  
  - 为每个本体绘制柱状图，展示显著的 GO 项。
  - 生成图片 `GO_BP_barplot.png `、 `GO_CC_barplot.png`、`GO_MF_barplot.png ` 并保存至当前文件夹。
- 绘制并保存<font color="red">圈(弦)图</font>：
  
  - 绘制圈(弦)图，展示基因与 GO 项的关系。
  - 生成图片 `GO_BP_circos.png `、 ` GO_CC_circos.png`、`GO_MF_circos.png ` 并保存至当前文件夹。
- 绘制并保存<font color="red">网络图</font>：
  
  - 绘制网络图，展示基因与 GO 项的关系。
  - 生成图片 `GO_BP_cnetplot.png `、 ` GO_CC_cnetplot.png`、`GO_MF_cnetplot.png ` 并保存至当前文件夹。
- 绘制并保存<font color="red">词云图</font>：
  
  - 绘制词云图，展示显著的 GO 项。
  - 生成图片 `GO_BP_wordcloud.png `、 ` GO_CC_wordcloud.png`、`GO_MF_wordcloud.png ` 并保存至当前文件夹。

## 6. KEGG 富集分析

- **定义 KEGG 富集分析函数**：
  
  - 定义一个函数 `perform_kegg_enrichment`，用于进行 KEGG 富集分析。
- **进行 KEGG 富集分析**：
  
  - 使用定义的函数进行 KEGG 富集分析。
- **保存结果**：
  
  - KEGG的分析结果保存至  `ego_result_kegg.csv` 文件。

## 7. KEGG 富集分析结果可视化

- 绘制并保存<font color="red">KEGG 通路图和通路背景图</font>：
  
  - 绘制前几个显著的 KEGG 通路图和通路背景图。
  - 生成图片 `hsa+五位数字（KEGG官方ID）.png` 并保存至当前文件夹。
- 绘制并保存<font color="red">气泡图</font>：
  
  - 绘制气泡图，展示显著的 KEGG 通路。
  - 生成图片 `kegg_bubbleplot.png` 并保存至当前文件夹。
- 绘制并保存<font color="red">柱状图</font>：
  
  - 绘制柱状图，展示显著的 KEGG 通路。
  - 生成图片 `kegg_barplot.png` 并保存至当前文件夹。
- 绘制并保存<font color="red">圈(弦)图</font>：
  
  - 绘制圈(弦)图，展示基因与 KEGG 通路的关系。
  - 生成图片 `kegg_circos.png` 并保存至当前文件夹。
- 绘制并保存<font color="red">网络图</font>：
  
  - 绘制网络图，展示基因与 KEGG 通路的关系。
  - 生成图片 `kegg_cnetplot.png` 并保存至当前文件夹。
- 绘制并保存<font color="red">词云图</font>：
  
  - 绘制词云图，展示显著的 KEGG 通路。
  - 生成图片 `kegg_wordcloud.png` 并保存至当前文件夹。


### 备注：

- 所有输出文件默认保存在当前工作目录。
- 如果需要调整图表样式或参数，请参考代码中的注释说明。
