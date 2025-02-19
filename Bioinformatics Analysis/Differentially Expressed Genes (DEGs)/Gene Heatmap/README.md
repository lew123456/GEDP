# 热力图绘制

## 1. 数据准备
- **数据源**：想要绘制热图的基因数据文件（示例：`GSE56081`差异分析得到的基因与`CellAge`数据库取交集）。
- **读取数据**：读取交叉后的基因表达数据 `DEGs_intersect.csv`。
- **转置数据**：将读取的数据进行转置数据矩阵。
- **创建分组**：创建样本分组信息。

---

## 2. 可视化
- **图表类型**：使用 `pheatmap` 包，绘制热力图（heatmap）。
- **绘图内容**：基于差异基因交集数据绘制聚类热图，行列分别按欧氏距离和Ward.D2方法聚类。
- **输出文件**：生成图片 `Gene_Cross_heatmap.png` 并保存至当前文件夹。


---

**备注**：所有输出文件默认保存在当前工作目录。