# 生物标志物ROC曲线验证分析
setwd("F:/R")          # 设置工作目录路径，用于指定R的工作文件夹 

# 加载所需的R包
library(pROC)
library(here)
library(tidyverse)

# 配置参数 ----------------------------------------------------------------

DATA <- "gene.csv"
ANALYSIS_GENES <- c("SP1", "FOXO1", "ESR1", "MAPK1")  # 可在此扩展分析基因
#前 15 个样本（索引 1 到 15）为一组，根据数据可修改。
CASE_NUM <- 15    
#后 8 个样本（索引 16 到 23）为一组，根据数据可修改。
CONTROL_NUM <- 8 
OUTPUT_FILE <- "roc_curves.png"

# 函数定义 ----------------------------------------------------------------

plot_single_roc <- function(data, 
                            gene_name, 
                            curve_color = "#AB274F",
                            auc_position = c(0.3, 0.4)) {
  
  # 参数验证
  stopifnot(
    gene_name %in% colnames(data),
    "label" %in% colnames(data),
    is.numeric(data[[gene_name]]),
    is.factor(data$label)
  )
  
  # 构建ROC对象
  roc_obj <- roc(response = data$label,
                 predictor = data[[gene_name]],
                 quiet = TRUE)
  
  # 绘制主曲线
  plot(roc_obj,
       main = paste("ROC Curve -", gene_name),
       col = curve_color,
       lwd = 2,
       print.auc = TRUE,
       print.auc.x = auc_position[1],
       print.auc.y = auc_position[2],
       legacy.axes = TRUE)  # 确保x轴显示1-特异性
  
  # 添加参考线
  abline(a = 0, b = 1, lty = 2, col = "gray60")
}

# 批量绘制ROC曲线
plot_multiple_roc <- function(data, genes, ncol = 2) {
  # 计算布局参数
  n <- length(genes)
  nrow <- ceiling(n / ncol)
  
  # 设置图形设备参数（结束后自动重置）
  op <- par(mfrow = c(nrow, ncol),
            mar = c(4, 4, 2.5, 1),
            family = "sans")
  on.exit(par(op))
  
  # 循环绘制
  walk(genes, ~plot_single_roc(data, .x))
}

# 动态验证数据格式
validate_data <- function(data, genes) {
  # 存在性检查
  missing_genes <- setdiff(genes, colnames(data))
  if (length(missing_genes) > 0) {
    stop("数据缺失必要基因列: ", paste(missing_genes, collapse = ", "))
  }
  
  # 类型检查
  non_numeric <- data[genes] %>% 
    select(where(Negate(is.numeric))) %>% 
    colnames()
  
  if (length(non_numeric) > 0) {
    stop("以下基因列非数值型: ", paste(non_numeric, collapse = ", "))
  }
  
  # 标签检查
  if (!"label" %in% colnames(data)) {
    stop("数据中缺少标签列 'label'")
  }
}

# 自动化数据预处理
prepare_data <- function(file_path, case_num, control_num) {
  # 安全读取
  if (!file.exists(file_path)) stop("文件不存在: ", file_path)
  
  read_csv(file_path, show_col_types = FALSE) %>% 
    mutate(
      label = factor(
        c(rep(1, case_num), rep(0, control_num)),
        levels = c(0, 1), 
        labels = c("Control", "Case")
      )
    ) %>% 
    drop_na() %>%  # 自动处理缺失值
    return()
}

# 执行流程 ----------------------------------------------------------------

analysis_data <- prepare_data(here(DATA), CASE_NUM, CONTROL_NUM)
validate_data(analysis_data, ANALYSIS_GENES)

# 绘图输出
png(OUTPUT_FILE, width = 1200, height = 900, res = 150)     #如果绘图出现乱码，修改图片的 width和height，可解决乱码
plot_multiple_roc(analysis_data, ANALYSIS_GENES, ncol = 2)  #ncol = 2 时，每行显示两条 ROC 曲线
dev.off()