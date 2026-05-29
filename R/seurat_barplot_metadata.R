library(Seurat)
library(ggplot2)
library(dplyr)

#' 动态绘制 Seurat 对象元数据两列组合的细胞数目柱状图
#'
#' @param seurat_obj 一个标准的 Seurat 对象。
#' @param x_var 字符串，指定作为横轴 (X轴) 的 meta.data 列名（如 "predicted.cell_subclass"）。
#' @param fill_var 字符串，指定作为柱状图拆分/颜色填充的 meta.data 列名（如 "Condition"）。
#' @param save_path 字符串，图片保存的路径和文件名，默认为 "cell_counts_barplot.png"。
#'
#' @return 返回一个 ggplot2 对象，并在指定路径保存一张高质量 PNG 图片。
#' @export
#'
#' @examples
#' # seurat_barplot_metadata(merged, x_var = "predicted.cell_subclass", fill_var = "Condition")
#' # seurat_barplot_metadata(merged, x_var = "Phase", fill_var = "Stage", save_path = "phase_by_stage.png")
seurat_barplot_metadata <- function(seurat_obj, x_var, fill_var, save_path = "cell_counts_barplot.png") {
  
  # 1. 检查输入的列名是否存在于元数据中
  if (!x_var %in% colnames(seurat_obj@meta.data)) {
    stop(paste("错误: 列名", x_var, "不存在于 seurat_obj@meta.data 中！"))
  }
  if (!fill_var %in% colnames(seurat_obj@meta.data)) {
    stop(paste("错误: 列名", fill_var, "不存在于 seurat_obj@meta.data 中！"))
  }

  # 2. 提取元数据并过滤掉 NA 值 (显式指定 dplyr::filter)
  meta_data <- seurat_obj@meta.data %>%
    dplyr::filter(!is.na(.data[[x_var]]), !is.na(.data[[fill_var]]))
  
  # 3. 统计每种组合的细胞数目 (全部显式指定 dplyr::)
  count_data <- meta_data %>%
    dplyr::group_by(.data[[x_var]], .data[[fill_var]]) %>%
    dplyr::tally(name = "Cell_Count") %>%
    dplyr::ungroup()

  # 计算总共有多少个不重复的 X 轴类别，用来动态调整图片宽度
  num_x_groups <- length(unique(count_data[[x_var]]))
  
  # 4. 动态计算图表的长宽 (单位：英寸)
  dynamic_width <- max(6, 4 + (num_x_groups * 1.2)) # 确保最少有 6 英寸宽
  dynamic_height <- 7
  
  # 5. 使用 ggplot2 绘图
  p <- ggplot(count_data, aes(x = .data[[x_var]], y = Cell_Count, fill = .data[[fill_var]])) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    # 在柱子上方标出细胞数目
    geom_text(
      aes(label = Cell_Count),
      position = position_dodge(0.8),
      vjust = -0.5,       
      size = 3.5,         
      fontface = "bold"
    ) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = paste("Cell Counts by", x_var, "and", fill_var),
      x = x_var,
      y = "Number of Cells",
      fill = fill_var
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"),
      panel.grid.major.x = element_blank(), 
      legend.position = "top"
    )
  
  # 6. 自动保存图片
  message(paste0("Saving plot with dynamic width: ", round(dynamic_width, 2), " inches"))
  ggsave(
    filename = save_path, 
    plot = p, 
    width = dynamic_width, 
    height = dynamic_height, 
    dpi = 300
  )
  
  return(p)
}
