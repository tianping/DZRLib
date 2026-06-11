library(Seurat)
library(ggplot2)
library(dplyr)

#' Plot Dynamic Cell Count Barplots for Two Combined Metadata Variables in a Seurat Object
#'
#' @param seurat_obj A standard Seurat object.
#' @param x_var Character string specifying the metadata column name to use as the X-axis (e.g., "predicted.cell_subclass").
#' @param fill_var Character string specifying the metadata column name to use for bar splitting and fill color (e.g., "Condition").
#' @param save_path Character string specifying the output path and file name for the plot. Defaults to "cell_counts_barplot.png".
#'
#' @return Returns a ggplot2 object and saves a high-resolution PNG image to the specified path.
#' @export
#'
#' @examples
#' # seurat_barplot_metadata(merged, x_var = "predicted.cell_subclass", fill_var = "Condition")
#' # seurat_barplot_metadata(merged, x_var = "Phase", fill_var = "Stage", save_path = "phase_by_stage.png")
seurat_barplot_metadata <- function(seurat_obj, x_var, fill_var, save_path = "cell_counts_barplot.png") {

  # 1. Verify if the specified variables exist in the metadata
  if (!x_var %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Error: Column name '", x_var, "' not found in seurat_obj@meta.data!"))
  }
  if (!fill_var %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Error: Column name '", fill_var, "' not found in seurat_obj@meta.data!"))
  }

  # 2. Extract metadata and filter out NA values
  meta_data <- seurat_obj@meta.data |>
    dplyr::filter(!is.na(.data[[x_var]]), !is.na(.data[[fill_var]]))

  # 3. Calculate cell count for each group combination
  count_data <- meta_data |>
    dplyr::group_by(.data[[x_var]], .data[[fill_var]]) |>
    dplyr::tally(name = "Cell_Count") |>
    dplyr::ungroup()

  # Count unique categories along X-axis to dynamically scale the plot width
  num_x_groups <- length(unique(count_data[[x_var]]))

  # 4. Dynamically compute dimensions (in inches)
  dynamic_width <- max(6, 4 + (num_x_groups * 1.2)) # Ensure a minimum width of 6 inches
  dynamic_height <- 7

  # 5. Generate plot using ggplot2
  p <- ggplot(count_data, aes(x = .data[[x_var]], y = Cell_Count, fill = .data[[fill_var]])) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    # Add text labels displaying cell counts above bars
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

  # 6. Automatically save the plot
  message(paste0("Saving plot with dynamic width: ", round(dynamic_width, 2), " inches"))
  ggplot2::ggsave(
    filename = save_path,
    plot = p,
    width = dynamic_width,
    height = dynamic_height,
    dpi = 300
  )

  return(p)
}
