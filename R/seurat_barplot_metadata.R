library(Seurat)
library(ggplot2)
library(dplyr)

#' Plot Cell Count Barplots for Metadata Variables in a Seurat Object
#'
#' This function sets the active identity of a Seurat object to a specified metadata
#' variable and generates a barplot showing cell counts. Optionally, it can split
#' each bar by an additional metadata variable with dodged bars.
#'
#' @param seurat_obj A standard Seurat object.
#' @param x_var Character string specifying the metadata column name to set as active
#'        identity and use as the X-axis (e.g., "predicted.cell_subclass").
#' @param fill_var Character string specifying the metadata column name to use for
#'        bar splitting and fill color (e.g., "Condition"). Defaults to NULL, which
#'        produces a simple barplot without splitting.
#' @param save_path Character string specifying the output path and file name for
#'        the plot. Defaults to "cell_counts_barplot.png".
#' @param width Numeric value specifying the width of the saved plot in inches.
#'        Defaults to 12.
#' @param height Numeric value specifying the height of the saved plot in inches.
#'        Defaults to 8.
#' @param dpi Integer value specifying the resolution of the saved plot. Defaults
#'        to 300.
#' @param angle Numeric value specifying the angle of x-axis text labels. Defaults
#'        to 45.
#' @param color_palette Optional character vector of colors to use for fill_var
#'        categories. If NULL, default ggplot2 colors are used.
#'
#' @return Returns a ggplot2 object. The original Seurat object's active identity
#'         is restored after plotting (if it came from a metadata column).
#' @export
#'
#' @examples
#' \dontrun{
#' # Simple barplot with only x_var
#' seurat_barplot_metadata(seurat_obj, x_var = "predicted.cell_subclass")
#'
#' # Barplot with fill_var for dodged bars
#' seurat_barplot_metadata(seurat_obj, x_var = "predicted.cell_subclass",
#'                         fill_var = "Condition")
#'
#' # Custom save path and appearance
#' seurat_barplot_metadata(seurat_obj, x_var = "Phase", fill_var = "Stage",
#'                         save_path = "phase_by_stage.png", width = 10, height = 6)
#' }
seurat_barplot_metadata <- function(seurat_obj, x_var, fill_var = NULL, save_path = "cell_counts_barplot.png", width = 12, height = 8, dpi = 300, angle = 45, color_palette = NULL) {

  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
  }

  # Validate seurat_obj
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }

  # Validate x_var
  if (!x_var %in% colnames(seurat_obj[[]])) {
    stop(paste("x_var '", x_var, "' not found in metadata columns. Available columns: ", paste(colnames(seurat_obj[[]]), collapse = ", ")))
  }

  # Store the original identity column names (if they come from metadata)
  original_ident_cols <- seurat_get_active_ident_column(seurat_obj)

  # Set active identity to x_var
  Seurat::Idents(seurat_obj) <- x_var

  # Extract metadata
  metadata <- seurat_obj[[]]

  # Create the plot based on whether fill_var is provided
  if (is.null(fill_var)) {
    # Simple barplot - count cells per identity class
    count_data <- as.data.frame(table(Seurat::Idents(seurat_obj)))
    colnames(count_data) <- c(x_var, "Count")

    # Create barplot
    p <- ggplot2::ggplot(count_data, ggplot2::aes(x = .data[[x_var]], y = Count)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue", color = "black", width = 0.7) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste("Cell Counts by", x_var),
                    x = x_var,
                    y = "Number of Cells") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle, hjust = 1),
                     plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

    # Add count labels on top of bars
    p <- p + ggplot2::geom_text(ggplot2::aes(label = Count),
                                vjust = -0.5,
                                size = 3.5)

  } else {
    # Validate fill_var
    if (!fill_var %in% colnames(metadata)) {
      # Restore original identity before error (if possible)
      if (length(original_ident_cols) > 0) {
        # If multiple columns match, use the first one
        Seurat::Idents(seurat_obj) <- original_ident_cols[1]
      }
      stop(paste("fill_var '", fill_var, "' not found in metadata columns. Available columns: ", paste(colnames(metadata), collapse = ", ")))
    }

    # Create contingency table for x_var and fill_var
    count_table <- table(metadata[[x_var]], metadata[[fill_var]])

    # Convert to data frame for ggplot2
    plot_data <- as.data.frame(count_table)
    colnames(plot_data) <- c(x_var, fill_var, "Count")

    # Create dodged barplot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_var]],
                                                  y = Count,
                                                  fill = .data[[fill_var]])) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.9),
                        color = "black") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste("Cell Counts:", x_var, "by", fill_var),
                    x = x_var,
                    y = "Number of Cells",
                    fill = fill_var) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle, hjust = 1),
                     plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                     legend.position = "right")

    # Apply custom color palette if provided
    if (!is.null(color_palette)) {
      p <- p + ggplot2::scale_fill_manual(values = color_palette)
    }

    # Add count labels on bars
    p <- p + ggplot2::geom_text(ggplot2::aes(label = Count),
                                position = ggplot2::position_dodge(width = 0.9),
                                vjust = -0.5,
                                size = 3)
  }

  # Save the plot
  ggplot2::ggsave(save_path, plot = p, width = width, height = height, dpi = dpi)

  # Restore original identity using the column name(s) (if they exist)
  if (length(original_ident_cols) > 0) {
    # If multiple columns match, use the first one (they are identical anyway)
    Seurat::Idents(seurat_obj) <- original_ident_cols[1]
    if (length(original_ident_cols) == 1) {
      message(paste("Restored original identity from column:", original_ident_cols[1]))
    } else {
      message(paste("Restored original identity from column:", original_ident_cols[1]))
      message(paste("Note: Multiple columns matched the original identity:",
                    paste(original_ident_cols, collapse = ", ")))
    }
  } else {
    # If original identity was set manually, we cannot restore it
    warning("Original active identity was set manually (not from a metadata column) and cannot be automatically restored. The Seurat object now uses '", x_var, "' as identity.")
  }

  message(paste("Plot saved to:", save_path))
  if (!is.null(fill_var)) {
    message(paste("Barplot created with", x_var, "on x-axis and", fill_var, "as fill variable (dodged bars)"))
  } else {
    message(paste("Simple barplot created with", x_var, "on x-axis"))
  }

  # Return the plot object
  return(p)
}
