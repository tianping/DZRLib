#' Generate and Save DimPlot for Seurat Object
#'
#' This function creates a dimensionality reduction plot (DimPlot) for a Seurat object
#' and saves it as a PDF file. It accepts all standard DimPlot parameters and allows
#' customization of the output file path, dimensions, and other plotting parameters.
#'
#' @param seurat_obj A standard Seurat object.
#' @param reduction Character string specifying the dimensionality reduction to use (e.g., "pca", "umap", "tsne").
#'                  Default is "umap".
#' @param output_dir Directory path where the plot will be saved. Default is current directory ("./").
#' @param output_file Name of the output file (without extension). Default is "DimPlot".
#' @param width Numeric value for plot width in inches. Default is 8.
#' @param height Numeric value for plot height in inches. Default is 6.
#' @param ... Additional arguments to pass to the DimPlot function (e.g., group.by, label, pt.size, etc.).
#'
#' @return Invisibly returns the path to the saved plot file.
#' @export
#'
#' @examples
#' # Basic usage
#' seurat_dimplot(my.seurat)
#'
#' # With custom reduction and output file
#' seurat_dimplot(my.seurat, reduction = "pca", output_file = "pca_plot")
#'
#' # With additional DimPlot parameters
#' seurat_dimplot(my.seurat, group.by = "cell.type", label = TRUE, pt.size = 0.5)
seurat_dimplot <- function(seurat_obj,
                           reduction = "umap",
                           output_dir = "./",
                           output_file = "DimPlot",
                           width = 8,
                           height = 6,
                           ...) {
  # Input validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }

  if (!is.numeric(width) || width <= 0) {
    stop("width must be a positive number")
  }

  if (!is.numeric(height) || height <= 0) {
    stop("height must be a positive number")
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate the DimPlot
  p <- do.call(DimPlot, c(list(seurat_obj, reduction = reduction), list(...)))

  # Save the plot as PDF
  output_path <- file.path(output_dir, paste0(output_file, ".pdf"))
  ggplot2::ggsave(output_path,
                 plot = p,
                 width = width,
                 height = height,
                 device = "pdf")

  message(paste("DimPlot saved to:", output_path))
  invisible(output_path)
}
