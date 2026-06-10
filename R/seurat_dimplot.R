#' Generate and Save DimPlot for Seurat Object
#'
#' This function creates a dimensionality reduction plot (DimPlot) for a Seurat object
#' and saves it as a file. It accepts all standard DimPlot parameters and allows
#' customization of the output file path, dimensions, and other plotting parameters.
#'
#' @param seurat_obj A standard Seurat object.
#' @param reduction Character string specifying the dimensionality reduction to use (e.g., "pca", "umap", "tsne").
#'                  Default is "umap".
#' @param ident Optional character string specifying the metadata column to set as the active identity.
#'              If provided, the function will set Idents(seurat_obj) to this value and include it in the filename.
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
#' # With identity parameter (includes in filename)
#' seurat_dimplot(my.seurat, ident = "cell.type", output_file = "dimplot_by_celltype")
#'
#' # With additional DimPlot parameters
#' seurat_dimplot(my.seurat, group.by = "cell.type", label = TRUE, pt.size = 0.5)
seurat_dimplot <- function(seurat_obj,
                           reduction = "umap",
                           ident = NULL,
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

  # Set active identity if provided
  if (!is.null(ident)) {
    if (!ident %in% colnames(seurat_obj@meta.data)) {
      stop(paste("Error: Identity variable '", ident, "' not found in seurat_obj@meta.data!"))
    }
    Seurat::Idents(seurat_obj) <- ident
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate the DimPlot
  p <- do.call(DimPlot, c(list(seurat_obj, reduction = reduction), list(...)))

  # Determine file extension from output_file
  file_parts <- strsplit(output_file, "\\.")[[1]]
  if (length(file_parts) > 1) {
    file_ext <- tolower(file_parts[length(file_parts)])
    if (file_ext %in% c("pdf", "png", "jpg", "jpeg")) {
      output_file_base <- paste(file_parts[1:(length(file_parts)-1)], collapse = ".")
    } else {
      output_file_base <- output_file
      file_ext <- "pdf"  # default to pdf if extension is not recognized
    }
  } else {
    output_file_base <- output_file
    file_ext <- "pdf"  # default extension
  }

  # Add ident to filename if provided
  if (!is.null(ident)) {
    output_file_base <- paste0(output_file_base, ".", ident)
  }

  # Save the plot
  output_path <- file.path(output_dir, paste0(output_file_base, ".", file_ext))
  ggplot2::ggsave(output_path,
                 plot = p,
                 width = width,
                 height = height,
                 device = file_ext)

  message(paste("DimPlot saved to:", output_path))
  invisible(output_path)
}
