library(Seurat)
library(ggplot2)

#' Generate Dynamic Validation DotPlots for Known Marker Genes
#'
#' This function takes a Seurat object, sets the active identity class, 
#' and generates two high-quality PDF DotPlots (scaled by radius and size). 
#' It automatically calculates dynamic plot dimensions based on the number of genes 
#' and clusters, allows manual fine-tuning scale factors, and accepts additional 
#' arguments for the underlying Seurat::DotPlot function.
#'
#' @param seurat_obj A standard Seurat object.
#' @param features A character vector of marker genes to visualize.
#' @param ident_var Character string specifying the metadata column to set as the active identity (e.g., "predicted.cell_subclass").
#' @param output_dir Character string specifying the directory to save the PDF files. Defaults to the current directory ("./").
#' @param file_prefix Character string used as the prefix for output file names. Defaults to "QC6.markers_provided".
#' @param width_scale Numeric factor to manually scale the dynamically calculated width. Defaults to 1.2.
#' @param height_scale Numeric factor to manually scale the dynamically calculated height. Defaults to 1.0.
#' @param ... Additional arguments passed directly to \code{\link[Seurat]{DotPlot}} (e.g., \code{assay = "RNA"}, \code{cols = c("blue", "yellow")}, etc.).
#'
#' @return Returns a list containing both ggplot2 DotPlot objects (\code{radius} and \code{size}).
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' marker_genes <- c("SOX2", "NES", "GAD1", "GAD2", "GFAP")
#' seurat_dotplot_markers(merged, features = marker_genes, ident_var = "predicted.cell_subclass")
#'
#' # Advanced usage with manual size fine-tuning and DotPlot parameter pass-through
#' seurat_dotplot_markers(
#'   seurat_obj = merged, 
#'   features = marker_genes, 
#'   ident_var = "predicted.cell_subclass",
#'   width_scale = 1.5,       # Increase width manually if gene names overlap
#'   height_scale = 0.8,      # Shrink height manually to compress layout
#'   assay = "SCT",           # Passed to Seurat::DotPlot
#'   col.min = -1,            # Passed to Seurat::DotPlot
#'   col.max = 2.5            # Passed to Seurat::DotPlot
#' )
#' }
seurat_dotplot_markers <- function(seurat_obj, features, ident_var, 
                                        output_dir = "./", 
                                        file_prefix = "QC.markers_provided",
                                        width_scale = 1,
                                        height_scale = 1,
                                        ...) {
  
  # 1. Sanity Check & Deduplication
  if (!ident_var %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Error: Identity variable '", ident_var, "' not found in seurat_obj@meta.data!"))
  }
  
  cleaned_features <- unique(features)
  # Filter out features that do not exist in the Seurat object's expression matrix
  existing_features <- cleaned_features[cleaned_features %in% rownames(seurat_obj)]
  
  if (length(existing_features) == 0) {
    stop("Error: None of the provided features were found in the Seurat object data matrix!")
  } else if (length(existing_features) < length(cleaned_features)) {
    missing_count <- length(cleaned_features) - length(existing_features)
    warning(paste(missing_count, "features were skipped because they are not present in the dataset."))
  }

  # 2. Set active identity
  Seurat::Idents(seurat_obj) <- ident_var
  
  # 3. Calculate Dynamic Layout Dimensions (Inches) and Apply Scale Factors
  base_width <- 4
  width_per_gene <- 0.2
  raw_width <- base_width + length(existing_features) * width_per_gene
  dynamic_width <- raw_width * width_scale
  
  base_height <- 2
  height_per_cluster <- 0.2
  num_clusters <- length(unique(seurat_obj[[ident_var, drop = TRUE]]))
  raw_height <- base_height + num_clusters * height_per_cluster
  dynamic_height <- raw_height * height_scale

  # 4. Define standard theme modification
  shared_theme <- ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
  
  # List to store plots for returning
  plot_list <- list()

  # 5. Generate and Save DotPlot (Scaled by Radius)
  message("Generating Marker DotPlot (scaled by radius)...")
  p_radius <- Seurat::DotPlot(seurat_obj, features = existing_features, scale.by = "radius", ...) + 
    shared_theme
  
  file_radius <- file.path(output_dir, paste0(file_prefix, ".DotPlot.radius.pdf"))
  ggplot2::ggsave(file_radius, plot = p_radius, width = dynamic_width, height = dynamic_height)
  plot_list$radius <- p_radius

  # 6. Generate and Save DotPlot (Scaled by Size)
  message("Generating Marker DotPlot (scaled by size)...")
  p_size <- Seurat::DotPlot(seurat_obj, features = existing_features, scale.by = "size", ...) + 
    shared_theme
  
  file_size <- file.path(output_dir, paste0(file_prefix, ".DotPlot.size.pdf"))
  ggplot2::ggsave(file_size, plot = p_size, width = dynamic_width, height = dynamic_height)
  plot_list$size <- p_size

  message(paste0("Success! Dynamic plots saved to '", output_dir, "'."))
  message(paste0("Final Dimensions -> Width: ", round(dynamic_width, 2), " in (scale: ", width_scale, "), Height: ", round(dynamic_height, 2), " in (scale: ", height_scale, ")."))
  
  return(plot_list)
}

