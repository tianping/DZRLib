#' Save Subsetted DimPlots Iteratively to Files
#'
#' Iterates through distinct levels or categories of a metadata column in a
#' Seurat object, generates subsetted dimension reduction plots (\code{DimPlot}),
#' and exports them to PDF files with automatically generated, descriptive file names.
#'
#' @param object A Seurat object.
#' @param split_by Character. Metadata column name used to split cells and iterate
#'   over (e.g., \code{"orig.ident"} or \code{"seurat_clusters"}).
#' @param group_by Character or \code{NULL}. Metadata column name used to color
#'   cells in \code{DimPlot}. If \code{NULL}, uses the default grouping variable.
#' @param reduction Character. Reduction technique to visualize (e.g., \code{"umap"},
#'   \code{"tsne"}, \code{"pca"}). Default is \code{"umap"}.
#' @param output_dir Character. Directory path where PDF files will be saved.
#'   Automatically created if it does not exist. Default is \code{"."}.
#' @param file_prefix Character. Prefix for generated PDF file names. Default is \code{"DimPlot"}.
#' @param width Numeric. Width of saved PDF files in inches. Default is \code{8}.
#' @param height Numeric. Height of saved PDF files in inches. Default is \code{7}.
#' @param ... Additional arguments passed directly to \code{\link[Seurat]{DimPlot}}
#'   (e.g., \code{pt.size}, \code{cols}, \code{label}, \code{label.size}, \code{raster}).
#'
#' @return Invisibly returns \code{NULL}. Plots are saved as PDF files in \code{output_dir}.
#'
#' @importFrom Seurat DimPlot
#' @importFrom ggplot2 ggtitle theme element_text ggsave
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#'
#' # Scenario 1: Split by sample, default cluster coloring
#' seurat_split_dimplots(
#'   object = my.exp,
#'   split_by = "orig.ident",
#'   reduction = "umap",
#'   output_dir = "results/umap_samples",
#'   label = TRUE,
#'   label.size = 6
#' )
#'
#' # Scenario 2: Split by cluster, colored by sample
#' seurat_split_dimplots(
#'   object = my.exp,
#'   split_by = "seurat_clusters",
#'   group_by = "orig.ident",
#'   reduction = "umap",
#'   output_dir = "results/umap_clusters"
#' )
#'
#' # Scenario 3: Split by cluster, colored by tissue origin with custom colors
#' seurat_split_dimplots(
#'   object = my.exp,
#'   split_by = "seurat_clusters",
#'   group_by = "origin_site",
#'   reduction = "umap",
#'   output_dir = "results/umap_origin",
#'   file_prefix = "TissueOrigin",
#'   cols = c("Ear" = "orange", "Back" = "purple"),
#'   pt.size = 0.8
#' )
#' }
seurat_split_dimplots <- function(object,
                                  split_by,
                                  group_by = NULL,
                                  reduction = "umap",
                                  output_dir = ".",
                                  file_prefix = "DimPlot",
                                  width = 8,
                                  height = 7,
                                  ...) {
  # 1. Input Validation
  if (!inherits(object, "Seurat")) {
    stop("Input 'object' must be a Seurat object.")
  }

  meta_data <- object@meta.data

  if (!split_by %in% colnames(meta_data)) {
    stop(sprintf("Column '%s' was not found in object metadata.", split_by))
  }

  if (!is.null(group_by) && !group_by %in% colnames(meta_data)) {
    stop(sprintf("Column '%s' was not found in object metadata.", group_by))
  }

  # 2. Ensure Output Directory Exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # 3. Obtain Unique Categories for Iteration
  split_factor <- meta_data[[split_by]]
  split_groups <- if (is.factor(split_factor)) levels(split_factor) else unique(split_factor)
  split_groups <- split_groups[!is.na(split_groups)]

  # 4. Iterative Plot Generation and Export
  for (grp in split_groups) {
    # Extract cell IDs for the current subset
    selected_cells <- rownames(meta_data)[which(meta_data[[split_by]] == grp)]

    # Skip iteration if no cells match
    if (length(selected_cells) == 0) next

    # Build DimPlot with forwarded '...' parameters
    p <- Seurat::DimPlot(
      object = object,
      cells = selected_cells,
      reduction = reduction,
      group.by = group_by,
      ...
    ) +
      ggplot2::ggtitle(as.character(grp)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    # Construct clean, informative filename
    grp_clean <- gsub("[^[:alnum:]_.-]", "_", as.character(grp))
    split_clean <- gsub("[^[:alnum:]_.-]", "_", split_by)

    if (is.null(group_by)) {
      file_name <- sprintf("%s_%s_split-%s_%s.pdf",
                           file_prefix, reduction, split_clean, grp_clean)
    } else {
      group_clean <- gsub("[^[:alnum:]_.-]", "_", group_by)
      file_name <- sprintf("%s_%s_split-%s_%s_by-%s.pdf",
                           file_prefix, reduction, split_clean, grp_clean, group_clean)
    }

    full_path <- file.path(output_dir, file_name)

    # Save plot to PDF
    ggplot2::ggsave(
      filename = full_path,
      plot = p,
      width = width,
      height = height
    )
  }

  invisible(NULL)
}
