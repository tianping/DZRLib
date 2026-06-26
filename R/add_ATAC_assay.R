#' Add ATAC Assay to Seurat Object
#'
#' This function adds an ATAC assay to an existing RNA Seurat object, aligning cell barcodes
#' between the RNA and ATAC datasets. It identifies common cells between the two assays,
#' subsets both objects to these common cells, and adds the ATAC assay to the Seurat object.
#'
#' @param seurat_obj A Seurat object containing RNA data (must have an "RNA" assay).
#' @param atac_assay A Seurat object or assay matrix containing ATAC data.
#' @param assay_name Character string specifying the name for the ATAC assay in the Seurat object.
#'                   Default is "ATAC".
#'
#' @return Returns the updated Seurat object with the ATAC assay added.
#' @export
#'
#' @examples
#' # Basic usage
#' # Assuming 'rna_seurat' is your RNA Seurat object and 'atac_data' is your ATAC assay
#' updated_seurat <- add_ATAC_assay(rna_seurat, atac_data)
#'
#' # With custom assay name
#' updated_seurat <- add_ATAC_assay(rna_seurat, atac_data, assay_name = "ATAC_peaks")
#'
#' @seealso \code{\link{DefaultAssay}} for managing default assays in Seurat objects
#'
#' @details
#' The function performs the following steps:
#' 1. Identifies common cells between RNA and ATAC datasets
#' 2. Subsets both objects to these common cells
#' 3. Verifies cell order alignment after subsetting
#' 4. Adds the ATAC assay to the Seurat object
#' 5. Sets the default assay back to "RNA"
#'
#' @note
#' - The function requires that both RNA and ATAC objects have the same cell barcodes
#' - Cell order must be consistent between the two datasets
#' - The RNA assay is set as the default assay after adding the ATAC assay
#'
#' @md
add_ATAC_assay <- function(seurat_obj, atac_assay, assay_name = "ATAC") {

  # Find common cells between RNA and ATAC
  common_cells <- intersect(colnames(seurat_obj[["RNA"]]), colnames(atac_assay))
  message(length(common_cells), " cells in common between RNA and ATAC")

  # Subset both the Seurat object and the ATAC assay to the common cells
  seurat_obj <- subset(seurat_obj, cells = common_cells)
  atac_assay <- subset(atac_assay, cells = common_cells)

  # Check if the cell order is aligned
  if (!all(colnames(seurat_obj[["RNA"]]) == colnames(atac_assay))) {
    stop("Cell order mismatch after subsetting! Something went wrong.")
  }

  # Add the ATAC assay to the Seurat object
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj[[assay_name]] <- atac_assay

  message("ATAC assay added successfully as '", assay_name, "'")

  # Return the updated Seurat object
  return(seurat_obj)
}
