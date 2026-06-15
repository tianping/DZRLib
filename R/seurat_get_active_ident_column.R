#' Find Which Metadata Column(s) Correspond to Current Active Identity
#'
#' This function compares the current active identity of a Seurat object with all
#' metadata columns to identify which column(s) have identical values to the
#' current active identity.
#'
#' @param seurat_obj A Seurat object
#'
#' @return A character vector of column names that match the current active identity.
#'         Returns an empty character vector if no matches are found (which can happen
#'         if active identity was set manually to a custom factor not present in metadata).
#' @export
#'
#' @examples
#' \dontrun{
#' # Find current identity source(s)
#' matching_cols <- seurat_get_active_ident_column(seurat_obj)
#' if (length(matching_cols) > 0) {
#'   print(paste("Active identity matches column(s):", paste(matching_cols, collapse = ", ")))
#' } else {
#'   print("Active identity was set manually (not from metadata)")
#' }
#' }
seurat_get_active_ident_column <- function(seurat_obj) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }

  # Get current active identity values
  current_idents <- Seurat::Idents(seurat_obj)
  current_idents_char <- as.character(current_idents)

  # Get all metadata columns
  metadata <- seurat_obj[[]]
  metadata_cols <- colnames(metadata)

  # Find all matching columns
  matching_cols <- c()
  for (col in metadata_cols) {
    # Convert metadata column to character for comparison
    col_values <- as.character(metadata[[col]])

    # Check if the values (ordered by cell barcodes) match exactly
    if (identical(current_idents_char, col_values)) {
      matching_cols <- c(matching_cols, col)
    }
  }

  return(matching_cols)
}
