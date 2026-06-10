library(Seurat)
library(dplyr)

#' Save Seurat Object Metadata to TSV File
#'
#' This function saves the metadata from a Seurat object to a TSV file. By default,
#' it saves all columns from the metadata, but users can specify which columns to include.
#' The default filename is based on the Seurat object name with a ".metadata.tsv" extension.
#'
#' @param seurat_obj A standard Seurat object.
#' @param columns Optional character vector of column names to save. If NULL (default),
#'                all columns from the metadata will be saved.
#' @param filename Character string specifying the output file path. Defaults to
#'                 "[seurat-object-name].metadata.tsv".
#'
#' @return Returns the path to the saved file, invisibly.
#' @export
#'
#' @examples
#' # Save all metadata with default filename
#' seurat_save_metadata(merged)
#'
#' # Save specific columns with custom filename
#' seurat_save_metadata(merged, columns = c("nCount_RNA", "nFeature_RNA", "Condition"),
#'                      filename = "custom_metadata.tsv")
seurat_save_metadata <- function(seurat_obj, columns = NULL, filename = NULL) {
  # Get the metadata from the Seurat object
  meta_data <- seurat_obj@meta.data

  # Filter columns if specified
  if (!is.null(columns)) {
    # Check if all requested columns exist
    missing_cols <- setdiff(columns, colnames(meta_data))
    if (length(missing_cols) > 0) {
      stop(paste("Error: The following columns were not found in the metadata:",
                 paste(missing_cols, collapse = ", ")))
    }
    meta_data <- meta_data[, columns, drop = FALSE]
  }

  # Generate default filename if not provided
  if (is.null(filename)) {
    # Extract object name (handle cases where object might be unnamed)
    obj_name <- deparse(substitute(seurat_obj))
    # Clean the name to remove any non-alphanumeric characters except underscores
    clean_name <- gsub("[^a-zA-Z0-9_]", "_", obj_name)
    filename <- paste0(clean_name, ".metadata.tsv")
  }

  # Save to TSV file
  write.table(meta_data, file = filename, sep = "\t", row.names = TRUE,
              quote = FALSE, col.names = TRUE)

  message(paste("Metadata saved to:", filename))
  invisible(filename)
}
