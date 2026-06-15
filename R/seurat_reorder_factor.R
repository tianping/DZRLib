#' Reorder Factor Levels in Seurat Metadata
#'
#' This function reorders the factor levels of a specified metadata column in a 
#' Seurat object. It supports various ordering methods and can automatically sync 
#' the active identity if the modified column is currently used as identity.
#'
#' @param seurat_obj A Seurat object.
#' @param column Character string specifying the metadata column name to reorder.
#' @param order_by Character string specifying the ordering method. Options:
#'        "custom" (user-defined order), "alphabetical", "numeric", "cell_count", 
#'        or "original". Default is "custom".
#' @param custom_order Character vector specifying the desired level order when 
#'        order_by = "custom". Must include all levels exactly once.
#' @param decreasing Logical value indicating whether to sort in decreasing order. 
#'        Applies to "alphabetical", "numeric", and "cell_count" methods. 
#'        Default is FALSE.
#' @param na.last Character string or logical controlling the placement of NAs.
#'        Options: "keep" (keep NA as a level), "remove" (remove NA values), 
#'        "first" (place NA first), "last" (place NA last). Default is "keep".
#' @param drop_unused_levels Logical value indicating whether to remove factor 
#'        levels that have zero cells. Default is FALSE.
#' @param sync_ident Logical value indicating whether to automatically sync the 
#'        active identity if the modified column is currently used as identity. 
#'        Default is TRUE.
#'
#' @return A Seurat object with the modified metadata column. The object is also 
#'         modified in place.
#' @export
#'
#' @examples
#' \dontrun{
#' # Custom order
#' seurat_obj <- seurat_reorder_factor(seurat_obj, "cell_type", 
#'                                      order_by = "custom",
#'                                      custom_order = c("T cell", "B cell", "Monocyte"))
#' 
#' # Alphabetical order descending
#' seurat_obj <- seurat_reorder_factor(seurat_obj, "condition", 
#'                                      order_by = "alphabetical", 
#'                                      decreasing = TRUE)
#' 
#' # Order by cell count
#' seurat_obj <- seurat_reorder_factor(seurat_obj, "cluster", 
#'                                      order_by = "cell_count", 
#'                                      decreasing = TRUE)
#' 
#' # Remove unused levels and sync identity
#' seurat_obj <- seurat_reorder_factor(seurat_obj, "predicted.cell_subclass", 
#'                                      order_by = "alphabetical",
#'                                      drop_unused_levels = TRUE)
#' }
seurat_reorder_factor <- function(seurat_obj, 
                                  column, 
                                  order_by = "custom", 
                                  custom_order = NULL,
                                  decreasing = FALSE,
                                  na.last = "keep",
                                  drop_unused_levels = FALSE,
                                  sync_ident = TRUE) {
  
  # Check if Seurat is installed
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
  }
  
  # Validate seurat_obj
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  
  # Validate column exists
  if (!column %in% colnames(seurat_obj[[]])) {
    stop(paste("Column '", column, "' not found in metadata. Available columns: ", 
               paste(colnames(seurat_obj[[]]), collapse = ", ")))
  }
  
  # Get the column values
  col_values <- seurat_obj[[column]][[1]]
  
  # Check if column is factor
  if (!is.factor(col_values)) {
    stop(paste("Column '", column, "' is not a factor. Current class: ", 
               class(col_values)[1], 
               "\nPlease convert it to factor first using:\n",
               "seurat_obj$", column, " <- factor(seurat_obj$", column, ")", sep = ""))
  }
  
  # Get current levels and cell counts
  current_levels <- levels(col_values)
  
  # Calculate cell counts for each level
  cell_counts <- table(col_values)
  
  # Handle NA values based on na.last parameter
  if (na.last == "remove") {
    # Remove NA values from the column
    na_present <- any(is.na(col_values))
    if (na_present) {
      col_values <- col_values[!is.na(col_values)]
      seurat_obj[[column]] <- col_values
      message("Removed ", sum(is.na(col_values)), " NA values from column '", column, "'")
    }
  }
  
  # Determine new level order based on order_by method
  unique_levels <- unique(col_values[!is.na(col_values)])
  unique_levels <- as.character(unique_levels)
  
  if (order_by == "custom") {
    # Validate custom_order
    if (is.null(custom_order)) {
      stop("When order_by = 'custom', custom_order must be provided.")
    }
    
    missing_levels <- setdiff(unique_levels, custom_order)
    extra_levels <- setdiff(custom_order, unique_levels)
    
    if (length(missing_levels) > 0) {
      stop(paste("custom_order is missing the following levels:", 
                 paste(missing_levels, collapse = ", ")))
    }
    
    if (length(extra_levels) > 0) {
      stop(paste("custom_order contains levels not in the data:", 
                 paste(extra_levels, collapse = ", ")))
    }
    
    new_levels <- custom_order
    
  } else if (order_by == "alphabetical") {
    new_levels <- sort(unique_levels, decreasing = decreasing)
    
  } else if (order_by == "numeric") {
    # Try to convert to numeric
    numeric_levels <- suppressWarnings(as.numeric(unique_levels))
    
    if (all(!is.na(numeric_levels))) {
      new_levels <- as.character(sort(numeric_levels, decreasing = decreasing))
    } else {
      warning("Cannot convert levels to numeric. Falling back to alphabetical order.")
      new_levels <- sort(unique_levels, decreasing = decreasing)
    }
    
  } else if (order_by == "cell_count") {
    # Get cell counts for each level
    count_df <- data.frame(level = unique_levels, 
                           count = as.numeric(table(col_values)[unique_levels]))
    
    # Sort by count, then alphabetically for ties
    count_df <- count_df[order(count_df$count, count_df$level, decreasing = decreasing), ]
    new_levels <- as.character(count_df$level)
    
  } else if (order_by == "original") {
    # Keep original order (use current levels if they exist, otherwise unique appearance order)
    if (length(current_levels) > 0) {
      new_levels <- current_levels
    } else {
      new_levels <- unique(as.character(col_values))
    }
    
  } else {
    stop(paste("order_by must be one of: 'custom', 'alphabetical', 'numeric', 'cell_count', 'original'"))
  }
  
  # Handle NA placement
  if (na.last == "first" && any(is.na(col_values))) {
    new_levels <- c(NA, new_levels)
  } else if (na.last == "last" && any(is.na(col_values))) {
    new_levels <- c(new_levels, NA)
  } else if (na.last == "keep" && any(is.na(col_values))) {
    # Keep NA in its current position relative to new_levels
    # Since NA is not typically in levels, we add it at the end
    if (!any(is.na(new_levels))) {
      new_levels <- c(new_levels, NA)
    }
  }
  
  # Remove unused levels if requested
  if (drop_unused_levels) {
    used_levels <- unique(as.character(col_values[!is.na(col_values)]))
    new_levels <- new_levels[new_levels %in% used_levels]
  }
  
  # Apply the new factor levels
  seurat_obj[[column]] <- factor(seurat_obj[[column]], levels = new_levels)
  
  # Output information about the change
  message("\n=== Factor Reordering Summary ===")
  message("Column: ", column)
  message("Ordering method: ", order_by)
  if (order_by == "custom") {
    message("Custom order provided")
  }
  message("\nOriginal level order:")
  message(paste("  ", paste(current_levels, collapse = " -> ")))
  message("\nNew level order:")
  message(paste("  ", paste(new_levels, collapse = " -> ")))
  
  if (drop_unused_levels) {
    removed_levels <- setdiff(current_levels, new_levels)
    if (length(removed_levels) > 0) {
      message("\nRemoved unused levels:")
      message(paste("  ", paste(removed_levels, collapse = ", ")))
    }
  }
  
  message("\nCell counts per level:")
  new_counts <- table(seurat_obj[[column]])
  for (lev in new_levels) {
    if (is.na(lev)) {
      message(paste("  NA:", new_counts[is.na(names(new_counts))]))
    } else {
      message(paste("  ", lev, ":", new_counts[lev]))
    }
  }
  
  # Sync active identity if needed
  if (sync_ident) {
    matching_cols <- seurat_get_active_ident_column(seurat_obj)
    
    if (column %in% matching_cols) {
      Seurat::Idents(seurat_obj) <- column
      message("\n✓ Active identity synced to updated column '", column, "'")
    } else {
      message("\n※ Column '", column, "' is not currently used as active identity. No syncing needed.")
    }
  }
  
  message("\n=== End of Summary ===\n")
  
  # Return the modified Seurat object
  return(seurat_obj)
}

