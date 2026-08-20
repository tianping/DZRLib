#' Automated Pseudobulk Differential Expression Analysis for Seurat
#'
#' Performs pairwise differential expression (DE) analysis on pseudobulk objects
#' or raw single-cell Seurat objects. Supports both stratified/nested comparisons
#' (e.g., comparing conditions within each cell type) and global comparisons.
#' Automatically exports aggregated count tables if AggregateExpression is run.
#'
#' @param object A Seurat object (either a raw single-cell object or a pre-aggregated pseudobulk object).
#' @param compare_col Character. Metadata column name used for pairwise comparisons (e.g., "origin_site" or "stim").
#' @param split_col Character, optional. Metadata column name used for stratifying/nesting the analysis (e.g., "CellType"). If \code{NULL}, performs global comparison across \code{compare_col}.
#' @param group_by Character vector, optional. Metadata columns passed to \code{AggregateExpression} if \code{object} is a raw single-cell object (e.g., \code{c("orig.ident", "origin_site", "CellType")}).
#' @param target_splits Character vector, optional. Specific subsets of \code{split_col} to analyze (e.g., \code{c("BEC1", "FB1")}).
#' @param comparisons List of character vectors, optional. Specific pairwise comparisons to run, e.g., \code{list(c("Back", "Ear"))}. If \code{NULL}, all pairwise combinations of groups in \code{compare_col} will be analyzed.
#' @param agg_args List of additional arguments passed to \code{\link[Seurat]{AggregateExpression}}, e.g., \code{list(assays = "RNA", slot = "counts")}.
#' @param findmarkers_args List of additional arguments passed to \code{\link[Seurat]{FindMarkers}}, default is \code{list(test.use = "DESeq2")}.
#' @param out_dir Character, optional. Directory path to save output CSV files. If \code{NULL}, results are not written to disk.
#' @param existing_results List, optional. An existing result list to incrementally update/append new results.
#' @param min_cells Integer, default 2. Minimum sample/pseudobulk cell count required per comparison group.
#'
#' @return A named list of data.frames containing DE analysis results for each comparison, and optionally the aggregated count table under \code{"_aggregated_counts"}.
#' @export
#'
#' @examples
#' \dontrun{
#' # Scenario 1: Input raw single-cell Seurat object with stratification by CellType
#' de_results <- seurat_RunPseudobulkDE(
#'   object = my.exp,
#'   compare_col = "origin_site",
#'   split_col = "CellType",
#'   group_by = c("orig.ident", "origin_site", "CellType"),
#'   agg_args = list(assays = "RNA", slot = "counts"),
#'   findmarkers_args = list(test.use = "DESeq2"),
#'   out_dir = "./results_by_celltype"
#' )
#'
#' # Scenario 2: Input pre-aggregated pseudobulk object for global comparison
#' de_global <- seurat_RunPseudobulkDE(
#'   object = my.exp.pseudo,
#'   compare_col = "stim",
#'   split_col = NULL,
#'   findmarkers_args = list(test.use = "DESeq2"),
#'   out_dir = "./results_global"
#' )
#'
#' # Scenario 3: Subset specific cell types and comparison pairs with incremental update
#' de_results_updated <- seurat_RunPseudobulkDE(
#'   object = my.exp.pseudo,
#'   compare_col = "origin_site",
#'   split_col = "CellType",
#'   target_splits = c("BEC1", "FB1"),
#'   comparisons = list(c("Back", "Ear")),
#'   existing_results = de_results,
#'   out_dir = "./results_updated"
#' )
#'
#' # Access specific comparison result from list
#' head(de_results[["BEC1___Back_vs_Ear"]])
#' }
seurat_RunPseudobulkDE <- function(object,
                                   compare_col,
                                   split_col = NULL,
                                   group_by = NULL,
                                   target_splits = NULL,
                                   comparisons = NULL,
                                   agg_args = list(assays = "RNA", return.seurat = TRUE),
                                   findmarkers_args = list(test.use = "DESeq2"),
                                   out_dir = NULL,
                                   existing_results = list(),
                                   min_cells = 2) {
  
  # ---------------------------------------------------------------------------
  # 1. Pseudobulk aggregation & Save Aggregated Counts (if group_by is given)
  # ---------------------------------------------------------------------------
  seurat_obj <- object
  aggregated_counts <- NULL
  
  if (!is.null(group_by)) {
    message(">>> Running AggregateExpression to build Pseudobulk object...")
    
    # Extract raw matrix output for table export
    agg_params_mat <- c(list(object = seurat_obj, group.by = group_by, return.seurat = FALSE), agg_args)
    agg_params_mat <- agg_params_mat[!duplicated(names(agg_params_mat))]
    aggregated_counts <- do.call(Seurat::AggregateExpression, agg_params_mat)
    
    # Create aggregated Seurat object
    agg_params_obj <- c(list(object = seurat_obj, group.by = group_by, return.seurat = TRUE), agg_args)
    agg_params_obj <- agg_params_obj[!duplicated(names(agg_params_obj))]
    seurat_obj <- do.call(Seurat::AggregateExpression, agg_params_obj)
    
    # Write aggregated counts to disk
    if (!is.null(out_dir)) {
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      
      assay_name <- names(aggregated_counts)[1]
      cnt_mat <- aggregated_counts[[assay_name]]
      cnt_df <- as.data.frame(as.matrix(cnt_mat))
      cnt_df <- tibble::rownames_to_column(cnt_df, var = "gene")
      
      utils::write.csv(cnt_df, file = file.path(out_dir, "aggregated_pseudobulk_counts.csv"), row.names = FALSE)
      message(">>> Saved aggregated count matrix to: ", file.path(out_dir, "aggregated_pseudobulk_counts.csv"))
    }
  }
  
  if (!is.null(out_dir) && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  results_list <- existing_results
  if (!is.null(aggregated_counts)) {
    results_list[["_aggregated_counts"]] <- aggregated_counts
  }
  
  all_combined_df <- list()
  
  # ---------------------------------------------------------------------------
  # 2. Build identities and set up comparison mode
  # ---------------------------------------------------------------------------
  if (!is.null(split_col)) {
    # Mode A: Stratified comparison within each split_col level
    seurat_obj$de_group <- paste(seurat_obj[[compare_col]][, 1], 
                                 seurat_obj[[split_col]][, 1], 
                                 sep = "_")
    Seurat::Idents(seurat_obj) <- "de_group"
    
    all_splits <- unique(as.character(seurat_obj[[split_col]][, 1]))
    splits_to_run <- if (!is.null(target_splits)) intersect(target_splits, all_splits) else all_splits
  } else {
    # Mode B: Global comparison across compare_col
    Seurat::Idents(seurat_obj) <- compare_col
    splits_to_run <- "GLOBAL"
  }
  
  # ---------------------------------------------------------------------------
  # 3. Iterate through subsets and run pairwise FindMarkers
  # ---------------------------------------------------------------------------
  for (sp in splits_to_run) {
    
    if (sp == "GLOBAL") {
      sub_meta <- seurat_obj@meta.data
      message("\n>>> Running Global Pairwise DE Analysis [Column: ", compare_col, "]")
    } else {
      sub_meta <- seurat_obj@meta.data[seurat_obj[[split_col]][, 1] == sp, ]
      message("\n>>> Processing Subset [", split_col, " = ", sp, "]")
    }
    
    # Check sample counts per group
    counts_table <- table(sub_meta[[compare_col]])
    valid_groups <- names(counts_table[counts_table >= min_cells])
    
    if (length(valid_groups) < 2) {
      warning("Skipping ", sp, ": Less than 2 groups meet the minimum sample threshold (>= ", min_cells, ").")
      next
    }
    
    # Define pairwise comparison pairs
    curr_comps <- list()
    if (!is.null(comparisons)) {
      for (comp in comparisons) {
        if (all(comp %in% valid_groups)) {
          curr_comps[[paste(comp[1], comp[2], sep = "_vs_")]] <- comp
        }
      }
    } else {
      combos <- utils::combn(valid_groups, 2, simplify = FALSE)
      for (cb in combos) {
        curr_comps[[paste(cb[1], cb[2], sep = "_vs_")]] <- cb
      }
    }
    
    # Execute DE analysis
    for (comp_name in names(curr_comps)) {
      pair <- curr_comps[[comp_name]]
      
      if (sp == "GLOBAL") {
        id1 <- pair[1]
        id2 <- pair[2]
        list_key <- paste(compare_col, comp_name, sep = "___")
      } else {
        id1 <- paste(pair[1], sp, sep = "_")
        id2 <- paste(pair[2], sp, sep = "_")
        list_key <- paste(sp, comp_name, sep = "___")
      }
      
      message("  -> Running DE: ", id1, " vs ", id2)
      
      fm_params <- c(list(object = seurat_obj, ident.1 = id1, ident.2 = id2, verbose = FALSE), findmarkers_args)
      fm_params <- fm_params[!duplicated(names(fm_params))]
      
      tryCatch({
        de_res <- do.call(Seurat::FindMarkers, fm_params)
        
        # Add metadata columns
        de_res <- de_res %>%
          tibble::rownames_to_column(var = "gene") %>%
          dplyr::mutate(Comparison = comp_name)
        
        if (sp != "GLOBAL") {
          de_res[[split_col]] <- sp
        }
        
        # Save to memory list
        results_list[[list_key]] <- de_res
        all_combined_df[[list_key]] <- de_res
        
        # Output individual CSV file
        if (!is.null(out_dir)) {
          prefix <- if (sp == "GLOBAL") "DE_GLOBAL" else paste0("DE_", sp)
          file_path <- file.path(out_dir, paste0(prefix, "_", comp_name, ".csv"))
          utils::write.csv(de_res, file = file_path, row.names = FALSE)
        }
        
      }, error = function(e) {
        message("   [Error] Execution failed for (", list_key, "): ", e$message)
      })
    }
  }
  
  # ---------------------------------------------------------------------------
  # 4. Export combined summary CSV
  # ---------------------------------------------------------------------------
  if (!is.null(out_dir) && length(all_combined_df) > 0) {
    merged_df <- dplyr::bind_rows(all_combined_df)
    utils::write.csv(merged_df, file = file.path(out_dir, "ALL_combined_pseudobulk_DE.csv"), row.names = FALSE)
    message("\n>>> All tasks completed. Integrated results written to: ", file.path(out_dir, "ALL_combined_pseudobulk_DE.csv"))
  }
  
  return(results_list)
}
