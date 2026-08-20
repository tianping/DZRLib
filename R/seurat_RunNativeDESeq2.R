#' Automated Native DESeq2 Analysis for Single-Cell Pseudobulk Data
#'
#' Performs native DESeq2 differential expression (DE) analysis on pseudobulk objects
#' or raw single-cell Seurat objects. Bypasses Seurat's FindMarkers wrapper to allow
#' multi-factor designs (covariates), LFC shrinkage, and output of normalized count tables.
#'
#' @param object A Seurat object (either raw single-cell or pre-aggregated pseudobulk).
#' @param compare_col Character. Primary metadata column used for DE comparison (e.g., "origin_site").
#' @param split_col Character, optional. Metadata column used for subset stratification (e.g., "CellType"). If \code{NULL}, performs global comparison.
#' @param covariates Character vector, optional. Metadata columns for covariates / batch effect correction (e.g., \code{c("donor_id", "sex")}).
#' @param group_by Character vector, optional. Metadata columns passed to \code{AggregateExpression} if \code{object} is raw single-cell data.
#' @param target_splits Character vector, optional. Specific subsets of \code{split_col} to analyze.
#' @param comparisons List of character vectors, optional. Pairwise comparison pairs where first element is numerator (group1) and second is denominator/control (group2). E.g., \code{list(c("Back", "Ear"))}.
#' @param min_count_sum Numeric, default 0. Minimum total count threshold across samples for gene pre-filtering. Set to 0 for no filtering.
#' @param lfc_shrink Logical, default \code{TRUE}. Whether to apply Log2 Fold Change shrinkage using \code{DESeq2::lfcShrink}.
#' @param shrink_type Character, default "ashr". LFC shrinkage algorithm type ("ashr", "normal", or "apeglm").
#' @param agg_args List of additional arguments passed to \code{\link[Seurat]{AggregateExpression}}.
#' @param out_dir Character, optional. Directory path to save output CSV files and count tables. If \code{NULL}, no files are written to disk.
#' @param existing_results List, optional. Existing result list to incrementally update/append new DE tables.
#' @param min_cells Integer, default 2. Minimum pseudobulk sample count required per comparison group.
#'
#' @return A named list containing DE result data.frames (with normalized counts appended) and optionally the aggregated count table.
#' @export
#'
#' @examples
#' \dontrun{
#' # Scenario 1: Native DESeq2 with covariates and CellType stratification
#' de_results <- seurat_RunNativeDESeq2(
#'   object = my.exp,
#'   compare_col = "origin_site",
#'   split_col = "CellType",
#'   covariates = c("donor_id"),
#'   group_by = c("orig.ident", "origin_site", "CellType", "donor_id"),
#'   min_count_sum = 10,
#'   lfc_shrink = TRUE,
#'   out_dir = "./deseq2_results"
#' )
#'
#' # Scenario 2: Direct comparison on pre-aggregated object without gene filtering
#' de_results_global <- seurat_RunNativeDESeq2(
#'   object = my.exp.pseudo,
#'   compare_col = "stim",
#'   split_col = NULL,
#'   min_count_sum = 0,
#'   out_dir = "./deseq2_global"
#' )
#' }
seurat_RunNativeDESeq2 <- function(object,
                                   compare_col,
                                   split_col = NULL,
                                   covariates = NULL,
                                   group_by = NULL,
                                   target_splits = NULL,
                                   comparisons = NULL,
                                   min_count_sum = 0,
                                   lfc_shrink = TRUE,
                                   shrink_type = "ashr",
                                   agg_args = list(assays = "RNA", return.seurat = TRUE),
                                   out_dir = NULL,
                                   existing_results = list(),
                                   min_cells = 2) {
  
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required for this function. Please install it.")
  }
  
  # ---------------------------------------------------------------------------
  # 1. Pseudobulk Aggregation & Save Aggregated Counts (if group_by is given)
  # ---------------------------------------------------------------------------
  seurat_obj <- object
  aggregated_counts <- NULL
  
  if (!is.null(group_by)) {
    message(">>> Running AggregateExpression to build Pseudobulk object...")
    
    # Force return.seurat = FALSE first to extract raw aggregated matrix if needed
    agg_params_mat <- c(list(object = seurat_obj, group.by = group_by, return.seurat = FALSE), agg_args)
    agg_params_mat <- agg_params_mat[!duplicated(names(agg_params_mat))]
    aggregated_counts <- do.call(Seurat::AggregateExpression, agg_params_mat)
    
    # Create aggregated Seurat object
    agg_params_obj <- c(list(object = seurat_obj, group.by = group_by, return.seurat = TRUE), agg_args)
    agg_params_obj <- agg_params_obj[!duplicated(names(agg_params_obj))]
    seurat_obj <- do.call(Seurat::AggregateExpression, agg_params_obj)
    
    # Save aggregated counts table to output directory if specified
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
  # 2. Determine Subsets (Stratified vs Global)
  # ---------------------------------------------------------------------------
  if (!is.null(split_col)) {
    all_splits <- unique(as.character(seurat_obj[[split_col]][, 1]))
    splits_to_run <- if (!is.null(target_splits)) intersect(target_splits, all_splits) else all_splits
  } else {
    splits_to_run <- "GLOBAL"
  }
  
  # ---------------------------------------------------------------------------
  # 3. Native DESeq2 Workflow Iteration
  # ---------------------------------------------------------------------------
  for (sp in splits_to_run) {
    
    # Subset Seurat object and metadata
    if (sp == "GLOBAL") {
      sub_obj <- seurat_obj
      message("\n>>> Running Global DESeq2 Analysis [Column: ", compare_col, "]")
    } else {
      cells_keep <- colnames(seurat_obj)[seurat_obj[[split_col]][, 1] == sp]
      sub_obj <- subset(seurat_obj, cells = cells_keep)
      message("\n>>> Processing Subset [", split_col, " = ", sp, "]")
    }
    
    meta_df <- sub_obj@meta.data
    meta_df[[compare_col]] <- as.factor(meta_df[[compare_col]])
    
    # Validate group sizes
    counts_table <- table(meta_df[[compare_col]])
    valid_groups <- names(counts_table[counts_table >= min_cells])
    
    if (length(valid_groups) < 2) {
      warning("Skipping ", sp, ": Less than 2 groups meet the minimum sample threshold (>= ", min_cells, ").")
      next
    }
    
    # Extract raw counts matrix (compatible with Seurat v4 and v5)
    assay_use <- Seurat::DefaultAssay(sub_obj)
    counts_matrix <- tryCatch({
      Seurat::GetAssayData(sub_obj, assay = assay_use, layer = "counts")
    }, error = function(e) {
      Seurat::GetAssayData(sub_obj, assay = assay_use, slot = "counts")
    })
    
    counts_matrix <- as.matrix(counts_matrix)
    mode(counts_matrix) <- "integer"
    
    # Construct Design Formula
    if (!is.null(covariates) && length(covariates) > 0) {
      design_formula <- stats::as.formula(paste("~", paste(c(covariates, compare_col), collapse = " + ")))
      for (cov in covariates) {
        meta_df[[cov]] <- as.factor(meta_df[[cov]])
      }
    } else {
      design_formula <- stats::as.formula(paste("~", compare_col))
    }
    
    # Define pairwise comparison pairs
    curr_comps <- list()
    if (!is.null(comparisons)) {
      for (comp in comparisons) {
        if (all(comp %in% valid_groups)) {
          curr_comps[[paste(comp[1], "vs", comp[2], sep = "_")]] <- comp
        }
      }
    } else {
      combos <- utils::combn(valid_groups, 2, simplify = FALSE)
      for (cb in combos) {
        curr_comps[[paste(cb[1], "vs", cb[2], sep = "_")]] <- cb
      }
    }
    
    # -------------------------------------------------------------------------
    # Build DESeqDataSet
    # -------------------------------------------------------------------------
    tryCatch({
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts_matrix,
        colData = meta_df,
        design = design_formula
      )
      
      # Optional gene pre-filtering based on total counts
      if (min_count_sum > 0) {
        keep_genes <- rowSums(DESeq2::counts(dds)) >= min_count_sum
        dds <- dds[keep_genes, ]
        message("  -> Applied gene filter (total count >= ", min_count_sum, "): ", sum(keep_genes), " / ", length(keep_genes), " genes retained.")
      }
      
      # Run DESeq pipeline
      dds <- DESeq2::DESeq(dds, quiet = TRUE)
      
      # Extract Normalized Counts Matrix
      norm_counts <- DESeq2::counts(dds, normalized = TRUE)
      norm_counts_df <- as.data.frame(norm_counts)
      colnames(norm_counts_df) <- paste0("norm_count_", colnames(norm_counts_df))
      norm_counts_df <- tibble::rownames_to_column(norm_counts_df, var = "gene")
      
      # Process each contrast
      for (comp_name in names(curr_comps)) {
        pair <- curr_comps[[comp_name]]
        group1 <- pair[1] # Numerator
        group2 <- pair[2] # Denominator/Control
        
        contrast_vec <- c(compare_col, group1, group2)
        contrast_label <- paste0(group1, "_vs_", group2, " (ref: ", group2, ")")
        
        message("  -> Extracting DESeq2 Contrast: ", group1, " vs ", group2)
        
        # Extract DESeq2 results
        res_deseq <- DESeq2::results(dds, contrast = contrast_vec)
        
        # Apply LFC Shrinkage if enabled
        if (lfc_shrink) {
          if (shrink_type == "ashr" && requireNamespace("ashr", quietly = TRUE)) {
            res_deseq <- DESeq2::lfcShrink(dds, contrast = contrast_vec, res = res_deseq, type = "ashr", quiet = TRUE)
          } else {
            res_deseq <- DESeq2::lfcShrink(dds, contrast = contrast_vec, res = res_deseq, type = "normal", quiet = TRUE)
          }
        }
        
        res_df <- as.data.frame(res_deseq)
        res_df <- tibble::rownames_to_column(res_df, var = "gene")
        
        # Add metadata annotation columns
        res_df <- res_df %>%
          dplyr::mutate(
            Comparison = comp_name,
            Contrast_Details = contrast_label,
            Group1_Numerator = group1,
            Group2_Denominator = group2
          )
        
        if (sp != "GLOBAL") {
          res_df[[split_col]] <- sp
        }
        
        # Merge Normalized Counts into the DE table
        final_de_df <- dplyr::left_join(res_df, norm_counts_df, by = "gene")
        
        # Save to result list
        list_key <- if (sp == "GLOBAL") paste(compare_col, comp_name, sep = "___") else paste(sp, comp_name, sep = "___")
        results_list[[list_key]] <- final_de_df
        all_combined_df[[list_key]] <- final_de_df
        
        # Output individual CSV
        if (!is.null(out_dir)) {
          prefix <- if (sp == "GLOBAL") "DESeq2_GLOBAL" else paste0("DESeq2_", sp)
          file_path <- file.path(out_dir, paste0(prefix, "_", comp_name, ".csv"))
          utils::write.csv(final_de_df, file = file_path, row.names = FALSE)
        }
      }
      
    }, error = function(e) {
      message("   [Error] DESeq2 processing failed for subset (", sp, "): ", e$message)
    })
  }
  
  # ---------------------------------------------------------------------------
  # 4. Export Combined Integrated Summary Table
  # ---------------------------------------------------------------------------
  if (!is.null(out_dir) && length(all_combined_df) > 0) {
    merged_df <- dplyr::bind_rows(all_combined_df)
    utils::write.csv(merged_df, file = file.path(out_dir, "ALL_combined_native_DESeq2.csv"), row.names = FALSE)
    message("\n>>> Integrated DESeq2 results written to: ", file.path(out_dir, "ALL_combined_native_DESeq2.csv"))
  }
  
  return(results_list)
}