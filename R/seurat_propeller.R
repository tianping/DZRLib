#' Run Propeller Proportion Analysis and Plots on Seurat Objects
#'
#' @description
#' Runs differential cell type proportion analysis on a Seurat object using \code{speckle::propeller} 
#' and generates cell proportion distribution plots using \code{plot_contingency_distribution}. 
#' Supports splitting the dataset by metadata factors (e.g., region or tissue), filtering specific 
#' subsets, running multiple transformations, and exporting cleanly formatted TSV tables and PDF plots.
#'
#' @param object A Seurat object containing single-cell transcriptomics metadata.
#' @param propeller_clusters Character string. Meta.data column name representing cell clusters or types 
#'   (passed to \code{speckle::propeller} as \code{clusters}). Default is \code{"cell_IntGroup"}.
#' @param propeller_sample Character string. Meta.data column name representing individual biological samples 
#'   (passed to \code{speckle::propeller} as \code{sample}). Default is \code{"FinalName"}.
#' @param propeller_group Character string. Meta.data column name representing experimental groups 
#'   (passed to \code{speckle::propeller} as \code{group}). Default is \code{"Group"}.
#' @param propeller_trend Logical. Passed to \code{speckle::propeller} as \code{trend}. Whether to fit 
#'   a mean-variance trend. Default is \code{FALSE}.
#' @param propeller_robust Logical. Passed to \code{speckle::propeller} as \code{robust}. Whether to use 
#'   robust empirical Bayes moderation. Default is \code{TRUE}.
#' @param propeller_transform Character string. Passed to \code{speckle::propeller} as \code{transform}. 
#'   Transformation method for proportions (\code{"both"}, \code{"asin"}, \code{"logit"}, or \code{"alr"}). 
#'   If set to \code{"both"}, runs both \code{"asin"} and \code{"logit"} transformations and exports 
#'   separate TSV files for each. Default is \code{"both"}.
#' @param split_by Character string. Meta.data column name used to subset the dataset (e.g., \code{"Region"}). 
#'   If \code{NULL} or column does not exist, analysis runs on the entire dataset. Default is \code{"Region"}.
#' @param split_select Character vector. Optional subset of values from \code{split_by} to process. 
#'   If \code{NULL}, all unique non-NA values in \code{split_by} will be processed. Default is \code{NULL}.
#' @param out_dir Character string. Directory path where output TSV tables and PDF plot files will be saved. 
#'   Default is \code{"."}.
#' @param file_prefix Character string. Prefix prepended to ALL output files (both TSV tables and PDF plots). 
#'   Default is \code{NULL}.
#' @param plot_width Numeric. PDF plot width in inches. Default is \code{8}.
#' @param plot_height Numeric. PDF plot height in inches. Default is \code{6}.
#' @param ... Additional arguments passed to \code{plot_contingency_distribution}.
#'
#' @return A nested list containing \code{propeller} result data frames for each split group and 
#'   transformation method, returned invisibly.
#'
#' @import Seurat
#' @import speckle
#' @import ggplot2
#' @export
seurat_propeller <- function(
  object,
  propeller_clusters  = "cell_IntGroup",
  propeller_sample    = "FinalName",
  propeller_group     = "Group",
  propeller_trend     = FALSE,
  propeller_robust    = TRUE,
  propeller_transform = "both",
  split_by            = "Region",
  split_select        = NULL,
  out_dir             = ".",
  file_prefix         = NULL,
  plot_width          = 8,
  plot_height         = 6,
  ...
) {
  # Ensure the output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Helper function to process a single Seurat dataset/subset
  process_dataset <- function(sub_obj, curr_split_val = NULL) {
    # Combine file prefix and current split category
    prefix_parts <- c(file_prefix, curr_split_val)
    prefix_parts <- prefix_parts[!is.null(prefix_parts) & prefix_parts != ""]

    if (length(prefix_parts) > 0) {
      combined_prefix <- paste(prefix_parts, collapse = "_")
      file_prefix_str <- paste0(combined_prefix, "_")
    } else {
      combined_prefix <- NULL
      file_prefix_str <- ""
    }

    split_log_str <- if (!is.null(curr_split_val)) paste0("Subset [", curr_split_val, "] ") else ""

    # --- 1. Quality Control & Sample Checks ---
    # Check 1: Ensure at least 2 groups exist
    groups_present <- unique(sub_obj[[propeller_group, drop = TRUE]][!is.na(sub_obj[[propeller_group, drop = TRUE]])])
    if (length(groups_present) < 2) {
      warning(paste0(split_log_str, "contains fewer than 2 groups (", 
                     paste(groups_present, collapse = ", "), "). Skipping analysis."))
      return(NULL)
    }

    # Check 2: Ensure at least one group has >= 2 sample replicates
    sample_vec <- sub_obj[[propeller_sample, drop = TRUE]]
    group_vec  <- sub_obj[[propeller_group, drop = TRUE]]
    reps_per_group <- tapply(sample_vec, group_vec, function(x) length(unique(x)))
    reps_per_group <- reps_per_group[!is.na(reps_per_group)]

    if (!any(reps_per_group >= 2)) {
      warning(paste0(split_log_str, "does not have any group with >= 2 sample replicates. Skipping analysis."))
      return(NULL)
    }

    # Set active ident for Seurat object
    Idents(sub_obj) <- propeller_clusters

    # --- 2. Determine transformations to run ---
    if (identical(propeller_transform, "both")) {
      transforms_to_run <- c("asin", "logit")
    } else {
      transforms_to_run <- propeller_transform
    }

    prop_results <- list()

    # --- 3. Run propeller Differential Proportion Analysis ---
    for (tr in transforms_to_run) {
      prop_res <- speckle::propeller(
        sub_obj,
        clusters  = sub_obj[[propeller_clusters, drop = TRUE]],
        sample    = sub_obj[[propeller_sample, drop = TRUE]],
        group     = sub_obj[[propeller_group, drop = TRUE]],
        trend     = propeller_trend,
        robust    = propeller_robust,
        transform = tr
      )

      # Convert row names to an explicit "Cluster" column to avoid header shifting in TSV
      export_df <- data.frame(Cluster = rownames(prop_res), prop_res, check.names = FALSE)

      # Save tab-delimited TSV file
      tsv_filename <- paste0(file_prefix_str, "propeller_results_", tr, ".tsv")
      tsv_filepath <- file.path(out_dir, tsv_filename)
      
      write.table(
        export_df, 
        file      = tsv_filepath, 
        sep       = "\t", 
        quote     = FALSE, 
        row.names = FALSE, 
        col.names = TRUE
      )
      message(paste0("  - Saved propeller TSV [", tr, "]: ", tsv_filepath))

      prop_results[[tr]] <- prop_res
    }

    # --- 4. Generate and Save Distribution Plot ---
    plot_contingency_distribution(
      sub_obj,
      x_var       = propeller_sample,
      fill_var    = propeller_clusters,
      out_dir     = out_dir,
      file_prefix = combined_prefix,
      width       = plot_width,
      height      = plot_height,
      ...
    )

    message(paste0("  - Completed plot generation with prefix: '", 
                   ifelse(is.null(combined_prefix), "(none)", combined_prefix), "'"))

    if (length(prop_results) == 1) {
      return(prop_results[[1]])
    } else {
      return(prop_results)
    }
  }

  results_list <- list()

  # Determine if dataset should be split by split_by factor
  has_split_col <- !is.null(split_by) && (split_by %in% colnames(object@meta.data))

  if (has_split_col) {
    unique_splits <- unique(object[[split_by, drop = TRUE]])
    unique_splits <- unique_splits[!is.na(unique_splits)]

    # Filter to split_select if user provided specific categories
    if (!is.null(split_select)) {
      unique_splits <- intersect(unique_splits, split_select)
      if (length(unique_splits) == 0) {
        warning("None of the values in 'split_select' were found in column '", split_by, "'.")
        return(invisible(NULL))
      }
    }

    for (s_val in unique_splits) {
      s_str <- as.character(s_val)
      message(paste0("\n=== Processing Split [", split_by, "]: ", s_str, " ==="))
      sub_object <- object[, object[[split_by, drop = TRUE]] == s_val]

      res <- tryCatch({
        process_dataset(sub_object, curr_split_val = s_str)
      }, error = function(e) {
        warning(paste0("Error encountered in subset [", s_str, "]: ", e$message))
        return(NULL)
      })

      results_list[[s_str]] <- res
    }
  } else {
    message("\n=== 'split_by' variable not specified or missing; processing full dataset ===")
    res <- tryCatch({
      process_dataset(object, curr_split_val = NULL)
    }, error = function(e) {
      warning(paste0("Error encountered during analysis: ", e$message))
      return(NULL)
    })
    results_list[["all"]] <- res
  }

  invisible(results_list)
}
