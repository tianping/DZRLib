#' Predict sample sex based on marker gene expression in scRNA-seq data
#'
#' This function pools single-cell data by individual, calculates CPM values,
#' and uses male/female marker genes to predict the sex of each sample.
#'
#' @param object A Seurat object
#' @param group.by Character string specifying the metadata column that identifies
#'   individual samples. Default: "orig.ident"
#' @param male_markers Character vector of male marker genes. 
#'   Default: c("RPS4Y1", "RPS4Y2", "DDX3Y")
#' @param female_markers Character vector of female marker genes.
#'   Default: c("XIST", "LOC106995245")
#' @param assay Character string specifying which assay to use. Default: "RNA"
#' @param layer Character string specifying which layer to use (for Seurat v5). 
#'   Default: "counts"
#' @param heatmap_log Logical; whether to use log2(CPM+1) for heatmap. Default: TRUE
#' @param heatmap_scale_by_gene Logical; whether to scale by row (gene) in heatmap. 
#'   Default: TRUE
#' @param scatter_log Logical; whether to use log2(CPM+1) for scatter plot. Default: TRUE
#' @param sex_ratio_threshold Numeric; threshold for male/female expression ratio.
#'   Samples with ratio > threshold are classified as male, < 1/threshold as female.
#'   Default: 2
#' @param output_dir Character string; directory to save output files. If NULL,
#'   no files are saved. Default: NULL
#' @param verbose Logical; whether to print progress messages. Default: TRUE
#'
#' @return A list containing:
#'   \item{metadata}{Data frame with sample metadata, marker gene CPM values, and PredictedSex}
#'   \item{pseudo_bulk_cpm}{Matrix of CPM values for marker genes (pooled by sample)}
#'   \item{scatter_plot}{ggplot object showing average male vs female marker expression}
#'   \item{heatmap}{ggplot object showing expression of all marker genes across samples}
#'   \item{predicted_sex}{Named character vector of predicted sex for each sample}
#'
#' @examples
#' \dontrun{
#' result <- seurat_decide_sex(
#'   object = seurat_obj,
#'   group.by = "donor_id",
#'   output_dir = "results/sex_prediction"
#' )
#' }
#'
#' @import Seurat
#' @import ggplot2
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra grid.arrange
#' @export
seurat_decide_sex <- function(
  object,
  group.by = "orig.ident",
  male_markers = c("RPS4Y1", "RPS4Y2", "DDX3Y"),
  female_markers = c("XIST", "LOC106995245"),
  assay = "RNA",
  layer = "counts",
  heatmap_log = TRUE,
  heatmap_scale_by_gene = TRUE,
  scatter_log = TRUE,
  sex_ratio_threshold = 2,
  output_dir = NULL,
  verbose = TRUE
) {
  
  # ==================== Input Validation ====================
  
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  
  if (!group.by %in% colnames(object@meta.data)) {
    stop("group.by column '", group.by, "' not found in object metadata")
  }
  
  # Check Seurat version for compatibility
  seurat_version <- packageVersion("SeuratObject")
  use_layer <- seurat_version >= "5.0.0"
  
  # Check if markers exist in the object
  all_markers <- c(male_markers, female_markers)
  markers_in_data <- all_markers[all_markers %in% rownames(object)]
  missing_markers <- all_markers[!all_markers %in% rownames(object)]
  
  if (length(missing_markers) > 0) {
    warning("The following markers are not found in the object: ", 
            paste(missing_markers, collapse = ", "))
  }
  
  if (length(markers_in_data) == 0) {
    stop("No marker genes found in the object. Please check the marker gene names.")
  }
  
  # Filter to only available markers
  male_markers_avail <- male_markers[male_markers %in% rownames(object)]
  female_markers_avail <- female_markers[female_markers %in% rownames(object)]
  
  if (length(male_markers_avail) == 0) {
    warning("No male markers available. Sex prediction may be unreliable.")
  }
  if (length(female_markers_avail) == 0) {
    warning("No female markers available. Sex prediction may be unreliable.")
  }
  
  if (verbose) {
    message("Available male markers: ", paste(male_markers_avail, collapse = ", "))
    message("Available female markers: ", paste(female_markers_avail, collapse = ", "))
    message("Seurat version: ", packageVersion("SeuratObject"))
    message("Using layer: ", layer)
  }
  
  # ==================== Step 1: Pool by individual ====================
  
  if (verbose) message("Pooling scRNA-seq data by individual: ", group.by)
  
  # Get count matrix - compatible with both Seurat v4 and v5
  if (use_layer) {
    # Seurat v5: use layer argument
    counts_matrix <- GetAssayData(object, assay = assay, layer = layer)
  } else {
    # Seurat v4: use slot argument
    counts_matrix <- GetAssayData(object, assay = assay, slot = layer)
  }
  
  # Get sample identities
  sample_ids <- object@meta.data[[group.by]]
  unique_samples <- unique(sample_ids)
  
  if (verbose) message("Found ", length(unique_samples), " unique samples")
  
  # Pool counts by sample (sum across cells)
  pooled_counts <- matrix(
    0, 
    nrow = nrow(counts_matrix), 
    ncol = length(unique_samples),
    dimnames = list(rownames(counts_matrix), unique_samples)
  )
  
  for (sample in unique_samples) {
    cells_in_sample <- colnames(object)[sample_ids == sample]
    if (length(cells_in_sample) > 0) {
      pooled_counts[, sample] <- Matrix::rowSums(counts_matrix[, cells_in_sample, drop = FALSE])
    } else {
      warning("No cells found for sample: ", sample)
    }
  }
  
  # ==================== Step 2: Calculate CPM ====================
  
  if (verbose) message("Calculating CPM values (no gene length normalization)")
  
  # Standard CPM: counts / total_counts * 1e6
  total_counts_per_sample <- colSums(pooled_counts)
  cpm_matrix <- sweep(pooled_counts, 2, total_counts_per_sample, FUN = "/") * 1e6
  
  # ==================== Step 3: Extract marker gene CPM ====================
  
  marker_cpm <- cpm_matrix[markers_in_data, , drop = FALSE]
  
  # ==================== Step 4: Calculate average expression ====================
  
  if (verbose) message("Calculating average marker gene expression")
  
  # Calculate mean for male and female markers
  mean_male <- rep(NA, length(unique_samples))
  mean_female <- rep(NA, length(unique_samples))
  
  if (length(male_markers_avail) > 0) {
    male_cpm <- marker_cpm[male_markers_avail, , drop = FALSE]
    if (nrow(male_cpm) == 1) {
      mean_male <- as.vector(male_cpm)
    } else {
      mean_male <- colMeans(male_cpm)
    }
  }
  
  if (length(female_markers_avail) > 0) {
    female_cpm <- marker_cpm[female_markers_avail, , drop = FALSE]
    if (nrow(female_cpm) == 1) {
      mean_female <- as.vector(female_cpm)
    } else {
      mean_female <- colMeans(female_cpm)
    }
  }
  
  # ==================== Step 5: Predict sex ====================
  
  if (verbose) message("Predicting sex using ratio threshold: ", sex_ratio_threshold)
  
  predicted_sex <- rep("Unknown", length(unique_samples))
  names(predicted_sex) <- unique_samples
  
  for (i in seq_along(unique_samples)) {
    sample <- unique_samples[i]
    
    # Skip if markers are missing
    if (is.na(mean_male[i]) || is.na(mean_female[i])) {
      predicted_sex[i] <- "Unknown"
      next
    }
    
    # Use log2(CPM+1) for ratio if scatter_log is TRUE
    if (scatter_log) {
      male_val <- log2(mean_male[i] + 1)
      female_val <- log2(mean_female[i] + 1)
    } else {
      male_val <- mean_male[i]
      female_val <- mean_female[i]
    }
    
    # Predict sex based on ratio
    if (male_val > female_val * sex_ratio_threshold) {
      predicted_sex[i] <- "Male"
    } else if (female_val > male_val * sex_ratio_threshold) {
      predicted_sex[i] <- "Female"
    } else {
      predicted_sex[i] <- "Unknown"
    }
  }
  
  if (verbose) {
    message("Sex prediction summary:")
    print(table(predicted_sex))
  }
  
  # ==================== Step 6: Create metadata ====================
  
  metadata <- data.frame(
    Sample = unique_samples,
    MeanMaleCPM = mean_male,
    MeanFemaleCPM = mean_female,
    PredictedSex = predicted_sex,
    stringsAsFactors = FALSE
  )
  
  # Add individual marker gene CPM values to metadata
  for (gene in markers_in_data) {
    metadata[[paste0(gene, "_CPM")]] <- marker_cpm[gene, ]
  }
  
  rownames(metadata) <- unique_samples
  
  # ==================== Step 7: Create scatter plot ====================
  
  if (verbose) message("Creating scatter plot")
  
  # Prepare data for scatter plot
  scatter_data <- data.frame(
    Sample = unique_samples,
    MaleAvg = mean_male,
    FemaleAvg = mean_female,
    PredictedSex = predicted_sex,
    stringsAsFactors = FALSE
  )
  
  # Apply log transformation if requested
  if (scatter_log) {
    scatter_data$MaleAvg <- log2(scatter_data$MaleAvg + 1)
    scatter_data$FemaleAvg <- log2(scatter_data$FemaleAvg + 1)
    x_lab <- "log2(Mean Male Marker CPM + 1)"
    y_lab <- "log2(Mean Female Marker CPM + 1)"
  } else {
    x_lab <- "Mean Male Marker CPM"
    y_lab <- "Mean Female Marker CPM"
  }
  
  # Define colors for sex
  sex_colors <- c("Male" = "#2E86AB", "Female" = "#A23B72", "Unknown" = "#D3D3D3")
  
  scatter_plot <- ggplot2::ggplot(
    scatter_data, 
    aes(x = MaleAvg, y = FemaleAvg, color = PredictedSex, label = Sample)
  ) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
    scale_color_manual(values = sex_colors) +
    labs(
      title = "Sex Prediction Based on Marker Gene Expression",
      x = x_lab,
      y = y_lab,
      color = "Predicted Sex"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5)
  
  # Add threshold lines if using log scale
  if (scatter_log) {
    scatter_plot <- scatter_plot +
      geom_abline(
        intercept = log2(sex_ratio_threshold), 
        slope = 1, 
        color = "darkgreen", 
        linetype = "dotted", 
        alpha = 0.7
      ) +
      geom_abline(
        intercept = -log2(sex_ratio_threshold), 
        slope = 1, 
        color = "darkgreen", 
        linetype = "dotted", 
        alpha = 0.7
      )
  }
  
  # ==================== Step 8: Create heatmap ====================
  
  if (verbose) message("Creating heatmap")
  
  # Prepare data for heatmap
  heatmap_data <- marker_cpm
  
  # Apply log transformation if requested
  if (heatmap_log) {
    heatmap_data <- log2(heatmap_data + 1)
    heatmap_title <- "log2(CPM + 1) of Sex Marker Genes"
  } else {
    heatmap_title <- "CPM of Sex Marker Genes"
  }
  
  # Create annotation for samples
  annotation_df <- data.frame(
    PredictedSex = predicted_sex
  )
  rownames(annotation_df) <- unique_samples
  
  # Define annotation colors
  ann_colors <- list(
    PredictedSex = c("Male" = "#2E86AB", "Female" = "#A23B72", "Unknown" = "#D3D3D3")
  )
  
  # Create heatmap using pheatmap
  heatmap_obj <- pheatmap::pheatmap(
    mat = heatmap_data,
    scale = if (heatmap_scale_by_gene) "row" else "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = annotation_df,
    annotation_colors = ann_colors,
    main = heatmap_title,
    fontsize = 10,
    fontsize_row = 10,
    fontsize_col = 10,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    silent = TRUE
  )
  
  # Extract the heatmap grob
  heatmap_grob <- heatmap_obj[[4]]
  
  # ==================== Step 9: Save outputs if requested ====================
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    if (verbose) message("Saving outputs to: ", output_dir)
    
    # Save metadata
    write.csv(metadata, file = file.path(output_dir, "sex_prediction_metadata.csv"), row.names = FALSE)
    
    # Save marker CPM matrix
    write.csv(marker_cpm, file = file.path(output_dir, "marker_gene_cpm.csv"))
    
    # Save scatter plot
    ggsave(
      filename = file.path(output_dir, "sex_prediction_scatter.pdf"),
      plot = scatter_plot,
      width = 10,
      height = 8
    )
    
    # Save heatmap
    ggsave(
      filename = file.path(output_dir, "sex_prediction_heatmap.pdf"),
      plot = heatmap_grob,
      width = 10,
      height = 8
    )
    
    # Save sex prediction summary
    sink(file.path(output_dir, "sex_prediction_summary.txt"))
    cat("Sex Prediction Summary\n")
    cat("======================\n\n")
    cat("Threshold (ratio):", sex_ratio_threshold, "\n")
    cat("Male markers:", paste(male_markers_avail, collapse = ", "), "\n")
    cat("Female markers:", paste(female_markers_avail, collapse = ", "), "\n\n")
    cat("Sex prediction counts:\n")
    print(table(predicted_sex))
    cat("\n\nDetailed predictions:\n")
    print(metadata)
    sink()
    
    if (verbose) message("All outputs saved successfully")
  }
  
  # ==================== Step 10: Return results ====================
  
  return(list(
    metadata = metadata,
    pseudo_bulk_cpm = marker_cpm,
    scatter_plot = scatter_plot,
    heatmap = heatmap_grob,
    predicted_sex = predicted_sex
  ))
}
