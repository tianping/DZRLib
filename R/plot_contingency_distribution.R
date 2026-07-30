#' Plot Contingency Distribution with Custom Palettes and Multilevel Ordering
#'
#' @description
#' This function generates proportion and/or count barplots from a \code{Seurat} object, 
#' a \code{data.frame}, or a \code{table}. It supports advanced nested metadata sorting, 
#' manual level ordering, custom color palettes, automatic descriptive file naming,
#' and multi-panel combined layout via \code{patchwork}.
#' All external calls are explicitly namespaced to prevent function masking.
#'
#' @param data A \code{Seurat} object, a \code{data.frame}, or a \code{table}.
#' @param x_var Character. The column/metadata name for the X-axis (groups). Required for \code{Seurat} and \code{data.frame} inputs. Optional for \code{table}.
#' @param fill_var Character. The column/metadata name for the fill color (categories). Required for \code{Seurat} and \code{data.frame} inputs. Optional for \code{table}.
#' @param plot_type Character. One of \code{"merge"} (default), \code{"both"}, \code{"proportion"}, or \code{"count"}.
#' @param position Character. Bar positioning: \code{"stack"} (default) or \code{"dodge"}.
#' @param x_order Character vector. Custom explicit order for X-axis levels. Elements omitted here will be filtered out with a warning. Works for all input types.
#' @param fill_order Character vector. Custom explicit order for fill levels. Elements omitted here will be filtered out with a warning. Works for all input types.
#' @param x_sort_by Character vector. Rules to sort X-axis groups. Options: 
#'   \code{"none"} (default), \code{"Total"} (total counts), \code{"Prop"} (group proportion), 
#'   or a character vector of one or more metadata/column names in \code{data} to perform multilevel nested sorting (e.g., \code{c("Region", "BioRep")}).
#' @param decreasing Logical. Used when \code{x_sort_by} is active. Controls sorting direction. Default is \code{FALSE} (ascending).
#' @param palette Character vector. Custom colors for the fill categories. Can be a named vector mapping categories to colors (e.g., \code{c("A" = "red", "B" = "blue")}) or an unnamed vector of colors. If \code{NULL}, default colors are used.
#' @param x_text_angle Numeric. Angle of X-axis text labels. Default is \code{45}.
#' @param file_format Character. Output file extension, e.g., \code{"pdf"} (default), \code{"png"}, \code{"jpeg"}.
#' @param width Numeric. Width of the output plot in inches. If \code{NULL}, automatically calculated.
#' @param height Numeric. Height of the output plot in inches. If \code{NULL}, automatically calculated.
#' @param out_dir Character. Directory to save the plots. Default is \code{"."}.
#' @param file_prefix Character. Custom prefix for the output files. If \code{NULL}, auto-generated.
#'
#' @return Invisible \code{NULL}. Plots are saved directly to disk.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' # Example: Using a Seurat Object directly with multilevel sorting and custom colors
#' my_palette <- c("Excitatory" = "#3182bd", "Inhibitory" = "#31a354", "Astro" = "#9e9ac8")
#' plot_contingency_distribution(
#'   data = merged,
#'   x_var = "FinalName",
#'   fill_var = "predicted.cell_subclass",
#'   x_sort_by = c("Region", "BioRep"),
#'   decreasing = FALSE,
#'   palette = my_palette
#' )
plot_contingency_distribution <- function(
    data,
    x_var = NULL,
    fill_var = NULL,
    plot_type = c("merge", "both", "proportion", "count"),
    position = c("stack", "dodge"),
    x_order = NULL,
    fill_order = NULL,
    x_sort_by = "none",
    decreasing = FALSE,
    palette = NULL,
    x_text_angle = 45,
    file_format = "pdf",
    width = NULL,
    height = NULL,
    out_dir = ".",
    file_prefix = NULL
) {
  
  # 1. Validate arguments
  plot_type <- match.arg(plot_type)
  position <- match.arg(position)
  
  # 2. Standardize input data (Support Seurat, Data.frame, and Table)
  plot_df <- NULL
  guessed_x <- "X"
  guessed_fill <- "Fill"
  is_table_input <- inherits(data, "table")
  is_seurat_input <- inherits(data, "Seurat")
  
  # Temporary data frame to store full metadata for Seurat/data.frame column sorting
  meta_df <- NULL 
  
  if (is_seurat_input) {
    if (is.null(x_var) || is.null(fill_var)) {
      stop("When 'data' is a Seurat object, both 'x_var' and 'fill_var' must be specified.")
    }
    
    # Identify all required metadata columns for plotting and sorting
    required_cols <- unique(c(x_var, fill_var))
    if (!identical(x_sort_by, "none") && !any(x_sort_by %in% c("Total", "Prop"))) {
      required_cols <- unique(c(required_cols, x_sort_by))
    }
    
    # Retrieve data safely
    meta_df <- as.data.frame(Seurat::FetchData(data, vars = required_cols))
    
    # Build standard table frequency
    plot_df <- as.data.frame(table(meta_df[[x_var]], meta_df[[fill_var]]))
    colnames(plot_df) <- c("X_Group", "Fill_Category", "Count")
    
    guessed_x <- x_var
    guessed_fill <- fill_var
    internal_x <- "X_Group"
    internal_fill <- "Fill_Category"
    internal_weight <- "Count"
    
  } else if (inherits(data, "data.frame")) {
    if (is.null(x_var) || !x_var %in% colnames(data)) {
      stop("When 'data' is a data.frame, a valid 'x_var' must be provided.")
    }
    if (is.null(fill_var) || !fill_var %in% colnames(data)) {
      stop("When 'data' is a data.frame, a valid 'fill_var' must be provided.")
    }
    
    meta_df <- as.data.frame(data)
    
    # Standardize data.frame input by computing a frequency table
    plot_df <- as.data.frame(table(meta_df[[x_var]], meta_df[[fill_var]]))
    colnames(plot_df) <- c("X_Group", "Fill_Category", "Count")
    
    guessed_x <- x_var
    guessed_fill <- fill_var
    internal_x <- "X_Group"
    internal_fill <- "Fill_Category"
    internal_weight <- "Count"
    
  } else if (is_table_input) {
    tab_names <- names(dimnames(data))
    if (!is.null(tab_names) && length(tab_names) >= 2) {
      if (nzchar(tab_names[1])) guessed_x <- tab_names[1]
      if (nzchar(tab_names[2])) guessed_fill <- tab_names[2]
    }
    
    if (guessed_x == "X" && guessed_fill == "Fill") {
      raw_expr <- deparse(substitute(data))
      match_res <- regmatches(raw_expr, regexec("table\\(([^,]+),\\s*([^)]+)\\)", raw_expr))[[1]]
      if (length(match_res) >= 3) {
        guessed_x <- gsub("^.*\\$", "", trimws(match_res[2]))
        guessed_fill <- gsub("^.*\\$", "", trimws(match_res[3]))
      }
    }
    
    if (!is.null(x_var)) guessed_x <- x_var
    if (!is.null(fill_var)) guessed_fill <- fill_var
    
    plot_df <- as.data.frame(data)
    colnames(plot_df) <- c("X_Group", "Fill_Category", "Count")
    
    internal_x <- "X_Group"
    internal_fill <- "Fill_Category"
    internal_weight <- "Count"
    
  } else {
    stop("Input 'data' must be a Seurat object, a data.frame, or a table.")
  }
  
  # Ensure character conversion before factor manipulation
  plot_df[[internal_x]] <- as.character(plot_df[[internal_x]])
  plot_df[[internal_fill]] <- as.character(plot_df[[internal_fill]])
  
  # 3. Method 1: Manual Reordering and Filtering
  if (!is.null(x_order)) {
    actual_x_elements <- unique(plot_df[[internal_x]])
    omitted_x <- setdiff(actual_x_elements, x_order)
    if (length(omitted_x) > 0) {
      warning("The following elements in x_var were omitted based on 'x_order': ", 
              paste(omitted_x, collapse = ", "))
    }
    plot_df <- plot_df[plot_df[[internal_x]] %in% x_order, ]
    plot_df[[internal_x]] <- factor(plot_df[[internal_x]], levels = x_order)
  }
  
  if (!is.null(fill_order)) {
    actual_fill_elements <- unique(plot_df[[internal_fill]])
    omitted_fill <- setdiff(actual_fill_elements, fill_order)
    if (length(omitted_fill) > 0) {
      warning("The following elements in fill_var were omitted based on 'fill_order': ", 
              paste(omitted_fill, collapse = ", "))
    }
    plot_df <- plot_df[plot_df[[internal_fill]] %in% fill_order, ]
    plot_df[[internal_fill]] <- factor(plot_df[[internal_fill]], levels = fill_order)
  }
  
  # Remove zero-count rows after potential subsetting
  plot_df <- plot_df[plot_df[[internal_weight]] > 0, ]
  if (nrow(plot_df) == 0) {
    stop("No data remaining after applying the custom ordering filters.")
  }
  
  # 4. Calculation of baseline statistics
  df_sum <- stats::aggregate(plot_df[[internal_weight]], by = list(plot_df[[internal_x]]), FUN = sum)
  colnames(df_sum) <- c(internal_x, "Total")
  
  plot_df_prop <- merge(plot_df, df_sum, by = internal_x)
  plot_df_prop$Prop <- plot_df_prop[[internal_weight]] / plot_df_prop$Total
  
  # Method 2: Automatic or Multilevel Column Sorting
  if (is.null(x_order) && !identical(x_sort_by, "none")) {
    
    if (length(x_sort_by) == 1 && x_sort_by == "Total") {
      ordered_x_groups <- df_sum[order(df_sum$Total, decreasing = decreasing), internal_x]
      plot_df[[internal_x]] <- factor(plot_df[[internal_x]], levels = ordered_x_groups)
      plot_df_prop[[internal_x]] <- factor(plot_df_prop[[internal_x]], levels = ordered_x_groups)
      
    } else if (length(x_sort_by) == 1 && x_sort_by == "Prop") {
      first_fill_cat <- if (!is.null(fill_order)) fill_order[1] else unique(plot_df_prop[[internal_fill]])[1]
      sub_prop <- plot_df_prop[plot_df_prop[[internal_fill]] == first_fill_cat, ]
      ordered_x_groups <- sub_prop[order(sub_prop$Prop, decreasing = decreasing), internal_x]
      
      all_groups <- unique(plot_df_prop[[internal_x]])
      ordered_x_groups <- c(ordered_x_groups, setdiff(all_groups, ordered_x_groups))
      
      plot_df[[internal_x]] <- factor(plot_df[[internal_x]], levels = ordered_x_groups)
      plot_df_prop[[internal_x]] <- factor(plot_df_prop[[internal_x]], levels = ordered_x_groups)
      
    } else {
      # Multilevel sorting based on metadata columns
      if (is_table_input) {
        warning("Multilevel column sorting is not supported for 'table' inputs. Falling back to default ordering.")
        ordered_x_groups <- unique(plot_df[[internal_x]])
      } else {
        # Ensure all requested sorting columns exist
        missing_cols <- setdiff(x_sort_by, colnames(meta_df))
        if (length(missing_cols) > 0) {
          stop("The following sorting columns were not found in the input metadata: ", 
               paste(missing_cols, collapse = ", "))
        }
        
        # Select unique mappings
        meta_sub <- unique(meta_df[, c(x_var, x_sort_by), drop = FALSE])
        
        # Prepare sorting arguments
        sort_args <- lapply(x_sort_by, function(col) {
          val <- meta_sub[[col]]
          if (is.character(val)) val <- factor(val)
          return(val)
        })
        
        if (decreasing) {
          sort_args <- lapply(sort_args, function(x) {
            if (is.numeric(x)) -x else -xtfrm(x)
          })
        }
        
        sort_order <- do.call(order, sort_args)
        ordered_x_groups <- as.character(meta_sub[[x_var]][sort_order])
        ordered_x_groups <- ordered_x_groups[ordered_x_groups %in% unique(plot_df[[internal_x]])]
      }
      
      plot_df[[internal_x]] <- factor(plot_df[[internal_x]], levels = ordered_x_groups)
      plot_df_prop[[internal_x]] <- factor(plot_df_prop[[internal_x]], levels = ordered_x_groups)
    }
  } else {
    # Default levels
    if (!is.factor(plot_df[[internal_x]])) {
      plot_df[[internal_x]] <- factor(plot_df[[internal_x]])
      plot_df_prop[[internal_x]] <- factor(plot_df_prop[[internal_x]], levels = levels(plot_df[[internal_x]]))
    }
    if (!is.factor(plot_df[[internal_fill]])) {
      plot_df[[internal_fill]] <- factor(plot_df[[internal_fill]])
      plot_df_prop[[internal_fill]] <- factor(plot_df_prop[[internal_fill]], levels = levels(plot_df[[internal_fill]]))
    }
  }
  
  # 5. Handle file naming and metadata prefixes dynamically
  if (is.null(file_prefix)) {
    file_prefix <- paste0(guessed_x, "_vs_", guessed_fill)
  }
  file_format <- gsub("^\\.", "", tolower(file_format))
  
  # 6. Dynamic width and height estimation
  num_x_groups <- length(unique(plot_df[[internal_x]]))
  if (is.null(width)) {
    factor <- if (position == "dodge") 0.8 else 0.6
    width <- 3 + (num_x_groups * factor)
  }
  if (is.null(height)) {
    height <- if (plot_type %in% c("merge", "both")) 8 else 5
  }
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # 7. Theme Layer Setup with FULL explicit namespacing
  v_just <- if (x_text_angle == 90) 0.5 else 1
  h_just <- if (x_text_angle == 0) 0.5 else 1
  
  base_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = x_text_angle, hjust = h_just, vjust = v_just, color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  lbl_x <- guessed_x
  lbl_fill <- guessed_fill
  
  # 8. Define Custom Palette Scale
  fill_scale <- if (!is.null(palette)) {
    ggplot2::scale_fill_manual(values = palette)
  } else {
    ggplot2::scale_fill_discrete()
  }
  
  # 9. Core Plotting Logic
  
  # --- Plot 1: Proportion Plot ---
  if (plot_type %in% c("merge", "both", "proportion")) {
    if (position == "dodge") {
      p1 <- ggplot2::ggplot(plot_df_prop, ggplot2::aes(x = .data[[internal_x]], y = .data$Prop, fill = .data[[internal_fill]])) +
        ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(preserve = "single")) +
        ggplot2::scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1))
    } else {
      p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[internal_x]], fill = .data[[internal_fill]])) +
        ggplot2::geom_bar(ggplot2::aes(weight = .data[[internal_weight]]), position = ggplot2::position_fill()) +
        ggplot2::scale_y_continuous(labels = scales::percent_format())
    }
    
    p1 <- p1 + base_theme + fill_scale + ggplot2::labs(title = "Proportion Distribution", x = lbl_x, y = "Proportion", fill = lbl_fill)
    
    if (plot_type == "proportion") {
      ggplot2::ggsave(file.path(out_dir, paste0(file_prefix, "_proportion.barplot.", file_format)), plot = p1, width = width, height = height, dpi = 300)
    }
  }
  
  # --- Plot 2: Count Plot ---
  if (plot_type %in% c("merge", "both", "count")) {
    p2 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[internal_x]], fill = .data[[internal_fill]]))
    bar_pos <- if (position == "dodge") ggplot2::position_dodge(preserve = "single") else "stack"
    p2 <- p2 + ggplot2::geom_bar(ggplot2::aes(weight = .data[[internal_weight]]), position = bar_pos)
    
    p2 <- p2 + base_theme + fill_scale + ggplot2::labs(title = "Count Distribution", x = lbl_x, y = "Cell Count", fill = lbl_fill)
    
    if (plot_type == "count") {
      ggplot2::ggsave(file.path(out_dir, paste0(file_prefix, "_count.barplot.", file_format)), plot = p2, width = width, height = height, dpi = 300)
    }
  }
  
  # --- Plot 3: Combined / Merged Plot Options ---
  if (plot_type == "merge") {
    # Hide X-axis elements for the top plot (count plot) to share X-axis with the bottom plot
    p2_top <- p2 + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )
    
    # Combine plots vertically: Count on top, Proportion on bottom
    p_merged <- (p2_top / p1) + patchwork::plot_layout(guides = "collect")
    
    ggplot2::ggsave(file.path(out_dir, paste0(file_prefix, "_merged.barplot.", file_format)), plot = p_merged, width = width, height = height, dpi = 300)
    
  } else if (plot_type == "both") {
    ggplot2::ggsave(file.path(out_dir, paste0(file_prefix, "_proportion.barplot.", file_format)), plot = p1, width = width, height = height, dpi = 300)
    ggplot2::ggsave(file.path(out_dir, paste0(file_prefix, "_count.barplot.", file_format)), plot = p2, width = width, height = height, dpi = 300)
  }
  
  invisible(NULL)
}

# Internal helper function
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}