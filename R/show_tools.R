#' DZRLib CLI Function Browser (v3.2 FINAL)
#'
#' Non-destructive exact-match function browser with optional fuzzy ranking.
#' Lists exported functions from DZRLib and supports category grouping.
#'
#' @param ... character strings used for searching (optional)
#' @param category optional filter: "seurat", "rna", "atac", "plotting", "utils"
#' @param fuzzy logical; whether to enable fuzzy ranking (default FALSE)
#' @param silent logical; if TRUE, suppress console output and return data.frame only
#'
#' @return data.frame (invisible) with columns:
#' func, category, match, score
#'
#' @export
#'
#' @examples
#' show_tools()
#' show_tools("seu")
#' show_tools("qc", fuzzy = TRUE)
#' show_tools(category = "seurat")
show_tools <- function(..., category = NULL, fuzzy = FALSE, silent = FALSE) {

  # ----------------------------
  # 1. namespace-safe export list
  # ----------------------------
  all_funcs <- getNamespaceExports("DZRLib")
  user_funcs <- all_funcs[grepl("_", all_funcs)]

  if (length(user_funcs) == 0) {
    message("No exported user functions found in DZRLib.")
    return(invisible(NULL))
  }

  # ----------------------------
  # 2. base table
  # ----------------------------
  df <- data.frame(
    func = user_funcs,
    stringsAsFactors = FALSE
  )

  # ----------------------------
  # 3. category assignment
  # ----------------------------
  df$category <- "utils"
  df$category[grepl("^seurat_", df$func)] <- "seurat"
  df$category[grepl("^rna_", df$func)]    <- "rna"
  df$category[grepl("^atac_", df$func)]   <- "atac"
  df$category[grepl("^plot_", df$func)]   <- "plotting"

  # ----------------------------
  # 4. search input
  # ----------------------------
  patterns <- unlist(list(...))
  patterns <- tolower(patterns)

  has_search <- length(patterns) > 0

  df$match <- FALSE
  df$score <- NA_real_

  # ----------------------------
  # 5. EXACT MATCH LOGIC (DEFAULT)
  # ----------------------------
  if (has_search) {

    func_low <- tolower(df$func)

    df$match <- Reduce(`|`, lapply(patterns, function(p) {
      grepl(p, func_low, fixed = TRUE)
    }))

    # ----------------------------
    # fuzzy ranking (OPTIONAL ONLY)
    # ----------------------------
    if (fuzzy) {

      df$score <- sapply(func_low, function(x) {
        min(sapply(patterns, function(p) {
          utils::adist(x, p)
        }))
      })
    }
  }

  # ----------------------------
  # 6. FILTER CORE (CRITICAL FIX)
  # ----------------------------
  if (has_search) {
    df <- df[df$match, , drop = FALSE]
  }

  # ----------------------------
  # 7. optional category filter
  # ----------------------------
  if (!is.null(category)) {
    category <- tolower(category)
    df <- df[df$category == category, , drop = FALSE]
  }

  # ----------------------------
  # 8. ranking (only within matched set)
  # ----------------------------
  df$rank <- 0

  if (fuzzy && has_search && nrow(df) > 0) {
    df$rank <- df$score
  }

  df <- df[order(
    df$rank,
    df$category,
    df$func
  ), ]

  # ----------------------------
  # 9. CLI OUTPUT
  # ----------------------------
  if (!silent) {

    cat("\n==================================================\n")
    cat("        🛠️ DZRLib CLI Function Browser v3.2 FINAL 🛠️\n")
    cat("==================================================\n")

    if (has_search) {
      cat("Mode: SEARCH (exact match only)\n")
      cat("Fuzzy: ", if (fuzzy) "ON" else "OFF", "\n", sep = "")
      cat("Matched:", nrow(df), "\n\n")
    } else {
      cat("Mode: FULL LIST\n")
      cat("Total:", nrow(df), "\n\n")
    }

    # hard stop if nothing matched
    if (has_search && nrow(df) == 0) {
      cat("No matching functions found.\n")
      cat("==================================================\n")
      return(invisible(df))
    }

    for (catg in unique(df$category)) {

      sub <- df[df$category == catg, , drop = FALSE]
      if (nrow(sub) == 0) next

      cat("[", toupper(catg), "]\n", sep = "")

      for (i in seq_len(nrow(sub))) {
        cat("  - ", sub$func[i], "\n", sep = "")
      }

      cat("\n")
    }

    cat("==================================================\n")
  }

  # ----------------------------
  # 10. return invisible table
  # ----------------------------
  invisible(df)
}

