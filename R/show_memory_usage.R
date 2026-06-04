#' Display Current R Session Memory Usage
#'
#' This function provides a detailed report of the current R session's memory usage,
#' including total memory allocated, memory in use, and memory limits. It's particularly
#' useful for monitoring memory consumption during bioinformatics analyses.
#'
#' @return Returns a data frame with memory statistics and prints a formatted summary.
#' @export
#'
#' @examples
#' # Basic usage
#' show_memory_usage()
#'
#' # Use in a loop to monitor memory during analysis
#' for(i in 1:10) {
#'   # Some memory-intensive operation
#'   show_memory_usage()
#' }
show_memory_usage <- function() {
  # Get memory statistics with error handling
  tryCatch({
    gc_output <- gc(reset = TRUE)
    gc_output_after <- gc()

    # Handle different R versions where gc() output structure may vary
    if ("Vcells" %in% names(gc_output)) {
      total_allocated <- gc_output[["Vcells"]] * object.size(1)
      vectors <- gc_output[["Vcells"]]
    } else {
      # Fallback for older R versions
      total_allocated <- gc_output[["Ncells"]] * mean(sapply(gc_output[1:3], function(x) object.size(x)))
      vectors <- gc_output[["Ncells"]]
    }

    if ("Vcells" %in% names(gc_output_after)) {
      total_in_use <- gc_output_after[["Vcells"]] * object.size(1)
    } else {
      total_in_use <- gc_output_after[["Ncells"]] * mean(sapply(gc_output_after[1:3], function(x) object.size(x)))
    }

    objects <- gc_output_after[["Ncells"]]
    max_memory <- memory.limit() * object.size(1)
    memory_used_pct <- (total_in_use / max_memory) * 100

    mem_stats <- list(
      total_allocated = total_allocated,
      total_in_use = total_in_use,
      vectors = vectors,
      objects = objects,
      max_memory = max_memory,
      memory_used_pct = memory_used_pct
    )

    # Convert to data frame for return value
    mem_df <- as.data.frame(t(mem_stats))
    rownames(mem_df) <- names(mem_stats)

    # Print formatted output
    cat("\n==================================================\n")
    cat("         💾  R Session Memory Usage Report  💾       \n")
    cat("==================================================\n\n")

    cat(sprintf("Total Memory Allocated: %.2f MB\n", mem_stats$total_allocated / 1024^2))
    cat(sprintf("Memory Currently in Use: %.2f MB\n", mem_stats$total_in_use / 1024^2))
    cat(sprintf("Memory Limit: %.2f MB\n", mem_stats$max_memory / 1024^2))
    cat(sprintf("Memory Used: %.1f%%\n", mem_stats$memory_used_pct))
    cat(sprintf("Number of Vectors: %d\n", mem_stats$vectors))
    cat(sprintf("Number of Objects: %d\n", mem_stats$objects))

    if (mem_stats$memory_used_pct > 80) {
      cat("\n⚠️  WARNING: Memory usage is high (>80%)!\n")
      cat("Consider running gc() or breaking your analysis into smaller chunks.\n")
    }

    cat("\n==================================================\n")

    return(mem_df)

  }, error = function(e) {
    cat("\nError accessing memory statistics:\n")
    cat(e$message, "\n")
    cat("Returning basic memory information...\n\n")

    # Fallback to basic memory information
    mem_stats <- list(
      total_allocated = NA,
      total_in_use = NA,
      vectors = NA,
      objects = NA,
      max_memory = memory.limit() * object.size(1),
      memory_used_pct = NA
    )

    mem_df <- as.data.frame(t(mem_stats))
    rownames(mem_df) <- names(mem_stats)

    cat(sprintf("Memory Limit: %.2f MB\n", mem_stats$max_memory / 1024^2))
    cat("Memory statistics could not be determined.\n")

    cat("\n==================================================\n")

    return(mem_df)
  })
}
