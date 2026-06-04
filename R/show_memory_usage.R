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
  # Get memory statistics
  mem_stats <- list(
    total_allocated = gc(reset = TRUE)[["Vcells"]] * object.size(1),
    total_in_use = gc()[["Vcells"]] * object.size(1),
    vectors = gc()[["Vcells"]],
    objects = gc()[["Ncells"]],
    max_memory = memory.limit() * object.size(1),
    memory_used_pct = (gc()[["Vcells"]] / memory.limit()) * 100
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
}
