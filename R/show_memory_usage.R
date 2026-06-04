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
    # Check if gc() returns a data frame (newer R) or matrix (older R)
    if (is.data.frame(gc_output)) {
      # Newer R version with data frame output
      total_allocated <- gc_output$Ncells * object.size(1)
      vectors <- gc_output$Ncells
      total_in_use <- gc_output_after$Ncells * object.size(1)
      objects <- gc_output_after$Ncells
    } else if (is.matrix(gc_output)) {
      # Older R version with matrix output
      total_allocated <- gc_output[1,1] * object.size(1)
      vectors <- gc_output[1,1]
      total_in_use <- gc_output_after[1,1] * object.size(1)
      objects <- gc_output_after[1,1]
    } else {
      # Fallback for unexpected output
      stop("Unexpected gc() output format")
    }

    # Try memory.limit() but handle cases where it's not available
    max_memory <- tryCatch({
      memory.limit() * object.size(1)
    }, error = function(e) {
      NA_real_
    })

    memory_used_pct <- if (!is.na(max_memory) && max_memory > 0) {
      (total_in_use / max_memory) * 100
    } else {
      NA_real_
    }

    mem_stats <- list(
      total_allocated = total_allocated,
      total_in_use = total_in_use,
      vectors = vectors,
      objects = objects,
      max_memory = max_memory,
      memory_used_pct = memory_used_pct
    )

    # Convert to data frame for return value
    mem_df <- data.frame(
      total_allocated = mem_stats$total_allocated,
      total_in_use = mem_stats$total_in_use,
      vectors = mem_stats$vectors,
      objects = mem_stats$objects,
      max_memory = mem_stats$max_memory,
      memory_used_pct = mem_stats$memory_used_pct
    )

    # Print formatted output
    cat("\n==================================================\n")
    cat("         💾  R Session Memory Usage Report  💾       \n")
    cat("==================================================\n\n")

    cat(sprintf("Total Memory Allocated: %.2f MB\n", mem_stats$total_allocated / 1024^2))
    cat(sprintf("Memory Currently in Use: %.2f MB\n", mem_stats$total_in_use / 1024^2))

    if (!is.na(mem_stats$max_memory)) {
      cat(sprintf("Memory Limit: %.2f MB\n", mem_stats$max_memory / 1024^2))
      cat(sprintf("Memory Used: %.1f%%\n", mem_stats$memory_used_pct))
    } else {
      cat("Memory Limit: Not available\n")
      cat("Memory Used: Not available\n")
    }

    cat(sprintf("Number of Vectors: %d\n", mem_stats$vectors))
    cat(sprintf("Number of Objects: %d\n", mem_stats$objects))

    if (!is.na(mem_stats$memory_used_pct) && mem_stats$memory_used_pct > 80) {
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
      max_memory = NA,
      memory_used_pct = NA
    )

    mem_df <- data.frame(
      total_allocated = mem_stats$total_allocated,
      total_in_use = mem_stats$total_in_use,
      vectors = mem_stats$vectors,
      objects = mem_stats$objects,
      max_memory = mem_stats$max_memory,
      memory_used_pct = mem_stats$memory_used_pct
    )

    cat("Memory statistics could not be determined.\n")
    cat("\n==================================================\n")

    return(mem_df)
  })
}
