#' Get Available CPU Cores
#'
#' This function returns the number of available CPU cores for parallel processing.
#' It provides a cross-platform way to detect available cores, with fallback options
#' if the primary detection method fails.
#'
#' @param fallback Numeric value to use as fallback if core detection fails. Default is 4.
#'
#' @return Numeric value representing the number of available CPU cores.
#' @export
#'
#' @examples
#' # Get available cores
#' cores <- get_available_cores()
#' print(paste("Available cores:", cores))
#'
#' # Use with fallback
#' cores <- get_available_cores(fallback = 2)
get_available_cores <- function(fallback = 4) {
  # Try to detect available cores
  tryCatch({
    # First try parallel::detectCores()
    if (requireNamespace("parallel", quietly = TRUE)) {
      cores <- parallel::detectCores()
      if (cores > 0) {
        return(cores)
      }
    }

    # Fallback to Sys.getenv("NUM_CORES") if set
    env_cores <- as.integer(Sys.getenv("NUM_CORES"))
    if (env_cores > 0) {
      return(env_cores)
    }

    # Final fallback
    return(fallback)

  }, error = function(e) {
    warning(paste("Error detecting cores:", e$message))
    return(fallback)
  })
}
