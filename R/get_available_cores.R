get_available_cores <- function() {
  slurm_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
  if (!is.na(slurm_cores) && slurm_cores > 0) {
    return(slurm_cores)
  }

  # fallback: try future or parallel
  if (requireNamespace("future", quietly = TRUE)) {
    return(future::availableCores())
  } else {
    return(parallel::detectCores())
  }
}
