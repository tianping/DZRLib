set_r_console_width <- function(width = NULL) {
  # If width is provided, use it; otherwise detect console width
  if (!is.null(width)) {
    if (!is.numeric(width) || length(width) != 1 || width <= 0) {
      stop("Width must be a positive number")
    }
    console_width <- as.integer(width)
  } else {
    # Try to detect console width using tput
    tryCatch({
      console_width <- as.integer(system2("tput", "cols", stdout = TRUE))
    }, error = function(e) {
      # Fallback to default width if detection fails
      console_width <- 80
      warning("Could not detect console width: ", e$message)
    })
  }

  # Set the R console width
  options(width = console_width)
  invisible(console_width)
}
