#' Set R Console Width
#'
#' This function sets the R console width. If a width is provided, it uses that value.
#' If no width is provided, it attempts to detect the console width using system tools
#' and sets the console width accordingly.
#'
#' @param width Numeric value specifying the console width in characters. If NULL (default),
#'              the function will attempt to detect the console width automatically.
#'
#' @return Returns the console width that was set, invisibly.
#' @export
#'
#' @examples
#' # Set console width to 100 characters
#' set_r_console_width(100)
#'
#' # Auto-detect console width
#' set_r_console_width()
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

# Export the function to make it visible after loading the package
if (getRversion() >= "2.14.0") {
  utils::globalVariables(c("set_r_console_width"))
}
