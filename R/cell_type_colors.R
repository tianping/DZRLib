#' Get Standardized Cell Type Colors
#'
#' This function returns standardized colors for different cell types. If specific cell types are requested,
#' it returns only those colors; otherwise, it returns the complete color palette.
#'
#' @param cell_types Optional character vector of cell type names to get colors for. If NULL (default),
#'                   returns colors for all cell types.
#'
#' @return A named vector of colors where names are cell type names and values are hex color codes.
#' @export
#'
#' @examples
#' # Get colors for all cell types
#' get_cell_type_colors()
#'
#' # Get colors for specific cell types
#' get_cell_type_colors(c("Excitatory neurons", "Inhibitory neurons"))
get_cell_type_colors <- function(cell_types = NULL) {
  # Master color palette for all cell types
  all_colors <- c(
    `Patterning centers` = "#b03c64",
    `dorsal NSC` = "#f584e4",
    enIPC = "#7ca4f9",
    `Excitatory neurons` = "#2166ac",
    `CR` = "#bccf42",
    `GE NSC` = "#f1b6da",
    inIPC = "#7fe63e",
    `Inhibitory neurons` = "#0e9c23",
    gIPC = "#ffc277",
    Astro = "#e08214",
    `OPC-Oligo` = "#ad630a",
    Mes = "#6aada3",
    Immune = "#7a7878",
    `RB & Vas` = "#525759",
    `Early subtypes` = "#440a63"
  )

  # If specific cell types requested, return only those colors
  if (!is.null(cell_types)) {
    return(all_colors[cell_types])
  }

  return(all_colors)
}
