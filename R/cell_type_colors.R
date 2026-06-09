# Define a function that returns standardized cell type colors
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
