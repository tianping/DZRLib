#' Print a Categorized List of Available Functions in DZRLib
#'
#' @export
show_tools <- function() {
  all_funcs <- ls("package:DZRLib")
  
  # Robust filtering: only keep custom functions containing an underscore
  # This perfectly excludes R system builtins (e.g., system.file, library.dynam)
  user_funcs <- all_funcs[grep("_", all_funcs)]
  
  cat("==================================================\n")
  cat("         🛠️  DZRLib Bioinformatics Toolkit  🛠️       \n")
  cat("==================================================\n")
  cat("Only functions with _ in their names are shown\n\n")
  
  # Group functions by their prefixes
  seurat_f <- user_funcs[grep("^seurat_", user_funcs)]
  deseq2_f <- user_funcs[grep("^deseq2_", user_funcs)]
  other_f  <- user_funcs[!user_funcs %in% c(seurat_f, deseq2_f)]
  
  cat("[Seurat Helper Functions]:\n")
  if(length(seurat_f) > 0) cat(paste("  -", seurat_f, collapse = "\n"), "\n") else cat("  (None)\n")
  
  cat("\n[DESeq2 Helper Functions]:\n")
  if(length(deseq2_f) > 0) cat(paste("  -", deseq2_f, collapse = "\n"), "\n") else cat("  (None)\n")
  
  cat("\n[General Utilities]:\n")
  if(length(other_f) > 0) cat(paste("  -", other_f, collapse = "\n"), "\n") else cat("  (None)\n")
  cat("==================================================\n")
}

