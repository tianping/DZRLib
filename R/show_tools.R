#' 打印 DZRLib 工具箱的函数分类清单
#'
#' @export
show_tools <- function() {
  all_funcs <- ls("package:DZRLib")
  
  cat("==================================================\n")
  cat("       🛠️  DZRLib 专属生信工具箱清单 🛠️       \n")
  cat("==================================================\n\n")
  
  # 找出不同前缀的函数
  seurat_f <- all_funcs[grep("^seurat_", all_funcs)]
  deseq2_f <- all_funcs[grep("^deseq2_", all_funcs)]
  other_f  <- all_funcs[!all_funcs %in% c(seurat_f, deseq2_f, "show_tools")]
  
  cat("[Seurat helper functions]:\n")
  if(length(seurat_f) > 0) cat(paste("  -", seurat_f, collapse = "\n"), "\n") else cat("  (暂无)\n")
  
  cat("\n[Deseq2 helper functions]:\n")
  if(length(deseq2_f) > 0) cat(paste("  -", deseq2_f, collapse = "\n"), "\n") else cat("  (暂无)\n")
  
  cat("\n[其他通用工具]:\n")
  if(length(other_f) > 0) cat(paste("  -", other_f, collapse = "\n"), "\n") else cat("  (暂无)\n")
  cat("==================================================\n")
}
