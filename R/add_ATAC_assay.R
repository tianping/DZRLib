## add_ATAC_assay: add an ATAC assay to an existing RNA Seurat object, aligning cell barcodes
# eg. After obtaining the common ATAC peak set, add new ATAC assay to the Seurat object => update ATAC assay
add_ATAC_assay <- function(seurat_obj, atac_assay, assay_name = "ATAC") {
  
  # Find common cells between RNA and ATAC
  common_cells <- intersect(colnames(seurat_obj[["RNA"]]), colnames(atac_assay))
  message(length(common_cells), " cells in common between RNA and ATAC")
  
  # Subset both the Seurat object and the ATAC assay to the common cells
  seurat_obj <- subset(seurat_obj, cells = common_cells)
  atac_assay <- subset(atac_assay, cells = common_cells)
  
  # Check if the cell order is aligned
  if (!all(colnames(seurat_obj[["RNA"]]) == colnames(atac_assay))) {
    stop("Cell order mismatch after subsetting! Something went wrong.")
  }
  
  # Add the ATAC assay to the Seurat object
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj[[assay_name]] <- atac_assay
  
  message("ATAC assay added successfully as '", assay_name, "'")
  
  # Return the updated Seurat object
  return(seurat_obj)
}

