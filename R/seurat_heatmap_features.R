#' Generate DoHeatmap for Seurat object features                                                                                                         
#'                                                                                                                                                       
#' @description Creates a heatmap visualization for specified features in a Seurat object.                                                               
#' The plot is saved to the specified output directory with dynamic height based on                                                                      
#' the number of features.                                                                                                                               
#'                                                                                                                                                       
#' @param seurat_obj A Seurat object containing the single-cell RNA-seq data.                                                                            
#' @param features Character vector of features (genes) to plot.                                                                                         
#' @param output_dir Directory path where the plot will be saved. Default is current directory.                                                          
#' @param output_file Name of the output file (without extension). Default is "QC6.markers_provided.DoHeatmap".                                          
#' @param width_scale Numeric scale factor for plot width. Default is 1 (no scaling).                                                                    
#' @param height_scale Numeric scale factor for plot height. Default is 1 (no scaling).                                                                  
#' @param ... Additional arguments to pass to the DoHeatmap function.                                                                                    
#'                                                                                                                                                       
#' @return Invisibly returns the path to the saved plot file.                                                                                            
#' @export                                                                                                                                               
#'                                                                                                                                                       
#' @examples                                                                                                                                             
#' # Basic usage                                                                                                                                         
#' seurat_heatmap_features(my.seurat, features = c("PECAM1", "CDH5"))                                                                                    
#'                                                                                                                                                       
#' # With custom output filename                                                                                                                         
#' seurat_heatmap_features(my.seurat, features = my_genes,                                                                                               
#'                        output_file = "custom_heatmap")                                                                                                
#'                                                                                                                                                       
seurat_heatmap_features <- function(seurat_obj,                                                                                                          
                                   features,                                                                                                             
                                   output_dir = "./",                                                                                                    
                                   output_file = "QC6.markers_provided.DoHeatmap",                                                                       
                                   width_scale = 1,                                                                                                      
                                   height_scale = 1,                                                                                                     
                                   ...) {                                                                                                                
  # Input validation                                                                                                                                     
  if (!inherits(seurat_obj, "Seurat")) {                                                                                                                 
    stop("seurat_obj must be a Seurat object")                                                                                                           
  }                                                                                                                                                      
                                                                                                                                                         
  if (!is.character(features) || length(features) == 0) {                                                                                                
    stop("features must be a non-empty character vector")                                                                                                
  }                                                                                                                                                      
                                                                                                                                                         
  if (!dir.exists(output_dir)) {                                                                                                                         
    dir.create(output_dir, recursive = TRUE)                                                                                                             
  }                                                                                                                                                      
                                                                                                                                                         
  # Filter features that exist in the object                                                                                                             
  valid_features <- features[features %in% rownames(seurat_obj)]                                                                                         
  if (length(valid_features) == 0) {                                                                                                                     
    warning("No valid features found in the Seurat object")                                                                                              
    return(invisible(NULL))                                                                                                                              
  }                                                                                                                                                      
                                                                                                                                                         
  # Create the heatmap                                                                                                                                   
  do.call(DoHeatmap, c(list(seurat_obj, features = valid_features), list(...)))                                                                          
  # Calculate dynamic height                                                                                                                             
  base_height_per_gene <- 0.2                                                                                                                            
  min_height <- 5                                                                                                                                        
  max_height <- 20                                                                                                                                       
  heatmap_height <- max(min_height,                                                                                                                      
                       min(length(valid_features) * base_height_per_gene,                                                                                
                           max_height))                                                                                                                  
                                                                                                                                                         
  # Save the plot                                                                                                                                        
  output_path <- file.path(output_dir, paste0(output_file, ".pdf"))                                                                                      
  ggplot2::ggsave(output_path,                                                                                                                                    
         device = "pdf",                                                                                                                                 
         width = 16 * width_scale,                                                                                                                       
         height = heatmap_height * height_scale)                                                                                                         
                                                                                                                                                         
  invisible(output_path)                                                                                                                                 
}                                                                                                                                                        
