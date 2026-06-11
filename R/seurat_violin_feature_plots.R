#' Generate VlnPlot and FeaturePlot for individual genes                                                                                                 
#'                                                                                                                                                       
#' @description Creates violin plots and feature plots for individual genes in a Seurat object.                                                          
#' Each gene gets its own pair of plots saved to the specified output directory.                                                                         
#'                                                                                                                                                       
#' @param seurat_obj A Seurat object containing the single-cell RNA-seq data.                                                                            
#' @param features Character vector of genes to plot.                                                                                                    
#' @param output_dir Directory path where the plots will be saved. Default is current directory.                                                         
#' @param vlnplot_prefix Prefix for VlnPlot output files. Default is "gene.".                                                                            
#' @param featureplot_prefix Prefix for FeaturePlot output files. Default is "gene.".                                                                    
#' @param vlnplot_params List of parameters to pass to VlnPlot function.                                                                                 
#' @param featureplot_params List of parameters to pass to FeaturePlot function.                                                                         
#' @param width Numeric value for plot width. Default is 7.                                                                                              
#' @param height Numeric value for plot height. Default is 7.                                                                                            
#'                                                                                                                                                       
#' @return Invisibly returns a list of paths to the saved plot files.                                                                                    
#' @export                                                                                                                                               
#'                                                                                                                                                       
#' @examples                                                                                                                                             
#' # Basic usage                                                                                                                                         
#' seurat_violin_feature_plots(my.seurat, features = c("PECAM1", "CDH5"))                                                                                
#'                                                                                                                                                       
#' # With custom prefixes                                                                                                                                
#' seurat_violin_feature_plots(my.seurat, features = my_genes,                                                                                           
#'                            vlnplot_prefix = "vln_",                                                                                                   
#'                            featureplot_prefix = "feat_")                                                                                              
#'                                                                                                                                                       
seurat_violin_feature_plots <- function(seurat_obj,                                                                                                      
                                       features,                                                                                                         
                                       output_dir = "./",                                                                                                
                                       vlnplot_prefix = "gene.",                                                                                         
                                       featureplot_prefix = "gene.",                                                                                     
                                       vlnplot_params = list(),                                                                                          
                                       featureplot_params = list(),                                                                                      
                                       width = 7,                                                                                                        
                                       height = 7) {                                                                                                     
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
  usersMarkers <- features[features %in% rownames(seurat_obj)]                                                                                           
  if (length(usersMarkers) == 0) {                                                                                                                       
    warning("No valid features found in the Seurat object")                                                                                              
    return(invisible(NULL))                                                                                                                              
  }                                                                                                                                                      
                                                                                                                                                         
  output_paths <- list()                                                                                                                                 

                                                                                                                                                         
  for (gene in usersMarkers) {                                                                                                                           
    # Create and save VlnPlot                                                                                                                            
    vln_plot <- do.call(VlnPlot, c(list(seurat_obj, features = gene), vlnplot_params))                                                                   
    vln_path <- file.path(output_dir, paste0(vlnplot_prefix, gene, ".VlnPlot.pdf"))                                                                      
    ggplot2::ggsave(vln_path, width = width, height = height)                                                                                                     
    output_paths$vlnplot <- c(output_paths$vlnplot, vln_path)                                                                                            
                                                                                                                                                         
    # Create and save FeaturePlot                                                                                                                        
    feat_plot <- do.call(FeaturePlot, c(list(seurat_obj, features = gene, order = TRUE), featureplot_params))                                                          
    feat_path <- file.path(output_dir, paste0(featureplot_prefix, gene, ".FeaturePlot.pdf"))                                                             
    ggplot2::ggsave(feat_path, width = width, height = height)                                                                                                    
    output_paths$featureplot <- c(output_paths$featureplot, feat_path)                                                                                   
  }                                                                                                                                                      
                                                                                                                                                         
  invisible(output_paths)                                                                                                                                
}                                                                                                                                                        

