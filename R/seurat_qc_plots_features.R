#' Generate all QC plots for Seurat object features                                                                                            [112/1808]
#'                                                                                                                                                       
#' @description Convenience function that generates all QC plots (DotPlot, DoHeatmap,                                                                    
#' VlnPlot, and FeaturePlot) for specified features in a Seurat object. This function                                                                    
#' combines the functionality of the individual plot functions.                                                                                          
#'                                                                                                                                                       
#' @param seurat_obj A Seurat object containing the single-cell RNA-seq data.                                                                            
#' @param features Character vector of features (genes) to plot. Optional if feature_files is provided.                                                  
#' @param feature_files Character vector of file paths containing features to plot. Each file should contain one feature per line.                       
#' @param output_dir Directory path where the plots will be saved. Default is current directory.                                                         
#' @param dotplot_file Name of the DotPlot output file (without extension). Default is "QC6.markers_provided.DotPlot".                                   
#' @param heatmap_file Name of the DoHeatmap output file (without extension). Default is "QC6.markers_provided.DoHeatmap".                               
#' @param vlnplot_prefix Prefix for VlnPlot output files. Default is "gene.".                                                                            
#' @param featureplot_prefix Prefix for FeaturePlot output files. Default is "gene.".                                                                    
#' @param ident_var Character string specifying the cluster identifier variable in metadata.                                                             
#' If NULL, uses the current identification in the Seurat object.                                                                                        
#' @param dotplot_params List of parameters to pass to seurat_dotplot_features.                                                                          
#' @param heatmap_params List of parameters to pass to seurat_heatmap_features.                                                                          
#' @param vlnplot_params List of parameters to pass to seurat_violin_feature_plots for VlnPlot.                                                          
#' @param featureplot_params List of parameters to pass to seurat_violin_feature_plots for FeaturePlot.                                                  
#'                                                                                                                                                       
#' @return Invisibly returns a list of paths to all saved plot files.                                                                                    
#' @export                                                                                                                                               
#'                                                                                                                                                       
#' @examples                                                                                                                                             
#' # Basic usage with feature vector                                                                                                                     
#' seurat_qc_plots_features(my.seurat, features = c("PECAM1", "CDH5"))                                                                                   
#'                                                                                                                                                       
#' # With custom output filenames                                                                                                                        
#' seurat_qc_plots_features(my.seurat, features = my_genes,                                                                                              
#'                         dotplot_file = "my_dotplot",                                                                                                  
#'                         heatmap_file = "my_heatmap",                                                                                                  
#'                         vlnplot_prefix = "vln_",                                                                                                      
#'                         featureplot_prefix = "feat_")                                                                                                 
#'
seurat_qc_plots_features <- function(seurat_obj,                                                                                                         
                                    features = NULL,                                                                                                     
                                    feature_files = NULL,                                                                                                
                                    output_dir = "./",                                                                                                   
                                    dotplot_file = "QC6.markers_provided.DotPlot",                                                                       
                                    heatmap_file = "QC6.markers_provided.DoHeatmap",                                                                     
                                    vlnplot_prefix = "gene.",                                                                                            
                                    featureplot_prefix = "gene.",                                                                                        
                                    ident_var = NULL,                                                                                                    
                                    dotplot_params = list(),                                                                                             
                                    heatmap_params = list(),                                                                                             
                                    vlnplot_params = list(),                                                                                             
                                    featureplot_params = list()) {                                                                                       
  # Input validation                                                                                                                                     
  if (!inherits(seurat_obj, "Seurat")) {                                                                                                                 
    stop("seurat_obj must be a Seurat object")                                                                                                           
  }                                                                                                                                                      
                                                                                                                                                         
  if (is.null(features) && is.null(feature_files)) {                                                                                                     
    stop("Either features or feature_files must be provided")                                                                                            
  }                                                                                                                                                      
                                                                                                                                                         
  # Collect all features                                                                                                                                 
  all_features <- character(0)                                                                                                                           

  # Add features from vector                                                                                                                             
  if (!is.null(features)) {                                                                                                                              
    if (!is.character(features)) {                                                                                                                       
      stop("features must be a character vector")                                                                                                        
    }                                                                                                                                                    
    all_features <- c(all_features, features)                                                                                                            
  }                                                                                                                                                      
                                                                                                                                                         
  # Add features from files                                                                                                                              
  if (!is.null(feature_files)) {                                                                                                                         
    if (!is.character(feature_files)) {                                                                                                                  
      stop("feature_files must be a character vector")                                                                                                   
    }                                                                                                                                                    
                                                                                                                                                         
    for (file in feature_files) {                                                                                                                        
      if (file.exists(file)) {                                                                                                                           
        file_features <- readLines(file)                                                                                                                 
        if (length(file_features) > 0) {                                                                                                                 
          all_features <- c(all_features, file_features)                                                                                                 
        }                                                                                                                                                
      } else {                                                                                                                                           
        warning(paste("Feature file not found:", file))                                                                                                  
      }                                                                                                                                                  
    }                                                                                                                                                    
  }                                                                                                                                                      
                                                                                                                                                         
  all_features <- unique(all_features)                                                                                                                   
  if (length(all_features) == 0) {                                                                                                                       
    stop("No valid features found from provided sources")                                                                                                
  }                                                                                                                                                      

  # Generate each plot type                                                                                                                              
  output_paths <- list()                                                                                                                                 
                                                                                                                                                         
  # DotPlot                                                                                                                                              
  if (!is.null(ident_var)) {                                                                                                                             
    dotplot_params$ident_var <- ident_var                                                                                                                
  }                                                                                                                                                      
  output_paths$dotplot <- seurat_dotplot_features(seurat_obj,                                                                                            
                                                 all_features,                                                                                           
                                                 output_dir,                                                                                             
                                                 output_file = dotplot_file,                                                                             
                                                 ... = dotplot_params)                                                                                   
                                                                                                                                                         
  # Heatmap                                                                                                                                              
  output_paths$heatmap <- seurat_heatmap_features(seurat_obj,                                                                                            
                                                  all_features,                                                                                          
                                                  output_dir,                                                                                            
                                                  output_file = heatmap_file,                                                                            
                                                  ... = heatmap_params)                                                                                  
                                                                                                                                                         
  # Violin and Feature plots                                                                                                                             
  output_paths$violin_feature <- seurat_violin_feature_plots(seurat_obj,                                                                                 
                                                             all_features,                                                                               
                                                             output_dir,                                                                                 
                                                             vlnplot_prefix = vlnplot_prefix,                                                            
                                                             featureplot_prefix = featureplot_prefix,                                                    
                                                             vlnplot_params = vlnplot_params,                                                            
                                                             featureplot_params = featureplot_params)                                                    
                                                                                                                                                         
  invisible(output_paths)                                                                                                                                
}                                                                                                                                                        
