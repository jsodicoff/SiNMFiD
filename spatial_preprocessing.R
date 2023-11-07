# Functions required for utilization of the DUMfound deconvolution pipeline on
# fully annotated region atlases

#' @param filepath Path to directory within which the atlas structure was generated
#' @param region The anatomical region with which the analysis is associated
#' @param spatial.data.file A filepath to the spatial expression matrix
#' @param coords An object or filepath to coordinates for the spatial samples
#' @param spatial.data.name A string, the name of the spatial dataset

save_spatial_data = function(filepath,
                        region,
                        spatial.data.file,
                        coords,
                        spatial.data.name
                        ){
  dir_new = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  if(!dir.exists(dir_new)){
    dir.create(dir_new)
    message("Created directory at ", dir_new)
  }
  file.copy(spatial.data.file, paste0(dir_new,"/",spatial.data.name,"_exp.RDS"))
  if(is.character(coords)){
    file.copy(coords, paste0(dir_new,"/",spatial.data.name,"_coords.RDS"))
  } else {
    saveRDS(coords, paste0(dir_new,"/",spatial.data.name,"_coords.RDS"))
  }
}

#' Generate plots of the provided spatial data on the XY, YZ, and XZ planes in
#' a reference grid to allow for indexing.
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param display.plots A logical, corresponding with if to display requested plots upon
#'    generation
#'
#' @return nothing
#'
#' @import ggplot2, cowplot, viridis
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'

reference_3d_coordinates = function(filepath,
                                    region,
                                    spatial.data.name,
                                    save.plots = FALSE){

  library(ggplot2)
  library(cowplot)
  coords = as.data.frame(readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_coords.RDS")))

  minmax = apply(coords, MARGIN = 2, function(x){range(x)})
  ranges = apply(minmax, MARGIN = 2, function(y){return(y[2]-y[1])})

  p1 = ggplot(coords, aes_string(x = colnames(coords)[1], y = colnames(coords)[2]))+
    geom_tile() +
    coord_fixed(ratio = 1)+
    geom_tile(aes(alpha = .1)) +
    coord_fixed(ratio = 1) +
    theme_minimal_grid() +
    #scale_y_continuous(minor_breaks = seq(minmax[1,2] , minmax[2,2], 1), breaks = seq(minmax[1,2] , minmax[2,2], 5)) +
    #scale_x_continuous(minor_breaks = seq(minmax[1,1] , minmax[2,1], 1), breaks = seq(minmax[1,1] , minmax[2,1], 5)) +
    xlab(colnames(coords)[1]) +
    ylab(colnames(coords)[2]) +
    ggtitle(paste0("X-Y silhouette of ",region))+
    theme(legend.position="none",
          text = ggplot2::element_text(size = 5),
          axis.text = ggplot2::element_text(size = 4),
          plot.title = element_text(size=6))
  p2 = ggplot(coords, aes_string(x = colnames(coords)[1], y = colnames(coords)[3])) +
    geom_tile() +
    coord_fixed(ratio = 1)+
    geom_tile(aes(alpha = .1)) +
    coord_fixed(ratio = 1) +
    theme_minimal_grid() +
    #scale_y_continuous(minor_breaks = seq(minmax[1,3] , minmax[2,3], 1), breaks = seq(minmax[1,3] , minmax[2,3], 5)) +
    #scale_x_continuous(minor_breaks = seq(minmax[1,1] , minmax[2,1], 1), breaks = seq(minmax[1,1] , minmax[2,1], 5)) +
    xlab(colnames(coords)[1]) +
    ylab(colnames(coords)[3]) +
    ggtitle(paste0("X-Z silhouette of ",region))+
    theme(legend.position="none",
          text = ggplot2::element_text(size = 5),
          axis.text = ggplot2::element_text(size = 4),
          plot.title = element_text(size=6))
  p3 = ggplot(coords, aes_string(x = colnames(coords)[2], y = colnames(coords)[3])) +
    geom_tile() +
    coord_fixed(ratio = 1)+
    geom_tile(aes(alpha = .1)) +
    coord_fixed(ratio = 1) +
    theme_minimal_grid() +
    #scale_y_continuous(minor_breaks = seq(minmax[1,3] , minmax[2,3], 1), breaks = seq(minmax[1,3] , minmax[2,3], 5)) +
    #scale_x_continuous(minor_breaks = seq(minmax[1,2] , minmax[2,2], 1), breaks = seq(minmax[1,2] , minmax[2,2], 5)) +
    xlab(colnames(coords)[2]) +
    ylab(colnames(coords)[3]) +
    ggtitle(paste0("Y-Z silhouette of ",region))+
    theme(legend.position="none",
          text = ggplot2::element_text(size = 5),
          axis.text = ggplot2::element_text(size = 4),
          plot.title = element_text(size=6))
  print(p1)
  print(p2)
  print(p3)
  message("Plots generated")
  if(save.plots){
    plots_dir = paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots")
    if(!dir.exists(plots_dir)){
      dir.create(plots_dir)
      message("Created directory at ", plots_dir)
    }
    ggplot2::ggsave(paste0(plots_dir, "/x_y_coord_reference.PNG"),
                    p1,
                    width = 500,
                    height = 500,
                    units = "px")
    ggplot2::ggsave(paste0(plots_dir, "/x_z_coord_reference.PNG"),
                    p2,
                    width = 500,
                    height = 500,
                    units = "px")
    ggplot2::ggsave(paste0(plots_dir, "/y_z_coord_reference.PNG"),
                    p3,
                    width = 500,
                    height = 500,
                    units = "px")
    message("Plots saved")
  }
}


#' Subset a spatial dataset by coordinates for analysis
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region The anatomical region with which the analysis is associated
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param subset.specs A list with length equal to the number of axises, with
#'    each entry a vector of length two, with the first element being the
#'    minimum value to include and the second being the maximum, or NaN to
#'    indicate a missing value
#' @param out.filepath Path to directory to save subset data to, if NULL the
#'    expression and coordinates in the atlas directory structure are
#'    overwritten
subset_spatial_data = function(filepath,
                               region,
                               spatial.data.name,
                               subset.specs = list(c(NaN, NaN),
                                                   c(NaN, NaN),
                                                   c(NaN, NaN)),
                               out.filepath = NULL){
  deconv_dir = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/")
  coords = readRDS(paste0(deconv_dir,spatial.data.name,"/",spatial.data.name,"_coords.RDS"))
  spatial.data = readRDS(paste0(deconv_dir,spatial.data.name,"/",spatial.data.name,"_exp.RDS"))
  for(i in 1:ncol(coords)){
    subset.specs[[i]][is.nan(subset.specs[[i]])] = range(coords[,i])[is.nan(subset.specs[[i]])]
    coords = coords[coords[,i] >= subset.specs[[i]][1] & coords[,i] <= subset.specs[[i]][2], ]
  }
  if(nrow(coords) == 0){
    stop("No samples selected -- provide a new range of coordinates.")
  } else {
    message(paste0("Sample subset to ", nrow(coords), " samples."))
  }
  spatial.data = spatial.data[, rownames(coords)]
  new_spatial.data.name = paste0(spatial.data.name,"_subset_",nrow(coords))
  if(!is.null(out.filepath)){
    saveRDS(spatial.data, paste0(out.filepath, "/", new_spatial.data.name,"_exp.RDS"))
    saveRDS(coords, paste0(out.filepath, "/", new_spatial.data.name,"_coords.RDS"))
    message("Saved expression and coordinates to ", out.filepath)
  } else {
    new_dir = paste0(deconv_dir,new_spatial.data.name)
    dir.create(new_dir)
    message("Created directory at ", new_dir)
    saveRDS(spatial.data, paste0(new_dir,"/",new_spatial.data.name,"_exp.RDS"))
    saveRDS(coords, paste0(new_dir,"/",new_spatial.data.name,"_coords.RDS"))
    message("Saved expression and coordinates to ", new_dir)
  }
}

qc_spatial_data = function(
    filepath,
    region,
    spatial.data.name,
    slide.seq = FALSE,
    clusters.from.atlas = TRUE,
    naive.clusters = FALSE,
    z = 1,
    n.umi.thresh = 150,
    rand.seed = 123
){
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    
    dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
    dir_spatial_new = paste0(dir_spatial, "_naive")

    dir.create(dir_spatial_new)
    file.copy(paste0(dir_spatial,"/",spatial.data.name,"_exp.RDS"),dir_spatial_new)
    file.copy(paste0(dir_spatial,"/",spatial.data.name,"_coords.RDS"),dir_spatial_new)
    file.rename(from = paste0(dir_spatial_new,"/",spatial.data.name,"_exp.RDS"), to = paste0(dir_spatial_new,"/",spatial.data.name,"_naive_exp.RDS"))
    file.rename(from = paste0(dir_spatial_new,"/",spatial.data.name,"_coords.RDS"), to = paste0(dir_spatial_new,"/",spatial.data.name,"_naive_coords.RDS"))

    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  
  
  spatial.data = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_exp.RDS"))
  coords = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_coords.RDS"))

  gene_data = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/gene_selection_",descriptor,".RDS"))
  gene_vec = gene_data[[2]]
  
  
  if(!slide.seq){
    spatial.data[spatial.data == -1] = NA
    genes_NA = apply(spatial.data, MARGIN = 1, function(x){sum(is.na(x))})
    mean_genes_NA = mean(genes_NA)
    genes_use = rownames(spatial.data)[genes_NA < (mean_genes_NA + z * mean_genes_NA)]
    message(length(genes_use), " genes saved out of ", nrow(spatial.data), " total.")
    spatial.data = spatial.data[rownames(spatial.data) %in% intersect(genes_use, gene_vec),]
    spatial.data[is.na(spatial.data)] = 0
  } else {
    original_dim = dim(spatial.data)
    spatial.data = spatial.data[rownames(spatial.data) %in% gene_vec, ]
    spatial.data = spatial.data[,Matrix::colSums(spatial.data) > n.umi.thresh]
    message(nrow(spatial.data), " genes used out of ", original_dim[1], " and ", ncol(spatial.data), " cells used out of ", original_dim[2])
  }
  
  coords = coords[colnames(spatial.data),]
  
  saveRDS(spatial.data, paste0(dir_spatial,"/",spatial.data.name,"_exp_qc_",descriptor,".RDS"))
  saveRDS(coords, paste0(dir_spatial,"/",spatial.data.name,"_coords_qc_",descriptor,".RDS"))
  saveRDS(rownames(spatial.data), paste0(dir_spatial,"/gene_selection_qc_",descriptor,".RDS"))
}



