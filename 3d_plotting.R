#' Generate gifs of cell type distributions derived from deconvolution in space
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param mat.use A string, either "raw", "proportions", or "assignment",
#'    corresponding to the raw cell type loadings, the normalized loadings, or
#'    cell type assignments for single cell spatial modalities
#' @param cell.types.plot A character vector of cell types to plot
#' @param dims A vector of two integers, corresponding the dimensions of the
#'    gifs generated
#'
#' @return nothing
#'
#' @import rgl
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }

generate_loading_gifs = function(
  filepath,
  region,
  spatial.data.name,
  rand.seed = 123,
  clusters.from.atlas = TRUE,
  naive.clusters = FALSE,
  cell.size = FALSE,
  mat.use = "proportions",#raw, proportions, or assignment
  cell.types.plot = NULL,
  filter = NULL,
  dims = c(500, 500)
){
  set.seed(rand.seed)
  library(rgl)
  
  
  
  descriptor = as.character(rand.seed)
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  coords = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_coords_qc_",descriptor,".RDS"))

  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")

  dir_gifs = paste0(dir_output,"/plots")
  if(!dir.exists(dir_gifs)){
    dir.create(paste0(dir_gifs))
    message("Created directory at ", dir_gifs)
  }
  
  if(mat.use != "assignment"){
    loadings = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))
    cell_types = colnames(loadings[[1]])
  } else {
    assignments = readRDS(paste0(dir_spatial,"/deconvolution_max_factor_assignment_",descriptor,".RDS"))
    cell_types = levels(assignments)
  }
  
  if(!is.null(filter) & mat.use != "assignment"){
    loadings[loadings < filter] = 0
    loadings[loadings >= filter] = 1
    descriptor = paste0(descriptor, "_filter_",filter)
  }
  
  
  grDevices::palette(viridis::viridis(option="A",n=50,direction = -1))
  if(is.null(cell.types.plot)){
    cell.types.plot = cell_types
  } else {
    cell.types.plot = intersect(cell.types.plot, cell_types)
  }
  


  for(cell_type in cell.types.plot){
    cell_type_consistent = sub("/", ".",sub(" ", "_",cell_type))
    if(mat.use != "assignment"){
      colors_view = as.numeric(cut(loadings[[mat.use]][,cell_type],breaks=50))
      try(rgl.close(), silent = TRUE)
      open3d(windowRect = c(0,0, dims[1], dims[2]));
      plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, add = TRUE)
      decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
      axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
      movie3d(spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(dir_gifs,"/", region, "_",cell_type_consistent,"_spatial_summary_",descriptor))
    } else {
      colors_view = (assignments == cell_type)*20 + 1
      try(rgl.close(), silent = TRUE)
      open3d(windowRect = c(0,0, dims[1], dims[2]));
      plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, type = "p", add = TRUE)
      decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
      axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
      movie3d(spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(dir_gifs, "/", region, "_",cell_type_consistent,"_spatial_summary_",descriptor))
    }

  }
}

generate_label_gifs = function(
  filepath,
  region,
  spatial.data.name,
  labels.plot,
  dims = c(500, 500)
){
  library(rgl)
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  dir_gifs = paste0(dir_spatial,"/plots")
  if(!dir.exists(dir_gifs)){
    dir.create(paste0(dir_gifs))
    message("Created directory at ", dir_gifs)
  }
  
  grDevices::palette(viridis::viridis(option="A",n=50,direction = -1))
  coords = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_coords.RDS"))
  
  if(is.vector(labels.plot)){
    coords = coords[names(labels.plot),]
    labels_use = unique(labels.plot)
  } else {
    coords = coords[rownames(labels.plot),]
    labels_use = colnames(labels.plot)
  }
  for(label_unique in labels_use){
    if(is.vector(labels.plot)){
      colors_view = (labels.plot == label_unique)*20 + 1
    } else {
      colors_view = labels.plot[,label_unique]*20 + 1
    }
    try(rgl.close(), silent = TRUE)
    open3d(windowRect = c(0,0, dims[1], dims[2]));
    plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, type = "p", add = TRUE)
    decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
    axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
    label_unique = sub("/", ".",sub(" ", "_",label_unique))
    movie3d(spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(dir_gifs, "/", region, "_",label_unique,"_spatial_summary"))
  }
}

overlay_subregion_gifs = function(
    filepath,
    region,
    spatial.data.name,
    rand.seed = 123,
    clusters.from.atlas = TRUE,
    naive.clusters = FALSE,
    cell.size = FALSE,
    mat.use = "proportions",#raw, proportions, or assignment
    cell.types.plot = NULL,
    subregions.plot = NULL,
    filter = NULL,
    dims = c(500, 500)
){
  set.seed(rand.seed)
  library(rgl)
  library(cocoframer)
  
  descriptor = as.character(rand.seed)
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
    
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  coords = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_coords_qc_",descriptor,".RDS"))
  
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")
  
  dir_gifs = paste0(dir_output,"/plots")
  if(!dir.exists(dir_gifs)){
    dir.create(paste0(dir_gifs))
    message("Created directory at ", dir_gifs)
  }
  
  if(mat.use != "assignment"){
    loadings = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))
    cell_types = colnames(loadings[[1]])
  } else {
    assignments = readRDS(paste0(dir_spatial,"/deconvolution_max_factor_assignment_",descriptor,".RDS"))
    cell_types = levels(assignments)
  }
  
  if(!is.null(filter) & mat.use != "assignment"){
    loadings[loadings < filter] = 0
    loadings[loadings >= filter] = 1
    descriptor = paste0(descriptor, "_filter_",filter)
  }
  
  
  grDevices::palette(viridis::viridis(option="A",n=50,direction = -1))
  if(is.null(cell.types.plot)){
    cell.types.plot = cell_types
  } else {
    cell.types.plot = intersect(cell.types.plot, cell_types)
  }  
  
  region_mesh = lapply(subregions.plot, function(x){try(ccf_2017_mesh(x), silent = T)})
  names(region_mesh) = subregions.plot
  region_mesh = region_mesh[sapply(region_mesh, function(x){!any(class(x) == "try-error")})]
  
  for(subregion in names(region_mesh)){
    region_mesh[[subregion]]$material$alpha = .1
    region_mesh[[subregion]]$material$color = "gray"
    
    for(cell_type in cell.types.plot){
      cell_type_consistent = sub("/", ".",sub(" ", "_",cell_type))
      if(mat.use != "assignment"){
        colors_view = as.numeric(cut(loadings[[mat.use]][,cell_type],breaks=50))
        try(rgl.close(), silent = TRUE)
        open3d(windowRect = c(0,0, dims[1], dims[2]));
        plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, add = TRUE)
        shade3d(region_mesh[[subregion]])
        decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
        axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
        movie3d(spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(dir_gifs,"/", region, "_", subregion,"_",cell_type_consistent,"_spatial_summary_",descriptor))
      } else {
        colors_view = (assignments == cell_type)*20 + 1
        try(rgl.close(), silent = TRUE)
        open3d(windowRect = c(0,0, dims[1], dims[2]));
        plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=1, type = "p", add = TRUE)
        shade3d(region_mesh[[subregion]])
        decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
        axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
        movie3d(spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(dir_gifs,"/", region, "_", subregion,"_",cell_type_consistent,"_spatial_summary_",descriptor))
      }
      
    }
  }
}
