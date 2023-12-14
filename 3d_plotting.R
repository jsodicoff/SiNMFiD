#' Routine for loading filepaths relevant to 3d plotting functionality
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#'
#' @return
#' @export
#'
#' @examples
three_d_plotting_routine = function(filepath,
                                    analysis.name,
                                    spatial.data.name,
                                    rand.seed){
  set.seed(rand.seed)

  dir_spatial <<- file.path(filepath,analysis.name,spatial.data.name,rand.seed)
  dir_output <<- file.path(dir_spatial,"downstream_output")


  dir_plots <<- file.path(dir_output,"plots")
  if(!dir.exists(dir_plots)){
    dir.create(paste0(dir_plots))
    message("Created directory at ", dir_plots)
  }
  dir_gifs <<- file.path(dir_plots,"gifs")
  if(!dir.exists(dir_gifs)){
    dir.create(paste0(dir_gifs))
    message("Created directory at ", dir_gifs)
  }
}

#' Clear relevant paths out of environment
#'
#' @return
#' @export
#'
#' @examples
three_d_plotting_routine_end = function(){
  rm(dir_spatial, dir_output, dir_plots, dir_gifs)
}

#' Title
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#' @param cell.type A string corresponding to one cell type found in
#'     the deconvolution results
#' @param mat.use A string, either "raw" or "proportions"
#'     referring to what version of the results to summarize
#' @param filter.samples Value for binarizing results, either presence above
#'     the provided threshold or absence below
#' @param dims Integer vector of length 2 corresponding to the width and
#'     height of the RGL window
#'
#' @return
#' @export
#'
#' @examples
view_in_rgl = function(
    filepath,
    analysis.name,
    spatial.data.name,
    rand.seed = 123,
    cell.type,
    mat.use = "proportions",
    filter.samples = NULL,
    dims = c(500, 500)
){
  three_d_plotting_routine(filepath,
                           analysis.name,
                           spatial.data.name,
                           rand.seed)

  coords = readRDS(file.path(dir_spatial,"coords_qc.RDS"))
  loadings = readRDS(file.path(dir_spatial,"deconvolution_output.RDS"))

  if(!is.null(filter.samples)){
    loadings[loadings < filter] = 0
    loadings[loadings >= filter] = 1
  }

  grDevices::palette(viridis::viridis(option="A",n=50,direction = -1))
  colors_view = as.numeric(cut(loadings[[mat.use]][,cell.type],breaks=50))
  try(rgl::close3d(), silent = TRUE)
  rgl::open3d(windowRect = c(0,0, dims[1], dims[2]));
  rgl::plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, add = TRUE)
  rgl::decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
  rgl::axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))

  three_d_plotting_routine_end()
}


#' Generate gifs of cell type distributions derived from deconvolution in space
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#' @param mat.use A string, either "raw" or "proportions"
#'     referring to what version of the results to summarize
#' @param cell.types.plot A character vector of cell types to plot
#' @param dims Integer vector of length 2 corresponding to the width and
#'     height of the RGL window
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
    analysis.name,
    spatial.data.name,
    rand.seed = 123,
    mat.use = "proportions",#raw, proportions
    cell.types.plot = NULL,
    filter = NULL,
    dims = c(500, 500)
){
  three_d_plotting_routine(filepath,
                           analysis.name,
                           spatial.data.name,
                           rand.seed)

  coords = readRDS(file.path(dir_spatial,"coords_qc.RDS"))
  loadings = readRDS(file.path(dir_spatial,"deconvolution_output.RDS"))
  cell_types = colnames(loadings[[1]])

  if(!is.null(filter)){
    loadings[loadings < filter] = 0
    loadings[loadings >= filter] = 1
  }

  if(is.null(cell.types.plot)){
    cell.types.plot = cell_types
  } else {
    cell.types.plot = intersect(cell.types.plot, cell_types)
  }

  grDevices::palette(viridis::viridis(option="A",n=50,direction = -1))
  for(cell_type in cell.types.plot){
    cell_type_consistent = sub("/", ".",sub(" ", "_",cell_type))
    colors_view = as.numeric(cut(loadings[[mat.use]][,cell_type],breaks=50))
    try(rgl::close3d(), silent = TRUE)
    rgl::open3d(windowRect = c(0,0, dims[1], dims[2]));
    rgl::plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, add = TRUE)
    rgl::decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
    rgl::axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
    rgl::movie3d(rgl::spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(dir_gifs,"/",cell_type_consistent,"_spatial_summary_",mat.use))
  }

  three_d_plotting_routine_end()
}

#' Title
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param labels.plot A named vector or matrix of labels to plot for the
#'     provided coordinates
#' @param dims Integer vector of length 2 corresponding to the width and
#'     height of the RGL window
#'
#' @return
#' @export
#'
#' @examples
generate_label_gifs = function(
    filepath,
    analysis.name,
    spatial.data.name,
    labels.plot,
    dims = c(500, 500)
){

  three_d_plotting_routine(filepath,
                           analysis.name,
                           spatial.data.name,
                           rand.seed)

  coords = readRDS(file.path(dir_spatial,"coords_qc.RDS"))


  if(is.vector(labels.plot)){
    coords = coords[names(labels.plot),]
    labels_use = unique(labels.plot)
  } else {
    coords = coords[rownames(labels.plot),]
    labels_use = colnames(labels.plot)
  }

  grDevices::palette(viridis::viridis(option="A",n=50,direction = -1))
  for(label_unique in labels_use){
    if(is.vector(labels.plot)){
      colors_view = (labels.plot == label_unique)*20 + 1
    } else {
      colors_view = labels.plot[,label_unique]*20 + 1
    }
    try(rgl::close3d(), silent = TRUE)
    rgl::open3d(windowRect = c(0,0, dims[1], dims[2]));
    rgl::plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, type = "p", add = TRUE)
    rgl::decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
    rgl::axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
    label_unique = sub("/", ".",sub(" ", "_",label_unique))
    rgl::movie3d(rgl::spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(dir_gifs, "/",label_unique,"_spatial_summary"))
  }

  three_d_plotting_routine_end()
}

#' Title
#'
#' @param filepath
#' @param analysis.name
#' @param spatial.data.name
#' @param rand.seed
#' @param mat.use
#' @param cell.types.plot
#' @param subregions.plot
#' @param filter
#' @param dims
#'
#' @return
#' @export
#'
#' @examples
overlay_subregion_gifs = function(
    filepath,
    analysis.name,
    spatial.data.name,
    rand.seed = 123,
    mat.use = "proportions",#raw, proportions, or assignment
    cell.types.plot = NULL,
    subregions.plot = NULL,
    filter = NULL,
    dims = c(500, 500)
){
  three_d_plotting_routine(filepath,
                           analysis.name,
                           spatial.data.name,
                           rand.seed)

  coords = readRDS(file.path(dir_spatial,"coords_qc.RDS"))
  loadings = readRDS(file.path(dir_spatial,"deconvolution_output.RDS"))
  cell_types = colnames(loadings[[1]])

  if(!is.null(filter)){
    loadings[loadings < filter] = 0
    loadings[loadings >= filter] = 1
  }


  if(is.null(cell.types.plot)){
    cell.types.plot = cell_types
  } else {
    cell.types.plot = intersect(cell.types.plot, cell_types)
  }

  region_mesh = lapply(subregions.plot, function(x){try(cocoframer::ccf_2017_mesh(x), silent = T)})
  names(region_mesh) = subregions.plot
  region_mesh = region_mesh[sapply(region_mesh, function(x){!any(class(x) == "try-error")})]

  grDevices::palette(viridis::viridis(option="A",n=50,direction = -1))
  for(subregion in names(region_mesh)){
    region_mesh[[subregion]]$material$alpha = .1
    region_mesh[[subregion]]$material$color = "gray"

    for(cell_type in cell.types.plot){
      cell_type_consistent = sub("/", ".",sub(" ", "_",cell_type))
      colors_view = as.numeric(cut(loadings[[mat.use]][,cell_type],breaks=50))
      try(close3d(), silent = TRUE)
      rgl::open3d(windowRect = c(0,0, dims[1], dims[2]));
      rgl::plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, add = TRUE)
      rgl::shade3d(region_mesh[[subregion]])
      rgl::decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
      rgl::axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
      rgl::movie3d(rgl::spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(dir_gifs,"/",subregion,"_",cell_type_consistent,"_spatial_summary_",mat.use))

    }
  }

  three_d_plotting_routine_end()
}
