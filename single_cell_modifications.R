assign_single_cells = function(
  filepath,
  region,
  spatial.data.name,
  rand.seed = 123,
  clusters.from.atlas = TRUE,
  naive.clusters = FALSE){
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  descriptor = as.character(rand.seed)
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
  }
 
  proportions = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))[[2]]
  cell_types = colnames(proportions)
  max = as.factor(apply(proportions, MARGIN = 1, function(x){cell_types[which.max(x)]}))
  names(max) = rownames(proportions)
  saveRDS(max, paste0(dir_spatial,"/deconvolution_max_factor_assignment_",descriptor,".RDS"))
}


#' Convert single cell spatial modalities into voxels by combining samples at
#' a preselected resolution
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param voxel.size A numeric, corresponding to the side-length of desired voxels
#' @param out.filepath A string corresponding to the directory to which to save generated data
#' @param verbose A logical, whether to print details on derived data
#' @return nothing
#'
#' @import
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }

voxelize_single_cells = function(
  filepath,
  region,
  spatial.data.name,
  voxel.size,
  out.filepath = NULL,
  verbose = TRUE
){
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  spatial.data = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_exp.RDS"))
  coords = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_coords.RDS"))


  minmax = apply(coords, MARGIN = 2, function(x){round(range(x),-1)})
  ranges = apply(minmax, MARGIN = 2, function(y){return(y[2]-y[1])})

  coords_adjusted = sapply(1:3, function(i){
    coords[,i] = voxel.size*round(coords[,i]/voxel.size)
  })
  colnames(coords_adjusted) = colnames(coords)

  coords_out = coords_adjusted[!duplicated(coords_adjusted),]
  voxel_exp = matrix(0L, nrow = nrow(spatial.data), ncol = nrow(coords_out))
  voxel_list = list()

  for(i in 1:nrow(coords_out)){
    voxels_use = coords_out[i,1] == coords_adjusted[,1]&
      coords_out[i,2] == coords_adjusted[,2]&
      coords_out[i,3] == coords_adjusted[,3]
    voxel_data = spatial.data[,voxels_use]
    voxel_list[[i]] = colnames(spatial.data)[voxels_use]
    if(!is.null(dim(voxel_data))){
      voxel_exp[,i] = Matrix::rowSums(voxel_data)
    } else {
      voxel_exp[,i] = voxel_data

    }
  }
  rownames(coords_out) = colnames(voxel_exp) = names(voxel_list) = paste0("voxel_",1:nrow(coords_out))
  rownames(voxel_exp) = rownames(spatial.data)


  if(is.null(out.filepath)){
    new_dir = paste0(dir_spatial, "_", voxel.size)
    dir.create(new_dir)
    message("Saving to new spatial data folder - " , new_dir)
    out.filepath = new_dir
  }
  saveRDS(voxel_exp, paste0(out.filepath,"/",spatial.data.name, "_", voxel.size, "_exp.RDS"))
  saveRDS(coords_out,paste0(out.filepath,"/",spatial.data.name, "_", voxel.size, "_coords.RDS"))
  saveRDS(voxel_list,paste0(out.filepath,"/",spatial.data.name, "_", voxel.size, "_voxels_to_samples.RDS"))

  if(verbose){
    message(paste0("Generated ", nrow(coords_out), " voxels at ", voxel.size, " cubed resolution." ))
    message(paste0("Mean samples per voxel: ", round(mean(sapply(voxel_list, length)))))
    message(paste0("Mean nUMI per voxel: ", round(mean(colSums(voxel_exp)))))
  }
}
