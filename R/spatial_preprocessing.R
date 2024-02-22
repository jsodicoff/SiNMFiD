#' Set up new analysis directory
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#'
#' @return nothing
#'
#' @export
start_analysis = function(filepath,
                          analysis.name){
  dir_new = file.path(filepath,analysis.name)
  if(!dir.exists(dir_new)){
    dir.create(dir_new)
    message("Created directory at ", dir_new)
  }
}

#' Add a new spatial dataset to the analysis directory
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.file Path to an RDS file containing desired expression data
#' @param coords.file Path to an RDS file containing desired coordinate data
#' @param spatial.data.name String identifying the spatial sample
#'
#' @return nothing
#' @export
save_spatial_data = function(filepath,
                             analysis.name,
                             spatial.data.file,
                             coords.file,
                             spatial.data.name
){
  dir_new = file.path(filepath,analysis.name,spatial.data.name)
  if(!dir.exists(dir_new)){
    dir.create(dir_new)
    message("Created directory at ", dir_new)
  }
  file.copy(spatial.data.file, file.path(dir_new,"exp.RDS"))
  file.copy(coords.file, file.path(dir_new,"coords.RDS"))
}

#' Generate silouhettes of the data along all three axes
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param save.plots A logical, corresponding with if to save requested plots upon
#'    generation
#'
#' @return nothing
#'
#' @import ggplot2
#' @import cowplot
#' @import viridis
#'
#' @export
reference_3d_coordinates = function(filepath,
                                    analysis.name,
                                    spatial.data.name,
                                    save.plots = FALSE){

  library(ggplot2)
  coords = as.data.frame(readRDS(file.path(filepath,analysis.name,spatial.data.name,"coords.RDS")))

  minmax = apply(coords, MARGIN = 2, function(x){range(x)})
  ranges = apply(minmax, MARGIN = 2, function(y){return(y[2]-y[1])})

  p1 = ggplot(coords, aes_string(x = colnames(coords)[1], y = colnames(coords)[2]))+
    geom_tile() +
    coord_fixed(ratio = 1)+
    geom_tile(aes(alpha = .1)) +
    coord_fixed(ratio = 1) +
    cowplot::theme_minimal_grid() +
    xlab(colnames(coords)[1]) +
    ylab(colnames(coords)[2]) +
    ggtitle(paste0("X-Y silhouette of ",analysis.name))+
    theme(legend.position="none",
          text = ggplot2::element_text(size = 5),
          axis.text = ggplot2::element_text(size = 4),
          plot.title = element_text(size=6))
  p2 = ggplot(coords, aes_string(x = colnames(coords)[1], y = colnames(coords)[3])) +
    geom_tile() +
    coord_fixed(ratio = 1)+
    geom_tile(aes(alpha = .1)) +
    coord_fixed(ratio = 1) +
    cowplot::theme_minimal_grid() +
    xlab(colnames(coords)[1]) +
    ylab(colnames(coords)[3]) +
    ggtitle(paste0("X-Z silhouette of ",analysis.name))+
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
    xlab(colnames(coords)[2]) +
    ylab(colnames(coords)[3]) +
    ggtitle(paste0("Y-Z silhouette of ",analysis.name))+
    theme(legend.position="none",
          text = ggplot2::element_text(size = 5),
          axis.text = ggplot2::element_text(size = 4),
          plot.title = element_text(size=6))

  message("Plots generated")
  print(p1)
  print(p2)
  print(p3)

  if(save.plots){
    plots_dir = file.path(filepath,analysis.name,spatial.data.name,"plots")
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
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param subset.specs A list with length equal to the number of axes, with
#'    each entry a vector of length two, with the first element being the
#'    minimum value to include and the second being the maximum, or NaN to
#'    indicate a missing value
#' @param new.spatial.data.name String, optional name for new analysis, otherwise
#'    the default "`spatial.data.name`_subset_`n.samples`" is used
#' @param out.filepath Path to directory to save subset data to
#'    if not within the analysis
#'
#' @return nothing
#' @export
subset_spatial_data = function(filepath,
                               analysis.name,
                               spatial.data.name,
                               subset.specs = list(c(NaN, NaN),
                                                   c(NaN, NaN),
                                                   c(NaN, NaN)),
                               new.spatial.data.name= NULL,
                               out.filepath = NULL){

  deconv_dir = file.path(filepath,analysis.name,spatial.data.name)
  coords = as.data.frame(readRDS(file.path(filepath,analysis.name,spatial.data.name,"coords.RDS")))
  spatial.data = readRDS(file.path(deconv_dir,"exp.RDS"))

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

  if(is.null(new.spatial.data.name)){
    new.spatial.data.name = paste0(spatial.data.name,"_subset_",nrow(coords))
  }
  if(!is.null(out.filepath)){
    saveRDS(spatial.data, paste0(out.filepath, "/", new.spatial.data.name,"_exp.RDS"))
    saveRDS(coords, paste0(out.filepath, "/", new.spatial.data.name,"_coords.RDS"))
    message("Saved expression and coordinates to ", out.filepath)
  } else {
    new_dir = file.path(filepath,analysis.name,new.spatial.data.name)
    dir.create(new_dir)
    message("Created directory at ", new_dir)
    saveRDS(spatial.data, paste0(new_dir,"/exp.RDS"))
    saveRDS(coords, paste0(new_dir,"/coords.RDS"))
    message("Saved expression and coordinates to ", new_dir)
  }
}

#' Quality-control spatial data
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param count.data Logical, if the spatial data is from a counts
#'     or intensity-based modality
#' @param z Double, the standard deviations above the mean that the number of
#'     NAs in a gene can be before the gene is removed, for intensity data
#' @param n.umi.thresh Integer number of counts below which to remove a sample,
#'     for counts based data
#' @param rand.seed Integer random seed
#'
#' @return nothing
#'
#' @import Matrix
#'
#' @export
qc_spatial_data = function(
    filepath,
    analysis.name,
    spatial.data.name,
    count.data = FALSE,
    z = 1,
    n.umi.thresh = 150,
    rand.seed = 123
){
  set.seed(rand.seed)
  dir_spatial = file.path(filepath,analysis.name,spatial.data.name)
  spatial.data = readRDS(file.path(dir_spatial,"exp.RDS"))
  coords = readRDS(file.path(dir_spatial,"coords.RDS"))
  gene_data = readRDS(file.path(filepath,analysis.name,rand.seed,"gene_selection.RDS"))
  gene_vec = gene_data[[2]]

  if(!count.data){
    spatial.data[spatial.data < 0] = NA
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

  dir_save = file.path(filepath,analysis.name,spatial.data.name,rand.seed)
  dir.create(dir_save)
  message("Created directory at ",dir_save)

  saveRDS(spatial.data, file.path(dir_save,"exp_qc.RDS"))
  saveRDS(coords, file.path(dir_save,"coords_qc.RDS"))
  saveRDS(rownames(spatial.data), file.path(dir_save,"gene_selection_qc.RDS"))
}

#' Coarse-grain spatial data to a predetermined resolution
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param voxel.size Integer, side length of one voxel
#' @param out.filepath Path to directory to save subset data to
#'    if not within the analysis
#' @param verbose Logical, if to print several lines of metadata on results
#'
#' @return nothing
#'
#' @import Matrix
#'
#' @export
voxelize_single_cells = function(
    filepath,
    analysis.name,
    spatial.data.name,
    voxel.size,
    out.filepath = NULL,
    verbose = TRUE
){
  dir_spatial = file.path(filepath,analysis.name,spatial.data.name)
  spatial.data = readRDS(file.path(dir_spatial,"exp.RDS"))
  coords = readRDS(file.path(dir_spatial,"coords.RDS"))

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
    new_dir = file.path(filepath,analysis.name,paste0(spatial.data.name,"_",voxel.size))
    dir.create(new_dir)
    message("Saving to new spatial data folder - " , new_dir)
    out.filepath = new_dir
  }

  if(verbose){
    message(paste0("Generated ", nrow(coords_out), " voxels at ", voxel.size, " cubed resolution." ))
    message(paste0("Mean samples per voxel: ", round(mean(sapply(voxel_list, length)))))
    message(paste0("Mean nUMI per voxel: ", round(mean(colSums(voxel_exp)))))
  }

  saveRDS(voxel_exp, file.path(out.filepath,"exp.RDS"))
  saveRDS(coords_out,file.path(out.filepath,"coords.RDS"))
  saveRDS(voxel_list,file.path(out.filepath,"voxels_to_samples.RDS"))
}

#' Transfer labels from coarse-grained sampled
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param labels.use Named vector of labels for the prevoxelized data
#' @param label.name String identifying the label set
#'
#' @return nothing
#' @export
register_voxel_to_label = function(filepath,
                                   analysis.name,
                                   spatial.data.name,
                                   labels.use,
                                   label.name){

  dir_spatial = file.path(filepath,analysis.name,spatial.data.name)
  voxels_to_samples = readRDS(file.path(dir_spatial,"voxels_to_samples.RDS"))

  unique_labels = unique(labels.use)

  voxel_to_subregion = sapply(unique_labels, function(unique_label){
    label_counts = sapply(names(voxels_to_samples), function(voxel_to_sample){
      sum(labels.use[voxels_to_samples[[voxel_to_sample]]] %in% unique_label)
    })
    return(unique_label[which.max(label_counts)])
  })
  names(voxel_to_subregion) = names(voxels_to_samples)

  saveRDS(proportion_loading_in_subregion, file.path(dir_spatial, paste0("voxel_assignment_by_label_",label.name,".RDS")))
}

#' Flip axes in spatial data
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param axes.flip A vector with three logicals, corresponding to which of
#'    the axes to invert
#' @param overwrite Logical, if the original data should be overwritten,
#'    otherwise "`spatial.data.name`_mirror`_x`/`_y`,`_z`is created
#'
#' @return nothing
#' @export
mirror_spatial_coords = function(filepath,
                                 analysis.name,
                                 spatial.data.name,
                                 axes.flip = c(FALSE,FALSE,FALSE),
                                 overwrite = T){

  deconv_dir = file.path(filepath,analysis.name,spatial.data.name)
  coords = readRDS(file.path(deconv_dir, "coords.RDS"))
  descriptor = "mirror"
  axis_des = c("_x","_y", "_z")

  for(i in 1:3){
    if(axes.flip[i]){
      coords_range = range(coords[,i])
      coords[,i] = -1*(coords[,i]-coords_range[1])+coords_range[2]
      descriptor = paste0(descriptor, axis_des[i])
    }
  }

  if(overwrite){
    saveRDS(coords, paste0(file.path(deconv_dir, "coords.RDS")))
  } else {
    coords_file = tempfile()
    saveRDS(coords, coords_file)
    save_spatial_data(filepath,
                      analysis.name,
                      paste0(file.path(deconv_dir, "exp.RDS")),
                      coords_file,
                      paste0(spatial.data.name,"_",descriptor))
  }

}

#' Use predefined transformations to match some modalities to the Allen CCF
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param ish Logical, if the data comes from the Allen Institute
#'    quantified ISH dataset
#'
#' @return nothing
#' @export
transform_coords_to_ccf = function(
    filepath,
    analysis.name,
    spatial.data.name,
    ish = T){

  deconv_dir = file.path(filepath,analysis.name,spatial.data.name)
  coords = file.path(deconv_dir, "coords.RDS")

  scale_factor = if(ish){200}else{25}
  coords = coords*scale_factor

  if(ish){
    coords[,3] = -coords[,3]+min(coords[,3])+max(coords[,3])+2400
    coords[,2] = -coords[,2]+min(coords[,2])+max(coords[,2])+5500
  } else {
    coords[,2] = -coords[,2]+min(coords[,2])+max(coords[,2])+3000
  }
  coords_file = tempfile()
  saveRDS(coords, coords_file)
  save_spatial_data(filepath,
                    analysis.name,
                    paste0(file.path(deconv_dir, "exp.RDS")),
                    coords_file,
                    paste0(spatial.data.name,"_ccf"))
}
