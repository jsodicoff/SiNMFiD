#' Sample from single cell reference datasets
#'
#' @param objs Either a list of matrices,
#'     a list of paths to RDS files,
#'     or liger objects
#' @param annotations Named vector of cell type assignments by sample
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param datasets.remove, Named list of datasets to remove, if a list of
#'     online liger objects is provided
#' @param n.cells Integer value corresponding to maximum number
#'     of samples per cell type
#' @param rand.seed Integer random seed
#'
#' @return nothing
#'
#' @import rhdf5
#'
#' @export
#' @examples

sample_single_cell = function(
    objs,
    annotations,
    filepath,
    analysis.name,
    datasets.remove = NULL,
    n.cells = 500,
    rand.seed = 123
){

  set.seed(rand.seed)

  message("Loading Data")

  objs = load_objs(objs)
  if(names(objs[[1]])[1] == "data"){
    online_obj = TRUE
  } else {
    online_obj = FALSE
  }

  if(!is.null(datasets.remove) & online_obj){
    datasets_remove = Reduce(c, lapply(names(datasets.remove), function(dataset_name){
      return(paste0(datasets.remove[[dataset_name]],"_",dataset_name))
    }))
  }
  objs = objs[!names(objs) %in% datasets_remove]

  if(online_obj){
    cells_in_sample = lapply(objs, function(i){
      rhdf5::h5read(i$file.path, "/matrix/barcodes")
    })
    gene_in_sample = lapply(objs, function(i){
      rhdf5::h5read(i$file.path, "/matrix/features/name")
    })
  } else {
    cells_in_sample = lapply(objs, colnames)
    gene_in_sample = lapply(objs, rownames)

  }
  reference_cell_names = Reduce(union,cells_in_sample)
  shared_genes = Reduce(intersect,gene_in_sample)


  annotated_cells = intersect(names(annotations),reference_cell_names)
  annotations = annotations[names(annotations) %in% annotated_cells]
  annotations = droplevels(annotations)
  freq_cells = table(annotations)
  freq_cells = freq_cells[names(freq_cells) != ""]


  message("Sampling from clusters")


  sample.cells = Reduce(c, lapply(names(freq_cells), function(cell_type){
    cells = names(annotations[annotations == cell_type])
    if(length(cells)> 0){
      return(sample(cells, min(length(cells), n.cells), replace =FALSE))
    } else {
      return(c())
    }
  }))

  message("Sample summary -- ")
  print(table(annotations[sample.cells]))

  dir_new = file.path(filepath,analysis.name,rand.seed)
  if(!dir.exists(dir_new)){
    dir.create(dir_new)
    message("Created directory at ", dir_new)
  }

  if(online_obj){
    norm.data = lapply(1:length(objs), function(i){
      mat_file = objs[[i]]$file.path
      out_mat = rliger:::Matrix.column_norm(Matrix::sparseMatrix(
        dims = c(length(gene_in_sample[[i]]), length(cells_in_sample[[i]])),
        i = as.numeric(rhdf5::h5read(mat_file, "/matrix/indices")+1),
        p = as.numeric(rhdf5::h5read(mat_file, "/matrix/indptr")),
        x = as.numeric(rhdf5::h5read(mat_file, "/matrix/data"))
      ))
      rownames(out_mat) = gene_in_sample[[i]]
      colnames(out_mat) = cells_in_sample[[i]]
      out_mat = out_mat[gene_in_sample[[i]] %in% shared_genes, cells_in_sample[[i]] %in% sample.cells]
      gene_means = rhdf5::h5read(mat_file, "gene_means")[gene_in_sample[[i]] %in% shared_genes]
      gene_sum_sq = rhdf5::h5read(mat_file, "gene_sum_sq")[gene_in_sample[[i]] %in% shared_genes]
      if(is.vector(out_mat)){
        out_mat = NULL
      } else {
        root_mean_sum_sq = sqrt(gene_sum_sq/(ncol(out_mat)-1))
        out_mat= sweep(out_mat, 1, root_mean_sum_sq, "/")
        out_mat[is.na(out_mat)] = 0
        out_mat[out_mat == Inf] = 0
      }
      return(out_mat)
    })
  } else {
    objs = lapply(1:length(objs), function(i){
      object = objs[[i]]
      object = object[rownames(object) %in% shared_genes,
                          colnames(object) %in% sample.cells]
      rm(object)
      object = object/colSums(object)
      if(is.vector(object)){
        object = NULL
      } else {
        gene_means = rowMeans(object)
        gene_sum_sq = rowSums(object^2)
        root_mean_sum_sq = sqrt(gene_sum_sq/(ncol(object)-1))
        object= sweep(object, 1, root_mean_sum_sq, "/")
        object[is.na(object)] = 0
        object[object == Inf] = 0
      }
      return(object)
    })
  }

  objs = objs[!sapply(objs, function(x){length(x) == 0})]
  saveRDS(objs, file.path(filepath,analysis.name,rand.seed,"norm_data.RDS"))
  saveRDS(sample.cells, file.path(filepath,analysis.name,rand.seed,"sampled_cells.RDS"))
  saveRDS(annotations,file.path(filepath,analysis.name,"source_annotations.RDS"))

}

#' Calculate cell sizes with all reference data
#'
#' @param objs Either a list of matrices,
#'    a list of paths to RDS files,
#'    a list of paths to H5 files,
#'    or liger objects
#' @param annotations Named vector of cell type assignments by sample
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param n.cells Integer value corresponding to maximum number
#'    of samples per cell type
#' @param rand.seed Integer random seed
#' @param plot.hist Logical, if to display and save histograms
#'    of nUMIs by cell type
#'
#' @return nothing
#'
#' @import
#'
#' @export
#' @examples
calculate_cell_sizes = function(
    objs,
    annotations,
    filepath,
    analysis.name,
    plot.hist = F
){

  objs = load_objs(objs)

  size_list = lapply(1:nlevels(annotations),
                     function(x){vector(mode = "integer")})
  names(size_list) = levels(annotations)

  for(object in objs){
    for(cell_type in names(size_list)){
      cells_subset = intersect(colnames(object),
                               names(annotations)[annotations == cell_type])
      if(length(cells_subset)>1){
        size_list[[cell_type]]  = c(size_list[[cell_type]],
                                    Matrix::colSums(object[,cells_subset]))
      } else if(length(cells_subset) == 1){
        size_list[[cell_type]] = c(size_list[[cell_type]],
                                   sum(object[,cells_subset]))
      }
    }
  }
  cell_type_mean = sapply(size_list, mean)

  if(plot.hist){
    pdf(file = file.path(filepath,analysis.name,"cell_size_histogram.pdf"),
        width = 6, height = 4)
    for(i in 1:length(size_list)){
      hist(size_list[[i]],
           main = paste0("Distribution of cell sizes for ",names(size_list)[i]),
           xlab = "Counts")
      abline(v = cell_type_mean[i], col = "red")
    }
    dev.off()
  }

  saveRDS(cell_type_mean, file.path(filepath,analysis.name,"cell_size.RDS"))
}

#' select variable genes with the Kruskal-Wallis test
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param deconv.gene.num Integer, the number of genes to select
#' @param gene.num.tol Integer, the maximum difference between the number of
#' genes selected and `deconv.gene.num`
#' @param rand.seed Integer random seed
#'
#' @return nothing
#'
#' @import matrixTests
#'
#' @export
#' @examples
select_defining_genes = function(
    filepath,
    analysis.name,
    deconv.gene.num = 2000,
    gene.num.tol = 50,
    rand.seed = 123
){

  objs = readRDS(file.path(filepath,analysis.name,rand.seed,"norm_data.RDS"))
  annotations = readRDS(file.path(filepath,analysis.name,"source_annotations.RDS"))
  annotations = as.vector(annotations)

  message("Selecting genes with the KW test")

  chisq_list = list()
  for(i in 1:length(objs)){
    cell_names = colnames(objs[[i]])
    chisq_list[[i]] = matrixTests::row_kruskalwallis(as.matrix(objs[[i]]),
                                                     annotations[names(annotations) %in% cell_names])$statistic
    names(chisq_list[[i]]) = rownames(objs[[i]])
  }
  shared_genes = Reduce(c, lapply(objs, rownames))

  var_thresh_start = 0.5
  var_thresh_old = var_thresh_start
  high = 1
  low = 0

  gene_vec = shared_genes
  for(i in 1:length(chisq_list)){
    chisq_list[[i]][is.na(chisq_list[[i]])] = 0
    gene_vec = intersect(gene_vec, names(chisq_list[[i]][chisq_list[[i]] > quantile(chisq_list[[i]], var_thresh_start)]))
  }

  while(abs(length(gene_vec)-deconv.gene.num) > gene.num.tol){
    if(length(gene_vec)>deconv.gene.num){
      var_thresh_new = (var_thresh_old+high)/2
      low = var_thresh_old
    } else {
      var_thresh_new = (var_thresh_old+low)/2
      high = var_thresh_old
    }
    gene_vec = shared_genes
    for(i in 1:length(chisq_list)){
      gene_vec = intersect(gene_vec,
                           names(chisq_list[[i]][chisq_list[[i]] > quantile(chisq_list[[i]], var_thresh_new, na.rm = TRUE)]))
    }
    var_thresh_old = var_thresh_new
  }

  message(paste0(length(gene_vec), " genes found with p = ",var_thresh_old))

  saveRDS(list(chisq_vals = chisq_list,
               genes_used = gene_vec),
          file.path(filepath,analysis.name,rand.seed,"gene_selection.RDS"))

}

#' Load data from one of multiple formats
#'
#' @param objs Either a list of matrices,
#' a list of paths to RDS files,
#' a list of paths to H5 files,
#' or liger objects
#'
#' @return nothing
#'
#' @import
#'
#' @export
#' @examples
load_objs = function(objs){
  if(!is.list(objs)){
    if(class(objs) == "liger"){
      if(length(objs@h5file.info) != 0){
        return(objs@h5file.info)
      } else {
        return(objs@raw.data)
      }
    } else {
      extension = lapply(objs, function(filepath_obj){
        return(strsplit(filepath_obj,"[.]")[[1]][2])
      })
      if(all(extension %in% "RDS")){
        objs = lapply(objs, function(object_path){
          out = readRDS(object_path)
          class(out) = "numeric"
          return(out)
        })
      }
    }
  } else if(class(objs[[1]]) %in% c("array","matrix")){
    return(objs)
  } else if(class(objs[[1]]) == "liger"){
      file_info = list()
      for(i in 1:length(objs)){
        file_info[[i]] = objs[[i]]@h5file.info
        names(file_info[[i]]) = paste0(names(file_info[[i]]),"_",names(objs)[i])
      }
      objs = Reduce(c, file_info)
      return(objs)
  } else {
    stop("Object input is not currently supported.")
  }
}
