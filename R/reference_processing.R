#' Sample from single cell reference datasets
#'
#' @param data.list Various formats are allowed, including 1. a liger object;
#' 2. a character vector containing file names to RDS/H5 files. 3.
#' Named list of liger object, RDS/H5 file name, matrix/dgCMatrix. List option
#' can have element types mixed. A liger object have to be of version older than
#' 1.99. RDS files must contain base dense matrix or dgCMatrix supported by
#' package "Matrix". H5 files must contain dataset processed by rliger < 1.99.
#' @param annotations Named factor of cell type assignments.
#' @param filepath Path to analysis directory where output sampling needs to be
#' stored.
#' @param analysis.name String identifying the analysis, used to make up a
#' sub-folder name.
#' @param datasets.remove Character vector of datasets to be excluded from
#' sampling if \code{data.list} is a liger object. Named list of dataset names
#' for exluding datasets in liger objects passed with a list \code{data.list}.
#' @param n.cells Integer value corresponding to maximum number of samples per
#' cell type. Default \code{500}.
#' @param rand.seed Integer random seed for reproducible sampling.
#' @param chunk Integer chunk size for processing sparse data stored in H5.
#' Number of cells to load into memory per iteration. Default \code{1000}.
#' @return Nothing is returned. File \code{"norm_data.RDS"} will be stored under
#' \code{"<filepath>/<analysis.name>/<rand.seed>/"}, containing a list of
#' downsampled scaled (not centered) data matrix. File
#' \code{"sampled_cells.RDS"} is stored at the same path, containing barcode
#' vector of the sampled cells. File \code{"source_annotations.RDS"} is stored
#' at \code{"<filepath>/<analysis.name>/"} which contains input
#' \code{annotations}.
#' @export
#' @examples
#' \dontrun{
#' # Explanation for how `datasets.remove` works with example:
#'
#' names(lig@raw.data)
#' # above should show "data1", "data2", "data3", ...
#' # Then when sampling from `lig`, the first two datasets can be excluded with
#' sample_single_cell(data.list = lig, datasets.remove = c("data1", "data2"))
#'
#' # If we got a list of liger object
#' sample_single_cell(data.list = list(human = lig1, mouse = lig2),
#'                    datasets.remove = list(human = c("data1", "data2"),
#'                                           mouse = c("10x1")))
#' }
sample_single_cell = function(
    data.list,
    # objs,
    annotations,
    filepath,
    analysis.name,
    datasets.remove = NULL,
    n.cells = 500,
    rand.seed = 123,
    chunk = 1000
){

  set.seed(rand.seed)

  message("Loading Data")

  objs <- load_objs(data.list, datasets.remove = datasets.remove)
  # if (inherits(objs[[1]], "list")) {
  #   online_obj = TRUE
  # } else {
  #   online_obj = FALSE
  # }
  # datasets_remove <- NULL
  # if (!is.null(datasets.remove) && isH5) {
  #   datasets_remove = Reduce(c, lapply(names(datasets.remove), function(dataset_name){
  #     return(paste0(datasets.remove[[dataset_name]],"_",dataset_name))
  #   }))
  # }
  # objs = objs[!names(objs) %in% datasets_remove]
  cells_in_sample <- lapply(objs, function(x) {
    if (inherits(x, "H5File")) x[["matrix/barcodes"]][]
    else colnames(x)
  })
  gene_in_sample <- lapply(objs, function(x) {
    if (inherits(x, "H5File")) x[["matrix/features/name"]][]
    else rownames(x)
  })
  # if(online_obj){
  #   cells_in_sample = lapply(objs, function(i){
  #     i$barcodes[]
  #     # rhdf5::h5read(i$file.path, "/matrix/barcodes")
  #   })
  #   gene_in_sample = lapply(objs, function(i){
  #     i$genes[]
  #     # rhdf5::h5read(i$file.path, "/matrix/features/name")
  #   })
  # } else {
  #   cells_in_sample = lapply(objs, colnames)
  #   gene_in_sample = lapply(objs, rownames)
  #
  # }
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
  dir_new <- file.path(filepath, analysis.name, rand.seed, fsep = )
  if (!dir.exists(dir_new)) {
    dir.create(dir_new, recursive = TRUE)
    message("Created directory at ", dir_new)
  }
  datasetNames <- names(objs)
  objs <- lapply(seq_along(objs), function(n) {
    x <- objs[[n]]
    barcodes <- cells_in_sample[[n]]
    features <- gene_in_sample[[n]]
    idx <- which(barcodes %in% sample.cells)
    if (inherits(x, "H5File")) {
      message("Reading subsample from H5 dataset: ", names(objs)[n])
      mat <- sampleH5Sparse(
        x, indicesPath = "matrix/indices", indptrPath = "matrix/indptr",
        valuesPath = "matrix/data", nrow = length(features),
        ncol = length(barcodes), idx = idx, chunk = chunk
      )
      mat <- Matrix.column_norm(mat)
      gene_means <- x[["gene_means"]][]
      gene_sum_sq <- x[["gene_sum_sq"]][]
      root_mean_sum_sq <- sqrt(gene_sum_sq/(length(barcodes) - 1))
      root_mean_sum_sq <- root_mean_sum_sq[features %in% shared_genes]
      mat <- mat[features %in% shared_genes, , drop = FALSE]
      # `mat@i + 1` gives the one-based gene index that select the
      # variance for `mat@x`.
      # Don't use sweep() because it coerce sparse to dense
      mat@x <- mat@x / root_mean_sum_sq[mat@i + 1]
      mat[is.na(mat)] <- 0
      mat[mat == Inf] <- 0
      dimnames(mat) <- list(features[features %in% shared_genes], barcodes[idx])
    } else {
      mat <- x[rownames(x) %in% shared_genes, , drop = FALSE]
      mat <- Matrix.column_norm(mat)
      gene_means <- rowMeans(mat)
      gene_sum_sq <- rowSums(mat^2)
      root_mean_sum_sq <- sqrt(gene_sum_sq/(ncol(mat) - 1))
      mat@x <- mat@x / root_mean_sum_sq[mat@i + 1]
      mat[is.na(mat)] <- 0
      mat[mat == Inf] <- 0
      mat <- mat[, idx, drop = FALSE]
    }
    return(mat)
  })
  names(objs) <- datasetNames

  # if(online_obj){
  #   norm.data = lapply(1:length(objs), function(i){
  #     mat_file = objs[[i]]$file.path
  #     mat.i <- objs[[i]]$indices[] + 1
  #     mat.p <- objs[[i]]$indptr[]
  #     mat.x <- objs[[i]]$data[]
  #     out_mat <- Matrix::sparseMatrix(i = mat.i, p = mat.p, x = mat.x,
  #                                     dims = c(length(gene_in_sample[[i]]),
  #                                              length(cells_in_sample[[i]])),
  #                                     dimnames = list(gene_in_sample[[i]],
  #                                                     cells_in_sample[[i]]))
  #     out_mat <- Matrix.column_norm(out_mat)
  #     # out_mat = Matrix.column_norm(Matrix::sparseMatrix(
  #     #   dims = c(length(gene_in_sample[[i]]), length(cells_in_sample[[i]])),
  #     #   i = as.numeric(rhdf5::h5read(mat_file, "/matrix/indices")+1),
  #     #   p = as.numeric(rhdf5::h5read(mat_file, "/matrix/indptr")),
  #     #   x = as.numeric(rhdf5::h5read(mat_file, "/matrix/data"))
  #     # ))
  #     # rownames(out_mat) = gene_in_sample[[i]]
  #     # colnames(out_mat) = cells_in_sample[[i]]
  #     out_mat = out_mat[gene_in_sample[[i]] %in% shared_genes, cells_in_sample[[i]] %in% sample.cells]
  #     gene_means = rhdf5::h5read(mat_file, "gene_means")[gene_in_sample[[i]] %in% shared_genes]
  #     gene_sum_sq = rhdf5::h5read(mat_file, "gene_sum_sq")[gene_in_sample[[i]] %in% shared_genes]
  #     if(is.vector(out_mat)){
  #       # possible condition branch?
  #       out_mat = NULL
  #     } else {
  #       root_mean_sum_sq = sqrt(gene_sum_sq/(ncol(out_mat)-1))
  #       out_mat= sweep(out_mat, 1, root_mean_sum_sq, "/")
  #       out_mat[is.na(out_mat)] = 0
  #       out_mat[out_mat == Inf] = 0
  #     }
  #     return(out_mat)
  #   })
  # } else {
  #   objs = lapply(1:length(objs), function(i){
  #     object = objs[[i]]
  #     object = object[rownames(object) %in% shared_genes,
  #                         colnames(object) %in% sample.cells]
  #     rm(object)
  #     object = object/colSums(object)
  #     if(is.vector(object)){
  #       object = NULL
  #     } else {
  #       gene_means = rowMeans(object)
  #       gene_sum_sq = rowSums(object^2)
  #       root_mean_sum_sq = sqrt(gene_sum_sq/(ncol(object)-1))
  #       object= sweep(object, 1, root_mean_sum_sq, "/")
  #       object[is.na(object)] = 0
  #       object[object == Inf] = 0
  #     }
  #     return(object)
  #   })
  # }

  objs <- objs[sapply(objs, function(x) ncol(x) > 0)]
  saveRDS(objs, file.path(filepath, analysis.name, rand.seed, "norm_data.RDS"))
  saveRDS(sample.cells, file.path(filepath, analysis.name, rand.seed, "sampled_cells.RDS"))
  saveRDS(annotations, file.path(filepath, analysis.name, "source_annotations.RDS"))
}

Matrix.column_norm <- function(A) {
  if (inherits(A, "dgTMatrix")) {
    temp = summary(A)
    A = Matrix::sparseMatrix(i = temp[, 1], j = temp[, 2], x = temp[, 3])
  }
  A@x <- A@x/rep.int(Matrix::colSums(A), diff(A@p))
  return(A)
}


#' Calculate cell sizes with all reference data
#'
#' @param data.list Various formats are allowed, including 1. a liger object;
#' 2. a character vector containing file names to RDS/H5 files. 3.
#' Named list of liger object, RDS/H5 file name, matrix/dgCMatrix. List option
#' can have element types mixed. A liger object have to be of version older than
#' 1.99. RDS files must contain base dense matrix or dgCMatrix supported by
#' package "Matrix". H5 files must contain dataset processed by rliger < 1.99.
#' @param annotations Named factor of all cell type assignments, should be
#' concatenated from all datasets.
#' @param filepath Path to analysis directory where output sampling needs to be
#' stored.
#' @param analysis.name String identifying the analysis, used to make up a
#' sub-folder name.
#' @param datasets.remove Character vector of datasets to be excluded from
#' sampling if \code{data.list} is a liger object. Named list of dataset names
#' for exluding datasets in liger objects passed with a list \code{data.list}.
#' See \code{\link{sample_single_cell}} examples.
#' @param plot.hist Logical, if to display and save histograms of nUMIs by cell
#' type
#' @param chunk Integer chunk size for processing sparse data stored in H5.
#' Number of cells to load into memory per iteration. Default \code{1000}.
#' @return Nothing is returned, but the following file will be stored to local:
#' \itemize{
#' \item{\code{"<filepath>/<analysis.name>/cell_size_histogram.pdf"} - A PDF
#' file for the histogram that shows nUMI per cell distribution for each cell
#' type}
#' \item{\code{"<filepath>/<analysis.name>/cell_size.RDS"} - RDS file of a
#' named numeric vector object, total number of counts per cell type across all
#' datasets.}
#' }
#' @export
calculate_cell_sizes = function(
    data.list,
    annotations,
    filepath,
    analysis.name,
    datasets.remove = NULL,
    plot.hist = FALSE,
    chunk = 1000
){
  objs <- load_objs(data.list, datasets.remove = datasets.remove)
  # objs = load_objs(objs)

  size_list = lapply(levels(annotations),
                     function(x) integer(0))
  names(size_list) = levels(annotations)

  for(object in objs){
    if (inherits(object, "dgCMatrix")) {
      isec <- intersect(colnames(object), names(annotations))
      object <- object[, isec, drop = FALSE]
      ann <- annotations[isec, drop = TRUE]
      for (ct in levels(ann)) {
        if (sum(ann == ct) == 0) next
        size_list[[ct]] <- c(
          size_list[[ct]],
          unname(colSums(object[, ann == ct, drop = FALSE]))
        )
      }
    } else if (inherits(object, "H5File")) {
      # Assuming we use 10X data structure
      # i.e. i - matrix/indices; p - matrix/indptr; x - matrix/data
      #      colnames - matrix/barcodes; rownames - matrix/features/name
      i.link <- object[["matrix/indices"]]
      p.link <- object[["matrix/indptr"]]
      x.link <- object[["matrix/data"]]
      barcodes <- object[["matrix/barcodes"]][]
      ncol <- length(barcodes)
      nrow <- object[["matrix/features/name"]]$dims
      nchunk <- ceiling(ncol / chunk)
      p.start <- NULL
      p.end <- 1
      mat.out <- NULL
      message("Counting for HDF5 data at: ", object$filename)
      pb <- utils::txtProgressBar(min = 0, max = nchunk, style = 3, file = stderr())
      for (i in seq_len(nchunk)) {
        p.start <- p.end
        p.end <- p.start + chunk
        if (p.end > (ncol + 1)) p.end <- ncol + 1
        p.orig <- p.link[p.start:p.end]
        p.new <- p.orig - p.orig[1]
        i.new <- i.link[(p.orig[1] + 1):p.orig[length(p.orig)]]
        x.new <- x.link[(p.orig[1] + 1):p.orig[length(p.orig)]]
        mat.chunk <- Matrix::sparseMatrix(
          i = i.new + 1, p = p.new, x = x.new, dims = c(nrow, p.end - p.start)
        )
        # idx.orig - The 1-based integer index of cells of the current chunk
        #            out of the whole dataset
        idx.orig <- seq(p.start, p.end - 1)
        isec <- intersect(barcodes[idx.orig], names(annotations))
        mat.chunk <- mat.chunk[, isec, drop = FALSE]
        ann <- annotations[isec, drop = TRUE]
        for (ct in levels(ann)) {
          if (sum(ann == ct) == 0) next
          size_list[[ct]] <- c(
            size_list[[ct]],
            unname(colSums(mat.chunk[, ann == ct, drop = FALSE]))
          )
        }
        utils::setTxtProgressBar(pb, i)
      }
      message()
    }
    # for(cell_type in names(size_list)){
    #   cells_subset = intersect(colnames(object),
    #                            names(annotations)[annotations == cell_type])
    #   if(length(cells_subset)>1){
    #     size_list[[cell_type]]  = c(size_list[[cell_type]],
    #                                 Matrix::colSums(object[,cells_subset]))
    #   } else if(length(cells_subset) == 1){
    #     size_list[[cell_type]] = c(size_list[[cell_type]],
    #                                sum(object[,cells_subset]))
    #   }
    # }
  }
  cell_type_mean = sapply(size_list, mean)

  if (isTRUE(plot.hist)) {
    pdf(file = file.path(filepath, analysis.name, "cell_size_histogram.pdf"),
        width = 6, height = 4)
    for (i in seq_along(size_list)) {
      hist(size_list[[i]],
           main = paste0("Distribution of cell sizes for ", names(size_list)[i]),
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
#' @param objs A named list of matrices (dgCMatrix), RDS file paths to matirces,
#' H5 file paths to LIGER analyzed datasets.
#' @return list object. List element type depends on input.
#' @export
load_objs <- function(objs, datasets.remove) {
  if (inherits(objs, "liger")) {
    # The `raw.data` slot stores what we need regardless of
    # whether it's in H5 or in memory.
    # If H5, it should be a list of `H5File` class object, which is a link to the
    # H5 file constructed by package hdf5r
    # If in memory, it is a list of `dgCMatrix`
    return(objs@raw.data[!names(objs@raw.data) %in% datasets.remove])
  }
  if (is.character(objs)) {
    # File name(S) in character vector
    data.list <- lapply(objs, function(path) {
      x.upper <- toupper(path)
      if (endsWith(x.upper, "RDS")) return(readRDS(path))
      if (endsWith(x.upper, "H5")) return(hdf5r::H5File$new(path, mode = "r"))
      else warning("Extension name not supported: ", path)
    })
    return(data.list)
  }
  # Other scenario requires a named list
  if (!is.list(objs) || is.null(names(objs))) {
    stop("`objs` has to be one of: \n",
         "1. a liger object (< 1.99)\n",
         "2. a character vector of RDS or H5 file paths\n",
         "3. a named list of matrix, dgCMatrix, string of RDS/H5 file path, liger object (< 1.99)")
  }
  data.list <- lapply(names(objs), function(n) {
    x <- objs[[n]]
    wrap <- list()
    if (inherits(x, "matrix") || inherits(x, "dgCMatrix")) {
      wrap[[n]] <- x
    } else if (is.character(x)) {
      x.upper <- toupper(x)
      if (endsWith(x.upper, "RDS")) wrap[[n]] <- readRDS(x)
      if (endsWith(x.upper, "H5")) wrap[[n]] <- hdf5r::H5File$new(x, mode = "r")
      else warning("Extension name not supported: ", x)
    }
    else if (inherits(x, "liger")) {
      keep <- !names(x@raw.data) %in% datasets.remove[[n]]
      wrap <- x@raw.data[keep]
      if (length(wrap) > 0) names(wrap) <- paste0(names(wrap), '_', n)
    }
    else warning("List element of class `", class(x)[1], "` is not supported.")
    return(wrap)
  })
  data.list <- Reduce(c, data.list)
  return(data.list)
}
# load_objs = function(objs){
#   if(!is.list(objs)){
#     if(class(objs) == "liger"){
#       if(length(objs@h5file.info) != 0){
#         return(objs@h5file.info)
#       } else {
#         return(objs@raw.data)
#       }
#     } else {
#       extension = lapply(objs, function(filepath_obj){
#         return(strsplit(filepath_obj,"[.]")[[1]][2])
#       })
#       if(all(extension %in% "RDS")){
#         objs = lapply(objs, function(object_path){
#           out = readRDS(object_path)
#           class(out) = "numeric"
#           return(out)
#         })
#       }
#     }
#   } else if(class(objs[[1]]) %in% c("array","matrix")){
#     return(objs)
#   } else if(class(objs[[1]]) == "liger"){
#     file_info = list()
#     for(i in 1:length(objs)){
#       file_info[[i]] = objs[[i]]@h5file.info
#       names(file_info[[i]]) = paste0(names(file_info[[i]]),"_",names(objs)[i])
#     }
#     objs = Reduce(c, file_info)
#     return(objs)
#   } else {
#     stop("Object input is not currently supported.")
#   }
# }

# Helper function that directly reads cell-subset H5 based CCS sparse matrix
# without loading the whole data into memory
# --------
# h5file - H5File class, created with package hdf5r
# xxPath - The path in the H5 file that points to the CCS sparse matrix
#          constructor vectors
# nrow   - Number of rows of the sparse matrix
# idx    - Integer subscriber, 1-based, which columns need to be extracted
# --------
# Returns column subset of original data in `dgCMatrix` class. Columns comes
# with exact order as `idx` specifies. Faster when `idx` asks for contiguous
# positions due to how HDF5 is designed.
sampleH5Sparse <- function(
    h5file,
    indicesPath,
    indptrPath,
    valuesPath,
    nrow,
    ncol,
    idx,
    chunk
) {
  i.link <- h5file[[indicesPath]]
  p.link <- h5file[[indptrPath]]
  x.link <- h5file[[valuesPath]]

  nchunk <- ceiling(ncol / chunk)
  p.start <- NULL
  p.end <- 1
  mat.out <- NULL
  pb <- utils::txtProgressBar(min = 0, max = nchunk, style = 3, file = stderr())
  for (i in seq_len(nchunk)) {
    p.start <- p.end
    p.end <- p.start + chunk
    if (p.end > (ncol + 1)) p.end <- ncol + 1
    idx.orig <- seq(p.start, p.end - 1)
    if (any(idx.orig %in% idx)) {
      p.orig <- p.link[p.start:p.end]
      p.new <- p.orig - p.orig[1]
      i.new <- i.link[(p.orig[1] + 1):p.orig[length(p.orig)]]
      x.new <- x.link[(p.orig[1] + 1):p.orig[length(p.orig)]]
      mat.chunk <- Matrix::sparseMatrix(
        i = i.new + 1, p = p.new, x = x.new, dims = c(nrow, p.end - p.start)
      )

      mat.chunk.sub <- mat.chunk[, idx.orig %in% idx, drop = FALSE]
      if (is.null(mat.out)) mat.out <- mat.chunk.sub
      else mat.out <- cbind(mat.out, mat.chunk.sub)
    }
    utils::setTxtProgressBar(pb, i)
  }
  message()
  return(mat.out)
}
