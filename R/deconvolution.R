#' Title
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#' @param lambda Double, regularization parameter for which increasing penalizes
#'    dataset-specific effects
#' @param thresh Double, minimum fractional change in objective function to
#'    continue iteration
#' @param max.iters Integer maximum of iterations to complete before pausing
#' @param nrep Number of random starts to complete
#' @param print.obj Logical, if to print current value of objective
#' @param verbose  Logical, if to print the final objective and best random seed
#'
#' @return nothing
#'
#' @import rliger
#'
#' @export
learn_gene_signatures =function(filepath,
                                analysis.name,
                                spatial.data.name,
                                rand.seed = 123,
                                lambda = 1,
                                thresh = 1e-8,
                                max.iters = 100,
                                nrep = 1,
                                print.obj = FALSE,
                                verbose = FALSE){
  set.seed(rand.seed)

  dir_spatial = file.path(filepath,analysis.name,spatial.data.name,rand.seed)
  norm.data = readRDS(file.path(filepath,analysis.name,rand.seed,"norm_data.RDS"))
  gene_vec = readRDS(file.path(filepath, analysis.name,spatial.data.name,rand.seed,"gene_selection_qc.RDS"))
  gene_vec = intersect(gene_vec, rownames(norm.data[[1]]))
  sample.cells = readRDS(file.path(filepath,analysis.name,rand.seed,"sampled_cells.RDS"))
  clusters = readRDS(file.path(filepath,analysis.name,"source_annotations.RDS"))

  clusters = clusters[sample.cells]

  message("Learning gene signatures")
  E = lapply(norm.data, function(x){t(as.matrix(x)[gene_vec,])})

  if (!all(sapply(X = E, FUN = is.matrix))) {
    stop("All values in 'object' must be a matrix")
  }
  N <- length(x = E)
  ns <- sapply(X = E, FUN = nrow)

  tmp <- gc()

  clust_levels = levels(clusters)
  numeric_clust = as.numeric(clusters)
  names(numeric_clust) = names(clusters)
  H_cells = lapply(norm.data, colnames)
  clust_vec = rep(0, length(clust_levels))
  N = length(E)

  H_indic <- lapply(
    X = H_cells,
    FUN = function(n) {
      numeric_clust_sub = numeric_clust[n]
      mat = do.call(rbind, lapply(X = 1:length(n),
                                  FUN = function(cell){
                                    return(replace(clust_vec,numeric_clust[n[cell]],1))
                                  }))
      rownames(mat) = n
      colnames(mat) = as.character(clust_levels)
      return(mat)
    }
  )
  g <- ncol(x = E[[1]])
  k = length(clust_levels)
  W <- matrix(
    data = abs(x = runif(n = g * k, min = 0, max = 2)),
    nrow = k,
    ncol = g
  )

  V <- lapply(
    X = 1:N,
    FUN = function(i) {
      return(matrix(
        data = abs(x = runif(n = g * k, min = 0, max = 2)),
        nrow = k,
        ncol = g
      ))
    }
  )
  tmp <- gc()
  best_obj <- Inf
  V_m = V
  W_m = W
  H_m = H_indic

  #E = lapply(E, t)
  for (rep_num in 1:nrep) {
    set.seed(seed = rand.seed + rep_num - 1)
    start_time <- Sys.time()
    delta <- 1
    iters <- 0
    pb <- txtProgressBar(min = 0, max = max.iters, style = 3)
    sqrt_lambda <- sqrt(x = lambda)
    obj0 <- sum(sapply(
      X = 1:N,
      FUN = function(i) {
        return(norm(x = E[[i]] - H_indic[[i]] %*% (W + V[[i]]), type = "F") ^ 2)
      }
    )) +
      sum(sapply(
        X = 1:N,
        FUN = function(i) {
          return(lambda * norm(x = H_indic[[i]] %*% V[[i]], type = "F") ^ 2)
        }
      ))
    tmp <- gc()

    while (delta > thresh & iters < max.iters) {
      W <- rliger:::solveNNLS(
        C = do.call(rbind, H_indic),
        B = do.call(rbind, lapply(X = 1:N,
                                  FUN = function(i){
                                    return(E[[i]] - H_indic[[i]] %*% (W + V[[i]]))})
        )
      )

      tmp <- gc()
      V <- lapply(
        X = 1:N,
        FUN = function(i) {
          return(rliger:::solveNNLS(
            C = rbind(H_indic[[i]], sqrt_lambda * H_indic[[i]]),
            B = rbind(E[[i]] - H_indic[[i]] %*% W, matrix(data = 0, nrow = ns[[i]], ncol = g))
          ))
        }
      )
      tmp <- gc()
      obj <- sum(sapply(
        X = 1:N,
        FUN = function(i) {
          return(norm(x = E[[i]] - H_indic[[i]] %*% (W + V[[i]]), type = "F") ^ 2)
        }
      )) +
        sum(sapply(
          X = 1:N,
          FUN = function(i) {
            return(lambda * norm(x = H_indic[[i]] %*% V[[i]], type = "F") ^ 2)
          }
        ))
      tmp <- gc()
      delta <- abs(x = obj0 - obj) / (mean(obj0, obj))
      obj0 <- obj
      iters <- iters + 1
      setTxtProgressBar(pb = pb, value = iters)
    }
    setTxtProgressBar(pb = pb, value = max.iters)

    if (obj < best_obj) {
      W_m <- W
      V_m <- V
      best_obj <- obj
      best_seed <- rand.seed + rep_num - 1
    }

    if (verbose) {
      if (print.obj) {
        cat("Objective:", obj, "\n")
      }
      cat("Best results with seed ", best_seed, ".\n", sep = "")
    }
  }
  out <- list()
  out$H <- H_m
  for (i in 1:length(E)) {
    rownames(x = out$H[[i]]) <- rownames(x = E[[i]])
  }

  out$V <- V_m
  out$W <- W_m

  names(x = out$V) <- names(x = out$H) <- names(x = E)

  rownames(out$W) = clust_levels
  colnames(out$W) = gene_vec

  saveRDS(out, file.path(dir_spatial, "gene_signature_output.RDS"))
}

#' Title
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#' @param cell.size Logical, if to scale gene signatures by cell sizes
#'
#' @return nothing
#'
#' @import rliger
#'
#' @export
deconvolve_spatial = function(filepath,
                              analysis.name,
                              spatial.data.name,
                              rand.seed = 123,
                              cell.size = T){
  set.seed(rand.seed)

  dir_spatial = file.path(filepath,analysis.name,spatial.data.name,rand.seed)
  spatial.data = readRDS(file.path(dir_spatial,"exp_qc.RDS"))
  out = readRDS(file.path(dir_spatial, "gene_signature_output.RDS"))

  W = out[["W"]]
  gene_vec = intersect(colnames(W), rownames(spatial.data))
  spatial.data = spatial.data[gene_vec,]
  W = W[,gene_vec]

  if(cell.size){
    cell_size_mean = readRDS(file.path(filepath,analysis.name,"cell_size.RDS"))
    W = t(sapply(1:nrow(W), function(i){return(cell_size_mean[i]*W[i,])}))
    rownames(W) = names(cell_size_mean)
  }

  message("Deconvolving spatial data")
  spatial.data = t(scale(t(as.matrix(spatial.data)), center = FALSE))
  spatial.data[spatial.data < 0 ] = 0
  spatial.data[is.nan(spatial.data)] = 0
  deconv_h = t(rliger:::solveNNLS(t(W),spatial.data))
  colnames(deconv_h) = rownames(W)
  deconv_frac = t(apply(deconv_h, MARGIN = 1, function(x){x/sum(x)}))
  rownames(deconv_frac) = rownames(deconv_h) = colnames(spatial.data)
  deconv_frac[is.nan(deconv_frac)] = 0
  message("Deconvolution completed")

  saveRDS(list(raw = deconv_h, proportions = deconv_frac), file.path(dir_spatial,"deconvolution_output.RDS"))

  dir_output = file.path(dir_spatial,"downstream_output")
  if(!dir.exists(dir_output)){
    dir.create(paste0(dir_output))
    message("Created directory at ", dir_output)
  }
}
