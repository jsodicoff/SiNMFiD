learn_gene_signatures =function(filepath,
                              region,
                              spatial.data.name,
                              lambda = 1,
                              thresh = 1e-8,
                              max.iters = 100,
                              nrep = 1,
                              rand.seed = 123,
                              print.obj = FALSE,
                              clusters.from.atlas = TRUE,
                              naive.clusters = FALSE,
                              verbose = FALSE){
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  norm.data = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",descriptor,"_norm_data.RDS"))
  gene_vec = readRDS(paste0(dir_spatial,"/gene_selection_qc_",descriptor,".RDS"))
  gene_vec = intersect(gene_vec, rownames(norm.data[[1]]))
  
  sample.cells = readRDS(paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/sampled_cells_", descriptor,".RDS")))
  
  if(clusters.from.atlas){
     clusters = readRDS(paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/clusters_",descriptor,".RDS")))
  } else {
     clusters = readRDS(paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/user_defined_clusters_",descriptor,".RDS")))
  }
  clusters = clusters[sample.cells]
  
  message("Learning gene signatures")
  E = lapply(norm.data, function(x){t(as.matrix(x)[gene_vec,])})
  
  if (!all(sapply(X = E, FUN = is.matrix))) {
    stop("All values in 'object' must be a matrix")
  }
  N <- length(x = E)
  ns <- sapply(X = E, FUN = nrow)
  #if (k >= min(ns)) {
  #  stop('Select k lower than the number of cells in smallest dataset: ', min(ns))
  #}
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
    # if (iters == max.iters) {
    #   print("Warning: failed to converge within the allowed number of iterations.
    #         Re-running with a higher max.iters is recommended.")
    # }
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
  
  out$H_refined <- lapply(X = 1:length(E),
                          FUN = function(i){
                            mat = t(x = rliger:::solveNNLS(
                              C = rbind(t(x = W) + t(x = V[[i]]), sqrt_lambda * t(x = V[[i]])),
                              B = rbind(t(x = E[[i]]), matrix(data = 0, nrow = g, ncol = ns[[i]]))
                            )
                            )
                            rownames(mat) = rownames(E[[i]])
                            return(mat)
                          })
  
  names(x = out$V) <- names(x = out$H) <- names(x = out$H_refined) <- names(x = E)
  
  rownames(out$W) = clust_levels
  colnames(out$W) = gene_vec
  
  saveRDS(out, paste0(dir_spatial, "/gene_signature_output_",descriptor,".RDS"))
}

deconvolve_spatial = function(filepath,
                              region,
                              spatial.data.name,
                              rand.seed = 123,
                              clusters.from.atlas = TRUE,
                              naive.clusters = FALSE,
                              cell.size = F,
                              W = NULL){
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  spatial.data = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_exp_qc_",descriptor,".RDS"))
  
  if(is.null(W)){
    out = readRDS(paste0(dir_spatial, "/gene_signature_output_",descriptor,".RDS"))
    W = out[["W"]]
  } else {
    W = readRDS(W)
    descriptor = paste0(rand.seed, "_custom_w")
  }
  
  gene_vec = intersect(colnames(W), rownames(spatial.data))
  spatial.data = spatial.data[gene_vec,]
  W = W[,gene_vec]
  
  if(cell.size){
    if(naive.clusters){
      cell_size_mean = readRDS(paste0(filepath,"/",region, "/", region, "_Deconvolution_Output/cell_size_",rand.seed,"_object_clusters_naive.RDS"))
    } else {
      cell_size_mean = readRDS(paste0(filepath, "/", region, "/", region, "_Deconvolution_Output/cell_size_",rand.seed,"_object_clusters.RDS"))
    }
    W = t(sapply(1:nrow(W), function(i){return(cell_size_mean[i]*W[i,])}))
    rownames(W) = names(cell_size_mean)
    descriptor = paste0(descriptor, "_size_scaled")
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
  saveRDS(list(raw = deconv_h, proportions = deconv_frac), paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))
  message("Deconvolution completed")
  
  dir.create(paste0(dir_spatial,"/",descriptor,"_output"))
  
}
