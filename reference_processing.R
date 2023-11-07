sample_single_cell = function(
    filepath,
    region,
    n.cells = 500,
    known.annotations = NULL,
    naive.clusters = FALSE,
    naive.clusters.remove = NULL,
    rand.seed = 123
){
  
  set.seed(rand.seed)
  
  message("Loading Data")
  
  object_paths = paste0(filepath,"/", region, "/Analysis",c(2,4,5),"_", region, "/onlineINMF_",region, "_object.RDS" )
  objects = lapply(object_paths, function(object_path){readRDS(object_path)})
  
  h5_files = Reduce(c, lapply(objects, function(object){sapply(1:length(object@norm.data), function(i){object@h5file.info[[i]]$file.path})}))
  rna_files = grep(paste0("_(sc10Xv3_|smartseq_|sn10Xv3_|sc10Xv2_)"), h5_files, value = TRUE)
  
  liger_cells = lapply(rna_files, function(i){
    rhdf5::h5read(i, "/matrix/barcodes")#change, extract from H5
  })
  
  liger_genes = lapply(rna_files, function(i){
    rhdf5::h5read(i, "/matrix/features")[[1]] #change, extract from H5
  })
  
  
  
  
  if(is.null(known.annotations)){
    clusters = c()
    for(analysis_num in c(2,4,5)){
      analysis_results = readRDS(paste0(filepath,"/",  region, "/Analysis", analysis_num, "_", region, "/Analysis", analysis_num, "_", region,"_Results_Table.RDS"))
      if("highRcluster" %in% colnames(analysis_results)){
        clust.use = "highRcluster"
      } else {
        clust.use = "lowRcluster"
      }
      if(naive.clusters){
        analysis_clusters = paste0(analysis_num ,"_",analysis_results[,clust.use])
        analysis_clusters[analysis_clusters %in% paste0(analysis_num, "_", naive.clusters.remove[[as.character(analysis_num)]])] = ""
      } else {
        analysis_clusters = as.character(analysis_results$highRAnnotations)
      }
      names(analysis_clusters) = analysis_results$Barcode
      clusters = c(clusters, analysis_clusters)
    }
    clusters = clusters[clusters != ""]
    clusters = as.factor(clusters)
  } else {
    clusters = known.annotations
  }
  

  annotated_cells = intersect(names(clusters),Reduce(union, liger_cells))
  
  clusters = clusters[names(clusters) %in% annotated_cells]
  clusters = droplevels(clusters)
  
  
  
  freq_cells = table(clusters)
  freq_cells = freq_cells[names(freq_cells) != ""]
  
  
  message("Sampling from clusters")
  
  liger_cells_combined = Reduce(c, liger_cells)
  
  sample.cells = Reduce(c, lapply(names(freq_cells), function(cell_type){
    cells = intersect(names(clusters[clusters == cell_type]), liger_cells_combined)
    if(length(cells)> 0){
      return(sample(cells, min(length(cells), n.cells), replace =FALSE))
    } else {
      return(c())
    }
  }))
  
  message("Sample summary -- ")
  print(table(clusters[sample.cells]))
  
  descriptor = as.character(rand.seed)
  
  if(is.null(known.annotations)){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
  }
  
  saveRDS(sample.cells, paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/sampled_cells_", descriptor,".RDS")))
  
  if(is.null(known.annotations)){
     saveRDS(clusters, paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/clusters_",descriptor,".RDS")))
  } else {
     saveRDS(clusters, paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/user_defined_clusters_",descriptor,".RDS")))
  }
          
  
  shared_genes = Reduce(intersect, liger_genes)

  norm.data = lapply(1:length(rna_files), function(i){
    n = rna_files[i]
    out_mat = rliger:::Matrix.column_norm(Matrix::sparseMatrix(
      dims = c(length(liger_genes[[i]]), length(liger_cells[[i]])),
      i = as.numeric(rhdf5::h5read(n, "/matrix/indices")+1),
      p = as.numeric(rhdf5::h5read(n, "/matrix/indptr")),
      x = as.numeric(rhdf5::h5read(n, "/matrix/data"))
    ))
    rownames(out_mat) = liger_genes[[i]]
    colnames(out_mat) = liger_cells[[i]]
    out_mat = out_mat[liger_genes[[i]] %in% shared_genes, liger_cells[[i]] %in% sample.cells]
    gene_means = rhdf5::h5read(rna_files[i], "gene_means")[liger_genes[[i]] %in% shared_genes]
    gene_sum_sq = rhdf5::h5read(rna_files[i], "gene_sum_sq")[liger_genes[[i]] %in% shared_genes]
    if(is.vector(out_mat)){
      out_mat = NULL
    } else {
      root_mean_sum_sq = sqrt(gene_sum_sq/(ncol(out_mat)-1))
      out_mat= sweep(out_mat, 1, root_mean_sum_sq, "/") #liger_cells[[i]] %in% sample.cells
      out_mat[is.na(out_mat)] = 0
      out_mat[out_mat == Inf] = 0
    }
    return(out_mat)
  })
  norm.data = norm.data[!sapply(norm.data, function(x){length(x) == 0})]
  saveRDS(norm.data, paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",descriptor,"_norm_data.RDS"))
}

calculate_cell_sizes = function(
    filepath,
    region,
    naive.clusters = FALSE,
    rand.seed = 123,
    plot.hist = F
    ){
      
    object_paths = paste0(filepath,"/", region, "/Analysis",c(2,4,5),"_", region, "/onlineINMF_",region, "_object.RDS" )
    objects = lapply(object_paths, function(object_path){readRDS(object_path)})
    
    h5_files = Reduce(c, lapply(objects, function(object){sapply(1:length(object@norm.data), function(i){object@h5file.info[[i]]$file.path})}))
    rna_files = grep(paste0("_(sc10Xv3_|smartseq_|sn10Xv3_|sc10Xv2_)"), h5_files, value = TRUE)
    
    liger_cells = lapply(rna_files, function(i){
      rhdf5::h5read(i, "/matrix/barcodes")#change, extract from H5
    })
     
    liger_genes = lapply(rna_files, function(i){
      rhdf5::h5read(i, "/matrix/features")[[1]] #change, extract from H5
    })
    
    if(naive.clusters){
      clusters = readRDS(paste0(filepath,"/", region, "/", region, "_Deconvolution_Output/clusters_",rand.seed,"_object_clusters_naive.RDS"))
    } else {
      clusters = readRDS(paste0(filepath,"/", region, "/", region, "_Deconvolution_Output/clusters_",rand.seed,"_object_clusters.RDS"))
    }
    
    size_list = lapply(1:nlevels(clusters),function(x){vector(mode = "integer")})
    names(size_list) = levels(clusters)
    
    for(i in 1:length(rna_files)){
      
      n = rna_files[i]
      raw_data = Matrix::sparseMatrix(
        dims = c(length(liger_genes[[i]]), length(liger_cells[[i]])),
        i = as.numeric(rhdf5::h5read(n, "/matrix/indices")+1),
        p = as.numeric(rhdf5::h5read(n, "/matrix/indptr")),
        x = as.numeric(rhdf5::h5read(n, "/matrix/data"))
      )
      
      colnames(raw_data) = liger_cells[[i]]
      
      for(cell_type in names(size_list)){
       cells_subset = intersect(liger_cells[[i]], names(clusters)[clusters == cell_type])
       if(length(cells_subset)>1){
         size_list[[cell_type]]  = c(size_list[[cell_type]],Matrix::colSums(raw_data[,cells_subset]))
       } else if(length(cells_subset) == 1){
         size_list[[cell_type]] = c(size_list[[cell_type]],sum(raw_data[,cells_subset]))
       }
      }
    }
    
    cell_type_mean = sapply(size_list, mean)
    
    if(plot.hist){
      if(naive.clusters){
        pdf(file = paste0(filepath,"/", region, "/", region, "_Deconvolution_Output/histogram_size_",rand.seed,"_object_clusters_naive.PDF"), width = 6, height = 4)
      } else {
        pdf(file = paste0(filepath,"/", region, "/", region, "_Deconvolution_Output/histogram_size_",rand.seed,"_object_clusters.PDF"), width = 6, height = 4)
      }
      for(i in 1:length(size_list)){
        hist(size_list[[i]],main = paste0("Distribution of cell sizes for ",names(size_list)[i]), xlab = "Counts")
        abline(v = cell_type_mean[i], col = "red")
      }
      dev.off()
    }
    
    if(naive.clusters){
      saveRDS(cell_type_mean, paste0(filepath,"/", region, "/", region, "_Deconvolution_Output/cell_size_",rand.seed,"_object_clusters_naive.RDS"))
    } else {
      saveRDS(cell_type_mean, paste0(filepath,"/",region,"/", region, "_Deconvolution_Output/cell_size_",rand.seed,"_object_clusters.RDS"))
    }
}

select_defining_genes = function(
    filepath,
    region,
    clusters.from.atlas = TRUE,
    naive.clusters = FALSE,
    deconv.gene.num = 2000,
    gene.num.tol = 50,
    rand.seed = 123
    ){
  set.seed(rand.seed)
    
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
  }
      
  norm.data = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",descriptor,"_norm_data.RDS"))
  
  if(clusters.from.atlas){
     clusters = readRDS(paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/clusters_",descriptor,".RDS")))
  } else {
     clusters = readRDS(paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/user_defined_clusters_",descriptor,".RDS")))
  }
  
  message("Selecting genes with the KW test")
  
  chisq_list = list()
  for(i in 1:length(norm.data)){
    chisq_list[[i]] = matrixTests::row_kruskalwallis(as.matrix(norm.data[[i]]),as.vector(clusters[names(clusters) %in%colnames(norm.data[[i]])]))$statistic
    names(chisq_list[[i]]) = rownames(norm.data[[i]])
  }
  shared_genes = Reduce(c, lapply(norm.data, rownames))
  
  
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
      gene_vec = intersect(gene_vec, names(chisq_list[[i]][chisq_list[[i]] > quantile(chisq_list[[i]], var_thresh_new, na.rm = TRUE)]))
    }
    var_thresh_old = var_thresh_new
  }
  
  message(paste0(length(gene_vec), " genes found with p = ",var_thresh_old))

  saveRDS(list(chisq_vals = chisq_list, genes_used = gene_vec), paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/gene_selection_",descriptor,".RDS"))
  
}

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
