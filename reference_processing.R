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
