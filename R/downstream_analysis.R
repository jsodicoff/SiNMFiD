#' Summarize cell-type and gene expression data by
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#' @param type A string, either "mean" or "sum", how results should be
#'     combined for summary
#' @param layer.list A named list of spatial samples by layer of interest
#' @param mat.use A string, either "raw", "proportions", or "assignments"
#'     referring to what version of the results to summarize
#' @param cell.types.use A string of cell type labels to include in the plot,
#'     by default all cell types present
#' @param genes.use A string of genes to include in a plot, by default none
#' @param return.objs Logical, whether to return a list of matrices of
#'     derived data
#'
#' @return cell-type and gene expression data summarized by layer in a named
#'     list, if `return.objs = TRUE`
#'
#' @export
summarize_by_layer = function(
    filepath,
    analysis.name,
    spatial.data.name,
    rand.seed = 123,
    layer.list,
    type = "mean",
    mat.use = "proportions",
    cell.types.use = NULL,
    genes.use = NULL,
    return.objs = FALSE){

  set.seed(rand.seed)

  dir_spatial = file.path(filepath,analysis.name,spatial.data.name)
  dir_output = file.path(filepath,analysis.name,rand.seed,"downstream_output")

  if(mat.use == "assignment"){
    #will be a whole thing
    assignments = readRDS(paste0(dir_spatial,"/deconvolution_max_factor_assignment_",descriptor,".RDS"))
    sub_vec = rep(0,nlevels(assignments))
    loadings = Reduce(rbind, lapply(assignments, function(x){
      subbed_vec = sub_vec; subbed_vec[as.numeric(x)] = 1; return(subbed_vec)
      }))
    colnames(loadings) = levels(assignments)
  } else {
    loadings = readRDS(file.path(filepath,analysis.name,rand.seed,"deconvolution_output.RDS"))[[mat.use]]
  }

  if(return.objs){
    return.obj = list()
  }

  if(use.cell.types){
    if(!is.null(cell.types.use)){
      cell.types.use = intersect(cell.types.use, colnames(loadings))
    } else {
      cell.types.use = colnames(loadings)
    }
    cell.types.use = cell.types.use[cell.types.use != ""]
    cell.type.matrix = matrix(0L, nrow = length(layer.list), ncol = length(cell.types.use))
    rownames(cell.type.matrix) = names(layer.list)
    colnames(cell.type.matrix) = cell.types.use
    for(i in 1:length(layer.list)){
      sub_loadings = as.matrix(loadings[rownames(loadings) %in% as.character(layer.list[[i]]),])
      if(ncol(sub_loadings) == 1){sub_loadings = t(sub_loadings)}
      for(j in 1:length(cell.types.use)){
        if(type == "mean"){
          cell.type.matrix[i,j] = mean(sub_loadings[ ,cell.types.use[j]])
        } else if(type == "sum"){
          cell.type.matrix[i,j] = sum(sub_loadings[ ,cell.types.use[j]])
        }
      }
    }
    saveRDS(cell.type.matrix, file.path(dir_output,"cell_type_layer_summary.RDS"))
    if(return.objs){
      return.obj[["cell.type.matrix"]] = cell.type.matrix
    }
  }
  if(!is.null(genes.use)){
    spatial.data = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_exp.RDS"))
    genes.use = intersect(genes.use, rownames(spatial.data))
    spatial.data[is.na(spatial.data)] = 0
    spatial.data = t(scale(t(spatial.data[genes.use,]), center = FALSE))
    spatial.data[spatial.data < 0 ] = 0
    gene.matrix = matrix(0L, nrow = length(layer.list), ncol = length(genes.use))
    rownames(gene.matrix) = names(layer.list)
    colnames(gene.matrix) = genes.use
    for(i in 1:length(layer.list)){
      sub_loadings = as.matrix(spatial.data[ ,colnames(spatial.data) %in% as.character(layer.list[[i]])])
      if(ncol(sub_loadings) == 1){sub_loadings = t(sub_loadings)}
      for(j in 1:length(genes.use)){
        if(type == "mean"){
          gene.matrix[i,j] = mean(sub_loadings[genes.use[j],])
        } else if(type == "sum"){
          gene.matrix[i,j] = sum(sub_loadings[genes.use[j],])
        }
      }
    }
    saveRDS(gene.matrix, file.path(dir_output,"gene_layer_summary.RDS"))
    if(return.objs){
      return.obj[["gene.matrix"]] = gene.matrix
    }
  }
  if(return.objs){
    return(return.obj)
  }
}


#' Plot results of `summarize_by_layer`
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#' @param print.plots Logical, whether to display results in the plots panel
#'
#' @return nothing
#'
#' @import ggplot2
#'
#' @export
plot_summarize_by_layer = function(
    filepath,
    analysis.name,
    spatial.data.name,
    rand.seed = 123,
    print.plots = T){

  dir_output = file.path(filepath,analysis.name,rand.seed,"downstream_output")
  dir_plots = file.path(dir_output,"plots")
  if(!dir.exists(dir_plots)){
    dir.create(paste0(dir_plots))
    message("Created directory at ", dir_plots)
  }
  ggplot2::theme_set(cowplot::theme_cowplot())

  cell_type = file.path(dir_output,"cell_type_layer_summary.RDS")

  if(file.exists(cell_type)){
    cell.type.matrix = readRDS(cell_type)
    labels.cell.type = expand.grid(rownames(cell.type.matrix), colnames(cell.type.matrix))
    cell.type.df = data.frame(Layers = as.character(labels.cell.type[,1]),
                              Cell_Types = as.character(labels.cell.type[,2]),
                              Values = as.vector(cell.type.matrix))

    if(ncol(cell.type.matrix) > 1){
      overall.cell.type.plot = ggplot2::ggplot(cell.type.df, ggplot2::aes(fill = Cell_Types, y = Values, x = Layers)) +
        ggplot2::theme(text = ggplot2::element_text(size = 10),
                       axis.text = ggplot2::element_text(size = 5),
                       legend.title = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(size = 3),
                       legend.key.height = ggplot2::unit(3, 'mm'),
                       legend.key.width = ggplot2::unit(1, 'mm')) +
        ggplot2::geom_bar(position = "dodge", stat = "identity") +
        ggplot2::xlab("Layer") +
        ggplot2::ylab("Value") +
        ggplot2::ggtitle(paste0("Distribution of cell types by layer"))

      if(print.plots){
        print(overall.cell.type.plot)
      }

      ggplot2::ggsave(file.path(dir_plots, "cell_type_layer_distribution.PNG"),
                      overall.cell.type.plot,
                      width = 1000,
                      height = 800,
                      units = "px")
    }

    for(i in colnames(cell.type.matrix)){
      cell_type_consistent = sub("/", ".",sub(" ", "_",i))
      message("Plotting ", cell_type_consistent)
      cell.type.df.sub = cell.type.df[cell.type.df$Cell_Types == i,]
      by.cell.type.plot = ggplot2::ggplot(cell.type.df.sub, ggplot2::aes(fill = Layers, x = Layers, y = Values)) +
        ggplot2::theme(text = ggplot2::element_text(size = 3),
                       axis.text = ggplot2::element_text(size = 2),
                       legend.position="none") +
        ggplot2::geom_bar(position = "dodge", stat = "identity") +
        ggplot2::xlab("Layer") +
        ggplot2::ylab("Value") +
        ggplot2::ggtitle(paste0("Distribution of ",cell_type_consistent, " cells by layer"))

      if(print.plots){
        print(by.cell.type.plot)
      }

      ggplot2::ggsave(file.path(dir_plots, paste0(cell_type_consistent,"_layer_distribution.PNG")),
                      by.cell.type.plot,
                      width = 500,
                      height = 400,
                      units = "px")
    }
  }

  gene = file.path(dir_output,"gene_layer_summary.RDS")

  if(file.exists(gene)){
    gene.matrix = readRDS(gene)
    labels.genes = expand.grid(rownames(gene.matrix), colnames(gene.matrix))
    gene.df = data.frame(Layers = as.character(labels.genes[,1]),
                         Genes = as.character(labels.genes[,2]),
                         Values = as.vector(gene.matrix))
    if(length(genes.use) > 1){
      overall.gene.plot = ggplot2::ggplot(gene.df, ggplot2::aes(fill = Genes, y = Values, x = Layers)) +
        ggplot2::theme(text = ggplot2::element_text(size = 8),
                       axis.text = ggplot2::element_text(size = 5),
                       legend.title = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(size = 3),
                       legend.key.height = ggplot2::unit(3, 'mm'),
                       legend.key.width = ggplot2::unit(1, 'mm')) +
        ggplot2::geom_bar(position = "dodge", stat = "identity") +
        ggplot2::xlab("Layer") +
        ggplot2::ylab("Value") +
        ggplot2::ggtitle(paste0("Distribution of gene expression by layer"))

      if(print.plots){
        print(overall.gene.plot)
      }

      ggplot2::ggsave(file.path(dir_plots, "gene_layer_distribution.PNG"),
                      overall.gene.plot,
                      width = 1000,
                      height = 800,
                      units = "px")
    }

    for(i in colnames(gene.matrix)){
      gene_consistent = sub("/", ".",sub(" ", "_",i))
      message("Plotting ", gene_consistent)
      gene.df.sub = gene.df[gene.df$Genes == i,]
      by.gene.plot = ggplot2::ggplot(gene.df.sub, ggplot2::aes(fill = Layers, x = Layers, y = Values)) +
        ggplot2::theme(text = ggplot2::element_text(size = 3),
                       axis.text = ggplot2::element_text(size = 2),
                       legend.position="none") +
        ggplot2::geom_bar(position = "dodge", stat = "identity") +
        ggplot2::xlab("Layer") +
        ggplot2::ylab("Value") +
        ggplot2::ggtitle(paste0("Distribution of ",i, " expression by layer"))

      if(print.plots){
        print(by.gene.plot)
      }

      ggplot2::ggsave(file.path(dir_plots, paste0(gene_consistent,"_gene_layer_distribution.PNG")),
                      by.gene.plot,
                      width = 500,
                      height = 400,
                      units = "px")

    }
  }
}

#' Calculate relationships between cell types
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#' @param cell.types.use A string of cell type labels to include in the plot,
#'     by default all cell types present
#' @param return.objs Logical, whether to return a list of matrices of
#'     derived data
#'
#' @return named list of cosine similarity matrix and hierarchical clustering,
#'     if `return.objs = TRUE`
#'
#' @import lsa
#'
#' @export
analyze_gene_signatures = function(filepath,
                                   analysis.name,
                                   spatial.data.name,
                                   rand.seed = 123,
                                   cell.types.use = NULL,
                                   return.objs = F){

  set.seed(rand.seed)

  dir_spatial = file.path(filepath,analysis.name,spatial.data.name,rand.seed)

  gene_sigs = readRDS(file.path(dir_spatial, "gene_signature_output.RDS"))
  genes = readRDS(file.path(dir_spatial,"gene_selection_qc.RDS"))

  dir_output = file.path(dir_spatial,"downstream_output")

  signatures = gene_sigs$W
  signatures = signatures[rowSums(signatures) != 0 & rownames(signatures) != "",]

  if(!is.null(cell.types.use)){
    signatures = signatures[intersect(rownames(signatures),cell.types.use),]
  }

  cos_sim = lsa::cosine(t(signatures))
  cos_dist= as.dist(1- cos_sim)
  hierarchical_clust <- hclust(cos_dist, method = "ward.D2")

  obj = list(cos_sim_signature = cos_sim, dendro_sig = hierarchical_clust)
  saveRDS(obj,
          file.path(dir_output,"gene_signature_analysis_summary.RDS"))

  if(return.objs){
    return(obj)
  }
}

#' Plot results of `analyze_gene_signatures`
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#' @param print.plots Logical, whether to display results in the plots panel
#'
#' @return nothing
#'
#' @import ggplot2
#' @import viridis
#'
#' @export
plot_analyze_gene_signatures = function(filepath,
                                   analysis.name,
                                   spatial.data.name,
                                   rand.seed = 123,
                                   print.plots = T){

  dir_spatial = file.path(filepath,analysis.name,spatial.data.name,rand.seed)
  dir_output = file.path(dir_spatial,"downstream_output")

  gene_sig_out = readRDS(file.path(dir_output,"gene_signature_analysis_summary.RDS"))
  cos_sim = gene_sig_out[[1]]
  hierarchical_clust = gene_sig_out[[2]]


  dir_plots = file.path(dir_output,"plots")
  if(!dir.exists(dir_plots)){
    dir.create(paste0(dir_plots))
    message("Created directory at ", dir_plots)
  }

  heatmap_df = data.frame(expand.grid(Cell_Type_1 = rownames(cos_sim),
                                      Cell_Type_2 = colnames(cos_sim)),
                          cos_sim = as.vector(cos_sim))
  heatmap_plot = ggplot2::ggplot(heatmap_df, ggplot2::aes(x = Cell_Type_1, y = Cell_Type_2, fill = cos_sim)) +
    ggplot2::labs(y = "Cell Types", fill = "", title = "Cosine Similarity for Cell Type Signatures") +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                   axis.title.x = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 8),
                   axis.text = ggplot2::element_text(size = 5),
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 3),
                   legend.key.height = ggplot2::unit(3, 'mm'),
                   legend.key.width = ggplot2::unit(1, 'mm')) +
    ggplot2::geom_tile() +
    viridis::scale_fill_viridis()

  dendro_plot = ggdendro::ggdendrogram(hierarchical_clust, rotate = FALSE, size = 2) +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 4),
                   axis.text = ggplot2::element_text(size = 4)) +
    ggplot2::ggtitle("Hierarchical Clustering", subtitle = "By Cell Type Signature")
  if(print.plots){
    print(heatmap_plot)
    print(dendro_plot)
  }

  ggplot2::ggsave(file.path(dir_plots, "cell_type_signature_heatmap.PNG"),
                  heatmap_plot,
                  width = 1500,
                  height = 1200,
                  units = "px")
  ggplot2::ggsave(file.path(dir_plots, "cell_type_signature_dendrogram.PNG"),
                  dendro_plot,
                  width = 1000,
                  height = 800,
                  units = "px")
}

#' Generate histograms of loading by cell type
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#' @param mat.use A string, either "raw" or "proportions"
#'     referring to what version of the results to summarize
#' @param cell.types.use A string of cell type labels to include in the plot,
#'     by default all cell types present
#' @param print.plots Logical, whether to display results in the plots panel
#' @param bin.num Integer number of bins to use in histogram
#'
#' @return nothing
#'
#' @import ggplot2
#'
#' @export
cell_type_loading_histogram = function(
    filepath,
    analysis.name,
    spatial.data.name,
    rand.seed = 123,
    mat.use = "proportions",
    cell.types.plot = NULL,
    print.plots = TRUE,
    bin.num = 30){

  set.seed(rand.seed)

  dir_spatial = file.path(filepath,analysis.name,spatial.data.name,rand.seed)
  dir_output = file.path(dir_spatial,"downstream_output")
  dir_gifs = file.path(dir_output,"plots")

  if(!dir.exists(dir_gifs)){
    dir.create(paste0(dir_gifs))
    message("Created directory at ", dir_gifs)
  }

  loadings = readRDS(filepath(dir_spatial,"deconvolution_output.RDS"))[[mat.use]]
  cell_types = colnames(loadings)

  if(!is.null(cell.types.plot)){
    cell.types.plot = intersect(cell.types.plot, cell_types)
  }
  cell.types.plot = gsub(" ", ".", gsub("-","_",cell.types.plot))
  colnames(loadings) = gsub(" ", ".", gsub("-","_",colnames(loadings)))

  loadings = as.data.frame(loadings[,cell.types.plot])

  hist_plot = list()
  for(cell_type in colnames(loadings)){
    range = loadings[,cell_type][2]
    hist_plot[[cell_type]] = ggplot(loadings, aes_string(x=cell_type)) +
      geom_histogram(bins = bin.num, fill = rainbow(bin.num)) +
      labs(y = "Count", x = "",title = paste0("Histogram of ",cell_type, " loading by voxel"))
    if(print.plots){
      print(hist_plot[[cell_type]])
    }
  }


  pdf(file = file.path(dir_gifs, paste0(mat.use,"_hist.PDF")), width = 7,height = 4)
  for(cell_type in names(hist_plot)){
    print(hist_plot[[cell_type]])
  }
  dev.off()
}

#' Calculate relationships between cell type distributions
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#' @param mat.use A string, either "raw" or "proportions"
#'     referring to what version of the results to summarize
#' @param cell.types.use A string of cell type labels to include in the plot,
#'     by default all cell types present
#' @param return.objs Logical, whether to return a list of matrices of
#'     derived data
#'
#' @return named list of pearson correlation matrix and hierarchical clustering,
#'     if `return.objs = TRUE`
#'
#' @import ggplot2
#'
#' @export
analyze_spatial_correlation = function(filepath,
                                       analysis.name,
                                       spatial.data.name,
                                       rand.seed = 123,
                                       mat.use = "proportions",
                                       cell.types.use = NULL,
                                       return.objs = F){

  set.seed(rand.seed)

  dir_spatial = file.path(filepath,analysis.name,spatial.data.name,rand.seed)
  dir_output = file.path(dir_spatial,"downstream_output")

  loadings = readRDS(file.path(dir_spatial,"deconvolution_output.RDS"))[[mat.use]]
  loadings = loadings[, colSums(loadings) != 0 & colnames(loadings)!=""]

  if(!is.null(cell.types.use)){
    loadings = loadings[,intersect(colnames(loadings),cell.types.use)]
  }

  corr_sim_dist = cor(loadings)

  cor_dist= as.dist(1- corr_sim_dist)
  hierarchical_clust_dist <- hclust(cor_dist, method = "ward.D2")

  obj = list(cos_sim_signature = corr_sim_dist, dendro_sig = hierarchical_clust_dist)

  saveRDS(obj,
          file.path(dir_output,"spatial_correlation_analysis_summary.RDS"))

  if(return.objs){
    return(obj)
  }

}

#' Plot results of `analyze_spatial_correlation`
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#' @param print.plots Logical, whether to display results in the plots panel
#'
#' @return nothing
#'
#' @import ggplot2
#' @import viridis
#'
#' @export
plot_analyze_spatial_correlation = function(filepath,
                                       analysis.name,
                                       spatial.data.name,
                                       rand.seed = 123,
                                       print.plots = TRUE){

  set.seed(rand.seed)

  dir_spatial = file.path(filepath,analysis.name,spatial.data.name,rand.seed)
  dir_output = file.path(dir_spatial,"downstream_output")
  dir_plots = file.path(dir_output,"plots")

  obj = readRDS(file.path(dir_output,"spatial_correlation_analysis_summary.RDS"))
  corr_sim_dist = obj[[1]]
  hierarchical_clust_dist = obj[[2]]

  if(!dir.exists(dir_plots)){
    dir.create(dir_plots)
    message("Created directory at ", dir_plots)
  }

  heatmap_dist_df = data.frame(expand.grid(Cell_Type_1 = rownames(corr_sim_dist),
                                           Cell_Type_2 = colnames(corr_sim_dist)),
                               cos_sim = as.vector(corr_sim_dist))
  heatmap_dist_plot = ggplot2::ggplot(heatmap_dist_df, ggplot2::aes(x = Cell_Type_1, y = Cell_Type_2, fill = cos_sim)) +
    ggplot2::labs(y = "Cell Types", fill = "", title = "Correlation for Cell Type Distribution") +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                   axis.title.x = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 8),
                   axis.text = ggplot2::element_text(size = 5),
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 3),
                   legend.key.height = ggplot2::unit(3, 'mm'),
                   legend.key.width = ggplot2::unit(1, 'mm')) +
    ggplot2::geom_tile() +
    viridis::scale_fill_viridis()

  dendro_dist_plot = ggdendro::ggdendrogram(hierarchical_clust_dist, rotate = FALSE, size = 2) +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 4),
                   axis.text = ggplot2::element_text(size = 4)) +
    ggplot2::ggtitle("Hierarchical Clustering", subtitle = "By Cell Type Distribution")

  if(print.plots){
    print(heatmap_dist_plot)
    print(dendro_dist_plot)
  }


  ggplot2::ggsave(file.path(dir_plots, "spatial_correlation_heatmap.PNG"),
                  heatmap_dist_plot,
                  width = 1500,
                  height = 1200,
                  units = "px")
  ggplot2::ggsave(file.path(dir_plots, "spatial_correlation_dendrogram.PNG"),
                  dendro_dist_plot,
                  width = 1000,
                  height = 800,
                  units = "px")
}

#' Calculate the Wasserstein distance between cell-types and genes
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#' @param mat.use A string, either "raw" or "proportions"
#'     referring to what version of the results to summarize
#' @param cell.types.use A string of cell type labels to include in the plot,
#'     by default all cell types present
#' @param genes.use A string of genes to include in a plot, by default none
#' @param p The `p` exponent used for the Minkowski distance
#' @param min.samples Integer value, the minimum number of samples a cell type
#'     can load on and be included in the analysis
#' @param return.objs Logical, whether to return a list of matrices of
#'     derived data
#'
#' @return matrix of pairwise Wasserstein distances if `return.objs = TRUE`
#'
#' @import transport
#'
#' @export
calculate_wasserstein = function(
    filepath,
    analysis.name,
    spatial.data.name,
    rand.seed = 123,
    mat.use = "proportions",
    cell.types.use = NULL,
    genes.use = NULL,
    p = 2,
    min.samples = 1,
    return.objs = F){

  set.seed(rand.seed)

  dir_spatial = file.path(filepath,analysis.name,spatial.data.name,rand.seed)
  dir_output = file.path(dir_spatial,"downstream_output")

  coords = readRDS(file.path(dir_spatial,"coords_qc.RDS"))

  deconv_out = readRDS(file.path(dir_spatial,"deconvolution_output.RDS"))
  loadings = deconv_out[[mat.use]]
  loadings = loadings[, colSums(loadings) != 0 & colnames(loadings)!=""]

  loadings = loadings[, colSums(loadings) != 0]
  if(!is.null(cell.types.use)){
    cell.types.use = intersect(cell.types.use, colnames(loadings))
  } else {
    cell.types.use = colnames(loadings)
  }
  cell.types.use = cell.types.use[cell.types.use != "" & apply(loadings, MARGIN = 2, function(x){sum(x != 0)}) > min.samples]
  loadings = loadings[rownames(loadings) %in% rownames(coords), cell.types.use]
  
  if(is.vector(loadings)){
    new_loadings = matrix(loadings)
    rownames(new_loadings) = names(loadings)
    colnames(new_loadings) = cell.types.use
    loadings = new_loadings
    rm(new_loadings)
  }

  colnames(loadings) = sub("/",".",sub(" ", "_", cell.types.use))

  coords = coords[rownames(loadings),]
  

  if(!is.null(genes.use)){
    exp = t(readRDS(file.path(dir_spatial, "exp_qc.RDS")))
    exp = exp[, colnames(exp) %in% genes.use]
    exp[exp < 0] = 0
    exp = exp[rownames(loadings),]

    loadings = cbind(loadings, exp)
  }
  loadings = scale(loadings, center = F)

  vars_1 = vars_2 = colnames(loadings)
  distance_mat = matrix(0L, nrow = ncol(loadings), ncol = ncol(loadings))
  colnames(distance_mat) = rownames(distance_mat) = vars_1
  loadings = apply(loadings, MARGIN = 2, function(x){x/sum(x)})
  for(var_1 in vars_1){
    mass_1 =  loadings[,var_1]
    distribution_1 = transport::wpp(coords,mass_1)
    for(var_2 in vars_2){
      mass_2 =  loadings[,var_2]
      distribution_2 = transport::wpp(coords,mass_2)
      distance_mat[var_1, var_2] =
        distance_mat[var_2, var_1] =
        transport::wasserstein(distribution_1, distribution_2, p)
    }
    vars_2 = vars_2[2:length(vars_2)]
  }

  saveRDS(distance_mat,
          file.path(dir_output,"wasserstein_distance_mat.RDS"))

  if(return.objs){
    return(distance_mat)
  }

}

#' Plot results of `calculate_wasserstein`
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param spatial.data.name String identifying the spatial sample
#' @param rand.seed Integer random seed
#' @param print.plots Logical, whether to display results in the plots panel
#'
#' @return nothing
#'
#' @import ggplot2
#' @import viridis
#'
#' @export
plot_calculate_wasserstein = function(
    filepath,
    analysis.name,
    spatial.data.name,
    rand.seed = 123,
    print.plots = T){

  set.seed(rand.seed)

  dir_spatial = file.path(filepath,analysis.name,spatial.data.name,rand.seed)
  dir_output = file.path(dir_spatial,"downstream_output")
  dir_plots = file.path(dir_output,"plots")

  distance_mat = readRDS(file.path(dir_output,"wasserstein_distance_mat.RDS"))

  if(!dir.exists(dir_plots)){
    dir.create(paste0(dir_plots))
    message("Created directory at ", dir_plots)
  }

  heatmap_wasserstein_df = data.frame(expand.grid(Cell_Type_1 = rownames(distance_mat),
                                                  Cell_Type_2 = colnames(distance_mat)),
                                      wasserstein_dist = as.vector(distance_mat))
  heatmap_wasserstein_plot = ggplot2::ggplot(heatmap_wasserstein_df, ggplot2::aes(x = Cell_Type_1, y = Cell_Type_2, fill = wasserstein_dist)) +
    ggplot2::labs(y = "Cell Types", fill = "", title = "Wasserstein Distance by Cell Type and Gene") +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                   axis.title.x = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 8),
                   axis.text = ggplot2::element_text(size = 5),
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 3),
                   legend.key.height = ggplot2::unit(3, 'mm'),
                   legend.key.width = ggplot2::unit(1, 'mm')) +
    ggplot2::geom_tile() +
    viridis::scale_fill_viridis()

  ggplot2::ggsave(paste0(dir_plots, "/wasserstein_heatmap_",descriptor,".PNG"),
                  heatmap_wasserstein_plot,
                  width = 1500,
                  height = 1200,
                  units = "px")
  if(print.plots){
    print(heatmap_wasserstein_plot)
  }

}

#' Summarize subregions of a vector of regions of interest
#'
#' @param regions A vector of region names
#' @param ontology.file A csv describing the Allen structure ontology
#' @param return.objs Logical, whether to return acronyms for all subregions found
#'
#' @return A vector of unique subregions within the provided regions,
#'     if `return.objs = TRUE`
#' @export
summarize_subregions = function(regions,
                                ontology.file = "Downloads/allen_structure_ontology.csv",
                                return.objs = F){
  allen_structure = read.csv(ontology.file,header=T)
  allen_structure_list = lapply(allen_structure$structure_id_path, function(x){strsplit(x, "/")}[[1]])

  subregion_vec = vector(mode = "character")
  for(region in regions){
    message(region)
    sub_allen = allen_structure[allen_structure$acronym == region, ]
    if(nrow(sub_allen)>1){
      id = sub_allen$id[which.max(sub_allen$depth)]
    } else {
      id = sub_allen$id[1]
    }
    if(length(id) != 0){
      subregions = allen_structure[sapply(allen_structure_list, function(x){id %in% x}) & grepl(region, allen_structure$acronym),]
      print(subregions[,c("acronym","name")])
      subregion_vec = c(subregion_vec, subregions[,"acronym"])
    } else {
      message("No subregions!")
    }
  }
  if(return.objs){
    return(unique(subregion_vec))
  }
}

#' Summarize cell types present in the source annotations
#'
#' @param filepath Path to analysis directory
#' @param analysis.name String identifying the analysis
#' @param return.objs  Logical, whether to return a vector of the names of clusters
#'
#' @return A vector of unique clusters in the source annotations,
#'     if `return.objs = TRUE`
#' @export
summarize_clusters = function(filepath,
                              analysis.name,
                              return.objs = F){
  clusts = readRDS(file.path(filepath,analysis.name,"source_annotations.RDS"))
  message(paste0(analysis.name, " cluster frequency"))
  print(table(clusts))
  if(return.objs){
    return(levels(clusts))
  }
}
