#' Convert proportions generated during the deconvolution into cell-type
#' assignments, for use with single cell spatial modalities (i.e. slideseq)
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param spatial.data.name A string, the name of the spatial dataset
#'
#' @return nothing
#'
#' @import
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }





#' Convert proportions generated during the deconvolution into cell-type
#' assignments, for use with single cell spatial modalities (i.e. slideseq)
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param layer.list A named list of spatial sample names, corresponding to
#'    groupings to be summarized by the function
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param plot A logical, corresponding with if to make bar plots corresponding
#'    to provided cell type and gene loading
#' @param type A string, either "mean" or "sum" corresponding with how to
#'    summarize across samples
#' @param mat.use A string, either "raw", "proportions", or "assignment",
#'    corresponding to the raw cell type loadings, the normalized loadings, or
#'    cell type assignments for single cell spatial modalities
#' @param use.cell.types A logical, if only the cell types provided to
#'    cell.types.use should be summarized, as opposed to all deconvolved.
#' @param cell.types.use A character vector of cell types to summarize
#' @param genes.use A character vector of genes to summarize
#'
#'
#' @return nothing
#'
#' @import
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
summarize_by_layer = function(
  filepath,
  region,
  layer.list,
  spatial.data.name,
  plot = FALSE,
  type = "mean",
  mat.use = "proportions",#"assignment
  clusters.from.atlas = T,
  naive.clusters = F,
  use.cell.types = TRUE,
  cell.types.use = NULL,
  cell.size = FALSE,
  genes.use = NULL,
  rand.seed = 123){
  
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")
  
  if(mat.use == "assignment"){
    assignments = readRDS(paste0(dir_spatial,"/deconvolution_max_factor_assignment_",descriptor,".RDS"))
    sub_vec = rep(0,nlevels(assignments))
    loadings = Reduce(rbind, lapply(assignments, function(x){subbed_vec = sub_vec; subbed_vec[as.numeric(x)] = 1; return(subbed_vec)}))
    colnames(loadings) = levels(assignments)
  } else {
    loadings = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))[[mat.use]]
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
      sub_loadings = as.matrix(loadings[rownames(loadings) %in% as.character(layer.list[[i]]), ])
      if(ncol(sub_loadings) == 1){sub_loadings = t(sub_loadings)}
      for(j in 1:length(cell.types.use)){
        if(type == "mean"){
          cell.type.matrix[i,j] = mean(sub_loadings[ ,cell.types.use[j]])
        } else if(type == "sum"){
          cell.type.matrix[i,j] = sum(sub_loadings[ ,cell.types.use[j]])
        }
      }
    }
    saveRDS(cell.type.matrix, paste0(dir_output,"/",spatial.data.name,"_cell_type_layer_summary_",descriptor,".RDS"))
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
    saveRDS(gene.matrix, paste0(dir_output,"/",spatial.data.name,"gene_layer_summary_",descriptor,".RDS"))
  }
  if(plot){
    dir_plots = paste0(dir_output,"/plots")
    if(!dir.exists(dir_plots)){
      dir.create(paste0(dir_plots))
      message("Created directory at ", dir_plots)
    }
    ggplot2::theme_set(cowplot::theme_cowplot())
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

      print(overall.cell.type.plot)

      ggplot2::ggsave(paste0(dir_plots, "/cell_type_layer_distribution_",descriptor,".PNG"),
                      overall.cell.type.plot,
                      width = 1000,
                      height = 800,
                      units = "px")
    }

    for(i in colnames(cell.type.matrix)){
      cell_type_consistent = sub("/", ".",sub(" ", "_",i))
      cell.type.df.sub = cell.type.df[cell.type.df$Cell_Types == i,]
      by.cell.type.plot = ggplot2::ggplot(cell.type.df.sub, ggplot2::aes(fill = Layers, x = Layers, y = Values)) +
        ggplot2::theme(text = ggplot2::element_text(size = 3),
                       axis.text = ggplot2::element_text(size = 2),
                       legend.position="none") +
        ggplot2::geom_bar(position = "dodge", stat = "identity") +
        ggplot2::xlab("Layer") +
        ggplot2::ylab("Value") +
        ggplot2::ggtitle(paste0("Distribution of ",cell_type_consistent, " cells by layer"))
        ggplot2::ggsave(paste0(dir_plots, "/",cell_type_consistent,"_layer_distribution_",descriptor,".PNG"),
                      by.cell.type.plot,
                      width = 500,
                      height = 400,
                      units = "px")
        message("Plotting ", cell_type_consistent)
    }
    if(!is.null(genes.use)){
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

        print(overall.gene.plot)

        ggplot2::ggsave(paste0(dir_plots, "/gene_layer_distribution_",descriptor,".PNG"),
                        overall.gene.plot,
                        width = 1000,
                        height = 800,
                        units = "px")
      }

      for(i in colnames(gene.matrix)){
        gene_consistent = sub("/", ".",sub(" ", "_",i))
        gene.df.sub = gene.df[gene.df$Genes == i,]
        by.gene.plot = ggplot2::ggplot(gene.df.sub, ggplot2::aes(fill = Layers, x = Layers, y = Values)) +
          ggplot2::theme(text = ggplot2::element_text(size = 3),
                         axis.text = ggplot2::element_text(size = 2),
                         legend.position="none") +
          ggplot2::geom_bar(position = "dodge", stat = "identity") +
          ggplot2::xlab("Layer") +
          ggplot2::ylab("Value") +
          ggplot2::ggtitle(paste0("Distribution of ",i, " expression by layer"))
        ggplot2::ggsave(paste0(dir_plots, "/",gene_consistent,"_gene_layer_distribution_",descriptor,".PNG"),
                        by.gene.plot,
                        width = 500,
                        height = 400,
                        units = "px")

      }
    }
  }
}

analyze_gene_signatures = function(filepath,
                                   region,
                                   spatial.data.name,
                                   clusters.from.atlas = TRUE,
                                   naive.clusters = FALSE,
                                   plot = FALSE,
                                   mat.use = "proportions",
                                   cell.types.use = NULL,
                                   cell.size = FALSE,
                                   rand.seed = 123){
  
  library(ggplot2)
  
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
    
  gene_sigs = readRDS(paste0(dir_spatial, "/gene_signature_output_",descriptor,".RDS"))
  genes = readRDS(paste0(dir_spatial,"/gene_selection_qc_",descriptor,".RDS"))
  
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
    dir_output = paste0(dir_spatial,"/",descriptor,"_output")

  
  
  signatures = gene_sigs$W
  signatures = signatures[rowSums(signatures) != 0 & rownames(signatures) != "",]
  
  if(!is.null(cell.types.use)){
    signatures = signatures[intersect(rownames(signatures),cell.types.use),]
  }
  
  cos_sim = lsa::cosine(t(signatures))
  
  
  cos_dist= as.dist(1- cos_sim)
  hierarchical_clust <- hclust(cos_dist, method = "ward.D2")
  
  saveRDS(list(cos_sim_signature = cos_sim, dendro_sig = hierarchical_clust),
          paste0(dir_output,"/gene_signature_analysis_summary_",descriptor,".RDS"))
  if(plot){
    
    dir_plots = paste0(dir_output,"/plots")
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
    
    print(heatmap_plot)
    print(dendro_plot)
   
    ggplot2::ggsave(paste0(dir_plots, "/cell_type_signature_heatmap_",descriptor,".PNG"),
                    heatmap_plot,
                    width = 1500,
                    height = 1200,
                    units = "px")
    ggplot2::ggsave(paste0(dir_plots, "/cell_type_signature_dendrogram_",descriptor,".PNG"),
                    dendro_plot,
                    width = 1000,
                    height = 800,
                    units = "px")
  }
}

cell_type_loading_histogram = function(
  filepath,
  region,
  spatial.data.name,
  rand.seed = 123,
  clusters.from.atlas = TRUE,
  naive.clusters = FALSE,
  cell.size = FALSE,
  mat.use = "proportions",
  cell.types.plot = NULL,
  print.plots = FALSE,
  bin.num = 30){

  library(ggplot2)
  
  set.seed(rand.seed)

  descriptor = as.character(rand.seed)
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }

  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")

  dir_gifs = paste0(dir_output,"/plots")
  if(!dir.exists(dir_gifs)){
    dir.create(paste0(dir_gifs))
    message("Created directory at ", dir_gifs)
  }

  loadings = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))[[mat.use]]
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


  pdf(file = paste0(dir_gifs, "/",descriptor,"_",mat.use,"_hist.PDF"), width = 7,height = 4)
  for(cell_type in names(hist_plot)){
    print(hist_plot[[cell_type]])
  }
  dev.off()
}

analyze_spatial_correlation = function(filepath,
                                       region,
                                       spatial.data.name,
                                       clusters.from.atlas = TRUE,
                                       naive.clusters = FALSE,
                                       cell.size = FALSE,
                                       plot = FALSE,
                                       mat.use = "proportions",
                                       cell.types.use = NULL,
                                       rand.seed = 123){
  
  library(ggplot2)
  
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")
  
  deconv_out = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))
  loadings = deconv_out[[mat.use]]
  loadings = loadings[, colSums(loadings) != 0 & colnames(loadings)!=""]
  
  if(!is.null(cell.types.use)){
    loadings = loadings[,intersect(colnames(loadings),cell.types.use)]
  }
  
  corr_sim_dist = cor(loadings)
  
  
  
  cor_dist= as.dist(1- corr_sim_dist)
  hierarchical_clust_dist <- hclust(cor_dist, method = "ward.D2")
  
  saveRDS(list(cos_sim_signature = corr_sim_dist, dendro_sig = hierarchical_clust_dist),
          paste0(dir_output,"/spatial_correlation_analysis_summary_",descriptor,".RDS"))
  
  if(plot){
    
    dir_plots = paste0(dir_output,"/plots")
    if(!dir.exists(dir_plots)){
      dir.create(paste0(dir_plots))
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

    print(heatmap_dist_plot)
    print(dendro_dist_plot)

    
    ggplot2::ggsave(paste0(dir_plots, "/spatial_correlation_heatmap_",descriptor,".PNG"),
                    heatmap_dist_plot,
                    width = 1500,
                    height = 1200,
                    units = "px")
    ggplot2::ggsave(paste0(dir_plots, "/spatial_correlation_dendrogram_",descriptor,".PNG"),
                    dendro_dist_plot,
                    width = 1000,
                    height = 800,
                    units = "px")
  }
}


calculate_wasserstein = function(
    filepath,
    region,
    spatial.data.name,
    clusters.from.atlas = TRUE,
    naive.clusters = FALSE,
    cell.size = FALSE,
    mat.use = "proportions",
    use.cell.types = TRUE,
    cell.types.use = NULL,
    genes.use = NULL,
    p = 2,
    min.samples = 1,
    plot = FALSE,
    rand.seed = 123){
  
  library(ggplot2)
  
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

  coords = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_coords_qc_",descriptor,".RDS"))
  
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")
  
  deconv_out = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))
  loadings = deconv_out[[mat.use]]
  loadings = loadings[, colSums(loadings) != 0 & colnames(loadings)!=""]
  
  
  loadings = loadings[, colSums(loadings) != 0]
  if(use.cell.types){
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
  }
  
  if(!is.null(genes.use)){
    exp = t(readRDS(paste0(dir_spatial,"/",spatial.data.name,"_exp_qc_",descriptor,".RDS")))
    exp = exp[, colnames(exp) %in% genes.use]
    exp[exp < 0] = 0
    exp = exp[rownames(loadings),]
    
    loadings = cbind(loadings, exp)
  }
  loadings = scale(loadings, center = F)
  library(transport)
  
  vars_1 = vars_2 = colnames(loadings)
  distance_mat = matrix(0L, nrow = ncol(loadings), ncol = ncol(loadings))
  colnames(distance_mat) = rownames(distance_mat) = vars_1
  loadings = apply(loadings, MARGIN = 2, function(x){x/sum(x)})
  for(var_1 in vars_1){
    mass_1 =  loadings[,var_1]
    distribution_1 = wpp(coords,mass_1)
    for(var_2 in vars_2){
      mass_2 =  loadings[,var_2]
      distribution_2 = wpp(coords,mass_2)
      distance_mat[var_1, var_2] = 
        distance_mat[var_2, var_1] = 
        wasserstein(distribution_1, distribution_2, p)
    }
    vars_2 = vars_2[2:length(vars_2)]
  }
  
  saveRDS(distance_mat,
          paste0(dir_output,"/wasserstein_distance_mat_",descriptor,".RDS"))
  
  if(plot){
    
    dir_plots = paste0(dir_output,"/plots")
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
  }
  
}
