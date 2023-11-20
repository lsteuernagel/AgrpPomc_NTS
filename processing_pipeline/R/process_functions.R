

##########
### Filter object
##########

#' Apply basic QC to seurat object
#'
#' Filter nCount_RNA, nFeature_RNA and percent_mt
#'
#' @param seurat_object Seurat object
#' @param parameter_list parameters list with all relevant parameters
#' @param features_exclude_list vector of gene names to exclude
#'
#' @return seurat object after filtering
#'

filter_object = function(seurat_object,parameter_list,features_exclude_list){

  # add exclude features
  #seurat_object[["percent_exclude_features"]] <- Seurat::PercentageFeatureSet(seurat_object,features=features_exclude_list[features_exclude_list %in% rownames(seurat_object)])
  seurat_object[["percent_mt"]]<- Seurat::PercentageFeatureSet(seurat_object, pattern = "^mt-")
  seurat_object[["percent_exclude_features"]] <- Matrix::colSums(seurat_object@assays$RNA@counts[features_exclude_list[features_exclude_list %in% rownames(seurat_object)],]) / seurat_object@meta.data$nCount_RNA

  # max umi
  #maxUMI_use = as.numeric(min(median(seurat_object@meta.data$nCount_RNA)*parameter_list$maxUMI_dynamic,parameter_list$maxUMI))
  # subset seurat:
  seurat_object = subset(seurat_object,subset = nCount_RNA <= parameter_list$maxUMI &
                           nCount_RNA >= parameter_list$minUMI &
                           nFeature_RNA >= parameter_list$minFeatures &
                           percent_mt <= parameter_list$max_mt)

  # if there are samples with very few cells: remove
  keep_samples=names(table(seurat_object@meta.data[,parameter_list$sample_column]))[table(seurat_object@meta.data[,parameter_list$sample_column]) > parameter_list$min_cells_sample]
  seurat_object@meta.data$tmp_id = seurat_object@meta.data[,parameter_list$sample_column]
  seurat_object = subset(seurat_object, subset = tmp_id %in% keep_samples)
  seurat_object@meta.data = seurat_object@meta.data[,!grepl("tmp_id",colnames(seurat_object@meta.data))]

  return(seurat_object)

}

##########
### Doublet detection
##########

#' Detect Doublets in seurat object
#'
#' Runs scDblFinder::scDblFinder per sample
#'
#' @param seurat_object Seurat object
#' @param parameter_list parameters list with all relevant parameters
#' @param features_exclude_list vector of gene names to exclude
#'
#' @return seurat object after doublet removal
#'

detect_doublets = function(seurat_object,parameter_list,features_exclude_list){

  if(parameter_list$doublet_rate == "default"){
    dbr = NULL
  }else{
    dbr = parameter_list$doublet_rate
  }

  all_samples = Seurat::SplitObject(seurat_object,split.by = parameter_list$sample_column)
  removed_cells = c()
  for(i in 1:length(all_samples)){
    sample_sce = Seurat::as.SingleCellExperiment(all_samples[[i]])
    sample_sce <- scDblFinder::scDblFinder(sample_sce, clusters = parameter_list$doublet_clusters ,includePCs = parameter_list$includePCs, dbr =dbr, nfeatures = parameter_list$nfeatures,verbose =FALSE) # dbr = parameter_list$doublet_rate,
    all_samples[[i]]@meta.data$scDblFinder = SummarizedExperiment::colData(sample_sce)$scDblFinder.class
    if(parameter_list$doublet_remove){
      all_samples[[i]] = subset(all_samples[[i]],cells = colnames(all_samples[[i]])[which(all_samples[[i]]@meta.data$scDblFinder == 'singlet')])
      removed_cells = c(removed_cells,length(which(all_samples[[i]]@meta.data$scDblFinder != 'singlet')))
    }
  }
  # merge
  seurat_object = merge(x = all_samples[[1]], y = all_samples[2:length(all_samples)])
  seurat_object@misc$removed_doublets = sum(removed_cells)

  return(seurat_object)
}


##########
### read_embedding
##########

#' Load an emebedding with cells x lowDims from flatfile, ensuring consistency with a Seurat object (or metadata only for faster usage)
#' @param filename_withpath filepath
#' @param seurat_object seuratobject associated with current embedding. If specified metadata does not have to be set explicitly.
#' @param seurat_object_metadata metadata only of seuratobject associated with current embedding
#' @return

read_embedding = function(filename_withpath,seurat_object=NULL,seurat_object_metadata=NULL){

  #get metadata
  if(!is.null(seurat_object_metadata)){
    metadata = seurat_object_metadata
  }else{
    if(!is.null(seurat_object)){
      metadata = seurat_object@meta.data
    }else{
      stop("Please provide either a dataframe with metadata or a seurat object with metadata that can be exctracted!")
    }
  }
  # load
  current_embedding = data.table::fread(filename_withpath,data.table = F)
  # use first col as rownames
  if(is.character(current_embedding[,1])){
    rnames = current_embedding[,1]
    current_embedding = current_embedding[,2:ncol(current_embedding)]
    rownames(current_embedding)=rnames
    # reorder to align with rest of object
    if(any(is.na(match(rownames(metadata),rownames(current_embedding))))){
      message("Found ",length(rnames)," rows in embedding and ",length(rownames(metadata))," rows in metadata.")
      stop("Cell names from loaded reduction and new object are not matching exactly. Stopping import.")
    }
    current_embedding = current_embedding[match(rownames(metadata),rownames(current_embedding)),]
  }else{
    warning("First column of loaded file is not of type character, using rownames of metadata as rownames of added reduction. This can induce bugs if the order changed due to split/merge of the Seurat object!")
    rownames(current_embedding) = rownames(metadata)
  }
  return(current_embedding)
}
