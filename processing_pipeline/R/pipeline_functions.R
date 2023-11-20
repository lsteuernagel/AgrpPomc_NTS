
# helper
writeList_to_JSON = function (list_with_rows, filename) {
  jsonfile = jsonlite::toJSON(list_with_rows, pretty = TRUE,
                              auto_unbox = TRUE, digits = NA)
  writeLines(jsonfile, con = paste0(filename))
}

##########
### scseq_parameters
##########

#' Set the parameters for the single cell pipeline
#'
#' Runs scDblFinder::scDblFinder per sample
#'
#' @param param_file where to export json version of parameter list to
#'
#' @return list of parameters for scseq_pipeline
#'

scseq_parameters = function(param_file = NULL,
                            project_path,
                            input_seurat,
                            project_name = "project",
                            genes_to_exclude_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/scseq_pipeline/data/features_exclude_list2.json",
                            n_cores = 50,
                            # filtering
                            feature_set_size = c(2000),
                            id_column = "Cell_ID",
                            global_seed = 123456,
                            minUMI = 1000,
                            minFeatures = 500,
                            maxUMI = 100000,
                            max_mt = 10,
                            min_cells_sample = 100,
                            sample_column = "Sample_ID",
                            k_param = 30,
                            dist_type="cosine",
                            # doublet
                            doublet_remove = TRUE,
                            doublet_clusters = TRUE, # vector of labels for each cell, or the name of a colData column of sce. Alternatively, if 'clusters=TRUE', fast clustering will be performed. If 'clusters' is a single integer, it will determine how many clusters to create (using k-means clustering)
                            includePCs = 40,
                            doublet_rate = "default", #0.03,  # provide "default" or NULL for 1% per thousand cells captured (so 4% among 4000 thousand cells)
                            nfeatures = 1000,
                            # seurat processing
                            run_seurat = TRUE,
                            npcs_PCA = 70, # ---> use same as doublet ?
                            # general scvi
                            run_scvi = TRUE,
                            batch_var = "Sample_ID" ,# Provide NULL
                            assay_name = "RNA",
                            integration_name = "scvi",
                            # scvi integration:
                            categorical_covariates =character(0), # batch_var#c("Dataset",batch_var)
                            continuous_covariates =character(0),
                            n_layers = 2,
                            n_latent = 50,
                            n_hidden = 256,
                            dropout_rate = 0.1,
                            max_epochs = 100, # set to 200-300
                            early_stopping = FALSE,
                            dispersion = "gene",
                            gene_likelihood = "zinb",
                            use_cuda =FALSE,
                            save_normalized = FALSE,
                            ## initial clustering
                            run_clustering = TRUE,
                            target_cluster_number = 100,
                            resolutions_to_try = c(0.5, 0.75, 1, 1.5, 2:10),
                            # annotation
                            run_annotation = TRUE,
                            celltype_signatures_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/scseq_pipeline/data/mouseHypoMap_celltype_signatures_mm.json",
                            # basic marker detection
                            run_markers = TRUE,
                            test.use = "wilcox",
                            logfc.threshold = 0.3,
                            min.pct = 0.1,
                            min.diff.pct = 0.05,
                            max.cells.per.ident = 20000,
                            min.cells.feature = 10,
                            min.cells.group =  10,
                            base = 2,
                            only.pos = TRUE
){

  ## This script creates the json with general parameters --> make otehr jsons or edit this script if other params should be used
  # requires some manually decisions which are added here.

  library(magrittr)

  parameter_list = list()

  # define these files:
  parameter_list$project_path = project_path
  parameter_list$input_seurat = input_seurat
  parameter_list$project_name =project_name

  # files that should not change:
  parameter_list$qc_path = paste0(parameter_list$project_path,"quality_control/")
  parameter_list$feature_set_file = paste0(parameter_list$project_path,"feature_set.json")
  parameter_list$genes_to_exclude_file = genes_to_exclude_file
  parameter_list$scvi_script_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/scseq_pipeline/python/process_scSeq_scvi.py"

  # processing
  parameter_list$n_cores = n_cores
  parameter_list$feature_set_size = feature_set_size
  parameter_list$id_column = id_column
  parameter_list$global_seed = global_seed
  parameter_list$minUMI = minUMI
  parameter_list$minFeatures = minFeatures
  parameter_list$maxUMI = maxUMI
  parameter_list$max_mt = max_mt
  parameter_list$min_cells_sample = min_cells_sample
  parameter_list$sample_column = sample_column
  parameter_list$k_param = k_param
  parameter_list$dist_type=dist_type

  # doublet
  parameter_list$doublet_remove = doublet_remove
  parameter_list$doublet_clusters = doublet_clusters
  parameter_list$includePCs = includePCs
  parameter_list$doublet_rate = doublet_rate
  parameter_list$nfeatures = nfeatures

  # seurat processing
  parameter_list$npcs_PCA = 70 # ---> use same as doublet ?

  # general scvi
  parameter_list$run_seurat = run_seurat
  parameter_list$run_scvi = run_scvi
  parameter_list$batch_var = batch_var
  parameter_list$feature_set_size = feature_set_size
  parameter_list$assay_name = assay_name
  parameter_list$integration_name = integration_name

  # scvi integration:
  parameter_list$categorical_covariates =character(0)
  parameter_list$continuous_covariates =character(0)
  parameter_list$n_layers = n_layers
  parameter_list$n_latent = n_latent
  parameter_list$n_hidden = n_hidden
  parameter_list$dropout_rate = dropout_rate
  parameter_list$max_epochs = max_epochs
  parameter_list$early_stopping = early_stopping
  parameter_list$dispersion = dispersion
  parameter_list$gene_likelihood = gene_likelihood
  parameter_list$use_cuda = use_cuda
  parameter_list$save_normalized = save_normalized

  ## initial clustering
  parameter_list$run_clustering = run_clustering
  parameter_list$target_cluster_number = target_cluster_number
  parameter_list$resolutions_to_try = resolutions_to_try
  parameter_list$run_annotation = run_annotation
  parameter_list$celltype_signatures_file = celltype_signatures_file

  # basic marker detection
  parameter_list$run_markers = run_markers
  parameter_list$test.use = test.use
  parameter_list$logfc.threshold = logfc.threshold
  parameter_list$min.pct = min.pct
  parameter_list$min.diff.pct = min.diff.pct
  parameter_list$max.cells.per.ident = max.cells.per.ident
  parameter_list$min.cells.feature = min.cells.feature
  parameter_list$min.cells.group =  min.cells.group
  parameter_list$base = base
  parameter_list$only.pos = only.pos

  # save
  if(!is.null(param_file)){
    dir.create(parameter_list$project_path,showWarnings = FALSE)
    message("Writing params to: ",paste0(parameter_list$project_path,param_file))
    writeList_to_JSON(parameter_list,filename = paste0(parameter_list$project_path,param_file))
    # message("Saving to: ", paste0(parameter_list$project_path,"parameters_processing.json"))
    #
  }
  # return
  return(parameter_list)
}

##########
### scseq_pipeline
##########

#' General single cell pipeline using a seurat as input
#'
#' Runs filtering, doublet detectiom, scvi, seurat, clustering, annotation and marker detection
#'
#' @param parameters a json or a list with all required parameters
#' @param verbose enable messages at each step
#' @param exec_singularity run scvi via a singularity exec call or directly via python3 (if function itself is already called from singularity)
#' @param singularity_image path to an image to use if exec_singularity is TRUE
#'
#' @return list of parameters for scseq_pipeline
#'

scseq_pipeline <- function(
  parameters,
  verbose = TRUE,
  exec_singularity = FALSE,
  singularity_image = "~/Documents/r_scvi_v3_42.simg"
){

  ##########
  ### Load parameters and packages
  ##########

  if(verbose) message(" ++++++++++ ",Sys.time(),": Load parameters and packages ")

  library(magrittr)
  library(scUtils)
  library(dplyr)
  library(Matrix)
  library(Seurat)

  if(length(parameters) == 1){
    if(file.exists(parameters)){
      # read all parameters and filepaths
      parameter_list = jsonlite::read_json(parameters)
      # if some fields are lists --> unlist
      parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
      # Set parameter_file to the file that was provided
      parameter_file <- parameters
    } else {
      stop("Cannot find parameters file. Either provide the full paramtere list or a valid json file.")
    }
  }else{
    # Parameters is provided as a list
    parameter_list <- parameters
    # Write out parameters as json to be able to pass it to python later on
    parameters_file <- tempfile(
      pattern = "scseq_pipeline_config",
      fileext = ".json"
    )
    scUtils::writeList_to_JSON(
      parameter_list,
      filename = parameters_file
    )
    # TODO: validate parameter file
  }
  # read features to excludes
  features_exclude_list= jsonlite::read_json(parameter_list$genes_to_exclude_file)
  features_exclude_list = unlist(lapply(features_exclude_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}}))

  system(paste0("mkdir -p ",parameter_list$project_path))

  ##########
  ### Load raw file
  ##########

  if(verbose) message(" ++++++++++ ",Sys.time(),": Load seurat object ")

  # load seurat
  seurat_object = readRDS(paste0(parameter_list$input_seurat))

  report_list = list()

  ##########
  ### execute QC and doublet substeps
  ##########

  # inital numbers
  report_list$project_name = parameter_list$project.name
  report_list$cell_number_raw = ncol(seurat_object)
  # filtering
  if(verbose) message(" ++++++++++ ",Sys.time(),": Filter seurat object ")
  seurat_object = filter_object(seurat_object,parameter_list,features_exclude_list)
  report_list$cell_number_filtering = ncol(seurat_object)
  # doublets
  if(verbose) message(" ++++++++++ ",Sys.time(),": Doublet detection ")
  seurat_object = detect_doublets(seurat_object,parameter_list,features_exclude_list)
  report_list$cell_number_doublets = ncol(seurat_object)
  # normalize data
  if(verbose) message(" ++++++++++ ",Sys.time(),": LogNormalize data ")
  seurat_object <- Seurat::NormalizeData(object = seurat_object,  verbose = F, assay = "RNA")
  # feature detection
  if(verbose) message(" ++++++++++ ",Sys.time(),": variable feature detection ")
  feature_set = scUtils::identify_variable_features(seurat_object,
                                                    n_hvgs_sizes = parameter_list$feature_set_size,
                                                    batch_var = parameter_list$batch_var,
                                                    assay_name = "RNA",
                                                    method = "vst",
                                                    ignore_genes_vector = features_exclude_list,
                                                    returnSeurat = FALSE,
                                                    seed = parameter_list$global_seed)
  seurat_object@assays$RNA@var.features = as.character(feature_set)
  scUtils::writeList_to_JSON(feature_set,filename = paste0(parameter_list$feature_set_file))
  # clean object
  message(" ++++++++++ ",Sys.time(),": Clean seurat object ")
  seurat_object@project.name = parameter_list$project_name
  seurat_object@misc= list()
  seurat_object@graphs= list()
  seurat_object@reductions= list()
  dummy=matrix(data = as.numeric())
  seurat_object@assays[["RNA"]]@scale.data <- dummy[,-1] # error is okay

  ##########
  ### scvi
  ##########

  # optionally run scvi:
  if(parameter_list$run_scvi){
    # export to anndata
    message(" ++++++++++ ",Sys.time(),": Save object to anndata" )
    merged_file_name = paste0(parameter_list$project_path,parameter_list$project_name,"_raw")
    # save h5seurat
    SeuratDisk::SaveH5Seurat(object = seurat_object,filename = paste0(merged_file_name,".h5seurat"), overwrite = TRUE, verbose = TRUE) %>% suppressMessages()
    # save to anndata
    SeuratDisk::Convert( paste0(merged_file_name,".h5seurat"), dest =  paste0(merged_file_name,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE) %>% suppressMessages()
    system(paste0("rm ",paste0(merged_file_name,".h5seurat")))
    # run scvi
    message(" ++++++++++ ",Sys.time(),": Starting scvi execution script in python")
    python_call <- paste0("python3 -u ",parameter_list$scvi_script_file," ",parameters_file)
    message(python_call)
    if(exec_singularity){
      system(paste0("singularity exec ",singularity_image," ", python_call))
    }else{
      system(python_call)
    }
    # add embedding manually
    message(" ++++++++++ ",Sys.time(),": Adding integrated scvi embedding ")
    current_embedding = read_embedding(paste0(parameter_list$project_path,parameter_list$project_name,"_scVI_reduction.txt"),seurat_object)
    # make dim red
    dimred <- Seurat::CreateDimReducObject(
      embeddings = as.matrix(current_embedding),
      stdev = as.numeric(apply(current_embedding, 2, stats::sd)),
      assay = "RNA",
      key = parameter_list$integration_name
    )
    # add  to object
    seurat_object@reductions[[parameter_list$integration_name]] = dimred
    # run umap and save model
    message(" ++++++++++ ",Sys.time(),": Build UMAP with ",parameter_list$k_param," n.neighbors" )
    seurat_object = Seurat::RunUMAP(seurat_object,
                                    reduction = parameter_list$integration_name,
                                    seed.use= parameter_list$global_seed,
                                    dims=1:ncol(seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
                                    reduction.name=paste0("umap_",parameter_list$integration_name),
                                    reduction.key = paste0("umap_",parameter_list$integration_name),
                                    verbose=F,
                                    n.neighbors = parameter_list$k_param,
                                    return.model = TRUE) %>% suppressMessages() %>% suppressWarnings()
    # run seurat SNN annoy
    message(" ++++++++++ ",Sys.time(),": Build SNN with ",parameter_list$k_param," n.neighbors" )
    seurat_object = Seurat::FindNeighbors(seurat_object,
                                          reduction=parameter_list$integration_name,
                                          dims = 1:ncol(seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
                                          k.param = parameter_list$k_param,
                                          nn.method="annoy",
                                          annoy.metric=parameter_list$dist_type,
                                          graph.name = c(paste0("NN_",parameter_list$integration_name),paste0("SNN_",parameter_list$integration_name)), verbose=FALSE) %>% suppressMessages() %>% suppressWarnings()
  }else{
    message(" ++++++++++ ",Sys.time(),": Not running scvi" )
  }
  ##########
  ### Seurat based processing
  ##########

  if(parameter_list$run_seurat){
    message(" ++++++++++ ",Sys.time(),": Running standard seurat Scaling+PCA+UMAP")
    seurat_object <- Seurat::ScaleData(object = seurat_object,
                                       assay = "RNA",
                                       features = feature_set,
                                       verbose = F)
    seurat_object <- Seurat::RunPCA(object = seurat_object,
                                    assay = "RNA",
                                    npcs = parameter_list$npcs_PCA,
                                    reduction.name = "pca",
                                    features = feature_set,
                                    reduction.key = "pca",
                                    verbose = F,
                                    seed.use = parameter_list$global_seed) %>% suppressMessages() %>% suppressWarnings()
    seurat_object <- Seurat::RunUMAP(object = seurat_object,
                                     assay = "RNA",
                                     reduction = "pca",
                                     reduction.name = paste0("umap_", "pca"),
                                     reduction.key = paste0("umap_", "pca"),
                                     dims = 1:parameter_list$npcs_PCA,
                                     n.neighbors = parameter_list$k_param,
                                     verbose = F,
                                     seed.use = parameter_list$global_seed) %>% suppressMessages()  %>% suppressWarnings()
  }else{
    message(" ++++++++++ ",Sys.time(),": Not running seurat processing" )
  }

  # check whether to continue and which SNN to use
  if(!parameter_list$run_seurat & !parameter_list$run_scvi){
    message(" ++++++++++ ",Sys.time(),": Completed QC. Set either run_scvi or run_seurat to TRUE to continue." )
    return(seurat_object)
  }
  if(parameter_list$run_scvi){
    clustering_SNN = paste0("SNN_",parameter_list$integration_name)
  }else{
    seurat_object = Seurat::FindNeighbors(seurat_object,
                                          reduction="pca",
                                          dims = 1:ncol(seurat_object@reductions[["pca"]]@cell.embeddings),
                                          k.param = parameter_list$k_param,
                                          nn.method="annoy",
                                          annoy.metric=parameter_list$dist_type,
                                          graph.name = c(paste0("NN_","pca"),paste0("SNN_","pca")), verbose=FALSE)
    clustering_SNN = paste0("SNN_","pca")
  }

  ##########
  ### clustering
  ##########

  if(parameter_list$run_clustering){
    # find "best" cluster resolution:
    message(" ++++++++++ ",Sys.time(),": Running preliminary louvain clustering")
    seurat_object = scUtils::determine_cluster_resolution(seurat_object = seurat_object,
                                                          target_cluster_number = parameter_list$target_cluster_number,
                                                          resolutions =  parameter_list$resolutions_to_try,
                                                          min_cells = 5,
                                                          graph_name = paste0("SNN_",parameter_list$integration_name),
                                                          cluster_col_name = "preliminary_clusters",
                                                          return_seurat = TRUE,
                                                          seed = parameter_list$global_seed)

  }

  ##########
  ### Preliminary cluster annotation
  ##########

  # celltype_signatures_file
  # run_annotation

  if(parameter_list$run_annotation & parameter_list$run_clustering){

    message(" ++++++++++ ",Sys.time(),": Run signature annotation")
    celltype_signatures = jsonlite::read_json(parameter_list$celltype_signatures_file)
    celltype_signatures = sapply(celltype_signatures,unlist)

    # add signature module scores
    ncol_before = ncol(seurat_object@meta.data)
    seurat_object = Seurat::AddModuleScore(seurat_object,features = celltype_signatures,name="Signature")
    colnames(seurat_object@meta.data)[grepl("Signature",colnames(seurat_object@meta.data))] = names(celltype_signatures)
    # make table with all signature scores
    signatures_per_cell = seurat_object@meta.data %>%
      dplyr::select(c("Cell_ID", names(celltype_signatures)))
    # add to misc
    seurat_object@misc$signatures_per_cell = signatures_per_cell

    # stats per preliminary cluster
    cluster_column = "preliminary_clusters"
    signature_per_cluster = seurat_object@meta.data %>%
      dplyr::select(c(cluster_column, names(celltype_signatures)))  %>%
      tidyr::gather(key="celltype_annotation",value="score", - !!sym(cluster_column)) %>%
      dplyr::group_by(!!sym(cluster_column),celltype_annotation) %>%
      dplyr::summarise(median_score = median(score)) %>%
      dplyr::arrange(desc(median_score)) %>%
      dplyr::ungroup()
    # get top result per cluster
    signature_per_cluster_top = signature_per_cluster %>%
      dplyr::group_by(!!sym(cluster_column)) %>%
      dplyr::top_n(n = 1,wt = median_score) %>%
      dplyr::distinct(celltype_annotation,!!sym(cluster_column))
    # add to meta data
    tempmeta = dplyr::left_join(seurat_object@meta.data[,1:ncol_before],signature_per_cluster_top, by=cluster_column)
    rownames(tempmeta) = tempmeta$Cell_ID

    seurat_object@meta.data = tempmeta

  }

  ##########
  ### Preliminary marker genes
  ##########

  if(parameter_list$run_markers & parameter_list$run_clustering){

    message(" ++++++++++ ",Sys.time(),": Preliminary marker detection" )
    # get params for markers
    test.use = parameter_list$test.use
    logfc.threshold = parameter_list$logfc.threshold
    min.pct =parameter_list$min.pct
    min.diff.pct = parameter_list$min.diff.pct
    max.cells.per.ident = parameter_list$max.cells.per.ident
    min.cells.feature = parameter_list$min.cells.feature
    min.cells.group =  parameter_list$min.cells.group
    only.pos = parameter_list$only.pos
    test.use = parameter_list$test.use

    cluster_column ="preliminary_clusters"
    Idents(seurat_object) = cluster_column
    # genes to test
    if(ncol(seurat_object@assays[['RNA']]@data) > 30000){
      set.seed(parameter_list$global_seed)
      subset = sample(colnames(seurat_object@assays[['RNA']]@data),size = 30000)
      gene_expr_dataset = seurat_object@assays[['RNA']]@data[,subset]
    }else{
      gene_expr_dataset = seurat_object@assays[['RNA']]@data
    }
    gene_expr_dataset[gene_expr_dataset != 0] <- 1 # set to 1 for occ
    gene_sums = Matrix::rowSums(gene_expr_dataset)
    rm(gene_expr_dataset)
    genes_to_include = names(gene_sums)[gene_sums > min.cells.feature ]
    genes_to_include = genes_to_include[! genes_to_include %in% features_exclude_list]
    message("Testing ",length(genes_to_include)," genes as cluster markers")
    counter=1

    # marker detection
    all_clusters = unique(seurat_object@meta.data[,cluster_column])
    all_markers_list = list()
    for(current_cluster in all_clusters){
      message("   Calculating cluster ",current_cluster," (",counter,"/",length(all_clusters),")")
      counter = counter+1
      all_markers_list[[paste0("c_",as.character(current_cluster))]] <- tryCatch({
        temp_df = Seurat::FindMarkers(object = seurat_object,
                                      ident.1 = current_cluster,
                                      assay = "RNA",
                                      logfc.threshold = logfc.threshold,
                                      features = genes_to_include,
                                      slot = "data",
                                      test.use =test.use,
                                      min.pct = min.pct,
                                      min.diff.pct = min.diff.pct,
                                      max.cells.per.ident=max.cells.per.ident,
                                      min.cells.feature = min.cells.feature,
                                      min.cells.group = min.cells.group,
                                      base = 2,
                                      verbose =FALSE,
                                      only.pos = only.pos)
        temp_df$cluster = current_cluster
        temp_df$gene = rownames(temp_df)
        temp_df
      },
      error=function(cond) {
        message("Cannot calculate markers. Error:",cond)
        return(NULL)
      })
    }
    all_markers = do.call(rbind, all_markers_list)

    all_markers$cluster = stringr::str_remove(all_markers$cluster,pattern = "c_")

    all_markers_filtered = all_markers %>% as.data.frame() %>% dplyr::filter(p_val_adj < 0.05)

    seurat_object@misc$marker_genes = all_markers_filtered

  }


  seurat_object@misc$report_list = report_list
  return(seurat_object)

  # ##########
  # ### Save results
  # ##########
  #
  # message(" ++++++++++ ",Sys.time(),": Save results" )
  # # report
  # report_df = do.call(rbind,report_list) %>% as.data.frame() %>% dplyr::mutate(info = names(report_list)) %>% dplyr::select(info, value=V1)
  # data.table::fwrite(report_df,file=paste0(parameter_list$project_path,parameter_list$project_name,"_report.txt"),sep="\t")
  # # save markers as txt
  # data.table::fwrite(all_markers,file=paste0(parameter_list$project_path,parameter_list$project_name,"_marker_genes.txt"),sep="\t")
  # # save RDS
  # saveRDS(seurat_object,file = paste0(parameter_list$project_path,parameter_list$project_name,"_processed.rds"))
  # message(" ++++++++++ ",Sys.time(),": Complete." )
  #

}

##########
### hypomap_projection paramaters
##########

projection_parameters = function(param_file = NULL,
                                 project_path,
                                 input_seurat,
                                 reference_seurat = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_publication/hypoMap.rds",
                                 reference_model = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_publication/hypoMap_model/",
                                 reference_reduction = "scvi",
                                 project_name = "project",
                                 query_assay = "RNA",
                                 n_cores = 50,
                                 id_column = "Cell_ID",
                                 max_epochs = 10,
                                 global_seed = 123456
){

  ## This script creates the json with general parameters --> make otehr jsons or edit this script if other params should be used
  # requires some manually decisions which are added here.

  library(magrittr)

  parameter_list = list()

  # define these files:
  parameter_list$project_path = project_path
  parameter_list$input_seurat = input_seurat
  parameter_list$project_name =project_name

  # processing
  parameter_list$n_cores = n_cores
  parameter_list$reference_seurat = reference_seurat
  parameter_list$reference_model = reference_model
  parameter_list$reference_reduction = reference_reduction
  parameter_list$query_assay = query_assay
  parameter_list$max_epochs = max_epochs
  parameter_list$id_column = id_column
  parameter_list$global_seed = global_seed

  # save
  if(!is.null(param_file)){
    dir.create(parameter_list$project_path,showWarnings = FALSE)
    message("Writing params to: ",paste0(parameter_list$project_path,param_file))
    writeList_to_JSON(parameter_list,filename = paste0(parameter_list$project_path,param_file))
    # message("Saving to: ", paste0(parameter_list$project_path,"parameters_processing.json"))
    #
  }
  # return
  return(parameter_list)

}

##########
### hypomap_projection
##########

#' Project a single cell dataset on hypomap
#'
#' Additional wrapper around mapscvi
#'
#' @param param_file where to export json version of parameter list to
#'
#' @return seurat object with additional reductions and metadata columns
#'

hypomap_projection = function(parameters,verbose = TRUE){

  ##########
  ### Load parameters and packages
  ##########

  if(verbose) message(" ++++++++++ ",Sys.time(),": Load parameters and packages ")

  library(magrittr)
  library(scUtils)
  library(dplyr)
  library(mapscvi) # use: https://github.com/lsteuernagel/mapscvi/blob/master/R/mapping_functions.R
  library(Seurat)

  if(length(parameters) == 1){
    if(file.exists(parameters)){
      # read all parameters and filepaths
      parameter_list = jsonlite::read_json(parameters)
      # if some fields are lists --> unlist
      parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
    }else{
      stop("Cannot find parameters file. Either provide the full paramtere list or a valid json file.")
    }
  }else{
    parameter_list =  parameters
    # TODO: validate parameter file
  }

  ##########
  ### Load raw file
  ##########

  if(verbose) message(" ++++++++++ ",Sys.time(),": Load seurat objects ")

  query_seurat_object = readRDS(parameter_list$input_seurat)
  message("reference_seurat path: ", parameter_list$reference_seurat)
  reference_seurat_object = readRDS(parameter_list$reference_seurat)

  ##########
  ### Project
  ##########

  if(verbose) message(" ++++++++++ ",Sys.time(),": Map to HypoMap ")
  if(verbose) message(" Add 'Batch_ID' to query seurat if you want to run batch-aware mapping.")
  query_seurat_object = mapscvi::map_new_seurat_hypoMap(query_seurat_object = query_seurat_object,
                                                        assay = parameter_list$query_assay,
                                                        label_col = "C465_named",
                                                        max_epochs = parameter_list$max_epochs,
                                                        reference_seurat = reference_seurat_object,
                                                        model_path = parameter_list$reference_model,
                                                        reference_reduction = parameter_list$reference_reduction)
  ## rename
  query_seurat_object@meta.data$C465_predicted = query_seurat_object@meta.data$predicted
  query_seurat_object@meta.data$C465_prediction_probability = query_seurat_object@meta.data$prediction_probability

  # propgatae manually for all others
  # C286
  prediction_C286_named= mapscvi::propagate_labels_prob(neighbors_object = query_seurat_object@neighbors$query_ref_nn,
                                                        label_vec = reference_seurat_object@meta.data$C286_named,
                                                        query_seurat_object = query_seurat_object)
  query_seurat_object@meta.data$C286_named_predicted = prediction_C286_named$predicted
  query_seurat_object@meta.data$C286_named_predicted_prob = prediction_C286_named$prediction_probability
  # C185
  prediction_C185_named= mapscvi::propagate_labels_prob(neighbors_object = query_seurat_object@neighbors$query_ref_nn,
                                                        label_vec = reference_seurat_object@meta.data$C185_named,
                                                        query_seurat_object = query_seurat_object)
  query_seurat_object@meta.data$C185_named_predicted = prediction_C185_named$predicted
  query_seurat_object@meta.data$C185_named_predicted_prob = prediction_C185_named$prediction_probability
  # C66
  prediction_C66_named= mapscvi::propagate_labels_prob(neighbors_object = query_seurat_object@neighbors$query_ref_nn,
                                                       label_vec = reference_seurat_object@meta.data$C66_named,
                                                       query_seurat_object = query_seurat_object)
  query_seurat_object@meta.data$C66_named_predicted = prediction_C66_named$predicted
  query_seurat_object@meta.data$C66_named_predicted_prob = prediction_C66_named$prediction_probability
  # C25
  prediction_C25_named= mapscvi::propagate_labels_prob(neighbors_object = query_seurat_object@neighbors$query_ref_nn,
                                                       label_vec = reference_seurat_object@meta.data$C25_named,
                                                       query_seurat_object = query_seurat_object)
  query_seurat_object@meta.data$C25_named_predicted = prediction_C25_named$predicted
  query_seurat_object@meta.data$C25_named_predicted_prob = prediction_C25_named$prediction_probability
  # C7
  prediction_C7_named= mapscvi::propagate_labels_prob(neighbors_object = query_seurat_object@neighbors$query_ref_nn,
                                                      label_vec = reference_seurat_object@meta.data$C7_named,
                                                      query_seurat_object = query_seurat_object)
  query_seurat_object@meta.data$C7_named_predicted = prediction_C7_named$predicted
  query_seurat_object@meta.data$C7_named_predicted_prob = prediction_C7_named$prediction_probability


  return(query_seurat_object)
}


##########
### deg_parameters
##########

deg_parameters = function(param_file = NULL,
                          project_path,
                          input_seurat,
                          project_name = "project",
                          cluster_variable,
                          res_name_add = "", # for slurm export
                          sample_variable="Sample_ID",
                          primary_variable="Condition",
                          reference_level = "Control",
                          other_variables=character(0),
                          nCounts = "nCount_RNA",
                          min_cells=10,
                          min_pct_valid_samples = 0.5,
                          min_pct = 0.05,
                          assay="RNA",
                          padjust_method = "bonferroni",
                          numCores = 1,
                          global_seed = 123456
){

  ## This script creates the json with general parameters --> make otehr jsons or edit this script if other params should be used
  # requires some manually decisions which are added here.

  library(magrittr)

  parameter_list = list()

  # define these files:
  parameter_list$project_path = project_path
  parameter_list$input_seurat = input_seurat
  parameter_list$project_name =project_name
  parameter_list$res_name_add = res_name_add # mostly for slurm

  # processing
  parameter_list$cluster_variable=cluster_variable
  parameter_list$sample_variable=sample_variable
  parameter_list$primary_variable=primary_variable
  parameter_list$reference_level = reference_level
  parameter_list$other_variables=other_variables
  parameter_list$nCounts = nCounts
  parameter_list$min_cells=min_cells
  parameter_list$min_pct_valid_samples = min_pct_valid_samples
  parameter_list$min_pct = min_pct
  parameter_list$assay=assay
  parameter_list$padjust_method=padjust_method
  parameter_list$numCores = numCores

  # save
  if(!is.null(param_file)){
    dir.create(parameter_list$project_path,showWarnings = FALSE)
    message("Writing params to: ",paste0(parameter_list$project_path,param_file))
    writeList_to_JSON(parameter_list,filename = paste0(parameter_list$project_path,param_file))
    # message("Saving to: ", paste0(parameter_list$project_path,"parameters_processing.json"))
    #
  }
  # return
  return(parameter_list)

}

##########
### differential_gene_expression
##########

#' Run differential gene expressin between conditions
#'
#' Requires a processed object with clustering.
#'
#' @param param_file where to export json version of parameter list to
#'
#' @return dataframe with differential gene expression results
#'
#'

differential_gene_expression = function(parameters, verbose =TRUE){

  library(magrittr)
  library(scUtils)
  library(dplyr)
  library(Seurat)

  if(length(parameters) == 1){
    if(file.exists(parameters)){
      # read all parameters and filepaths
      parameter_list = jsonlite::read_json(parameters)
      # if some fields are lists --> unlist
      parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
    }else{
      stop("Cannot find parameters file. Either provide the full parameter list or a valid json file.")
    }
  }else{
    parameter_list =  parameters
    # TODO: validate parameter file
  }

  ##########
  ### Load raw file
  ##########

  if(verbose) message(" ++++++++++ ",Sys.time(),": Load seurat object ")

  seurat_object = readRDS(parameter_list$input_seurat)

  ##########
  ### run DEG with nebula
  ##########

  if(verbose) message(" ++++++++++ ",Sys.time(),": Run DEG detection with nebula ")

  print(scUtils::FindDEG_nebula)

  print("----")

  print(FindDEG_nebula)

  nebula_result = FindDEG_nebula(seurat_object,
                                 cluster_variable= parameter_list$cluster_variable,
                                 sample_variable= parameter_list$sample_variable,
                                 primary_variable=parameter_list$primary_variable,
                                 reference_level = parameter_list$reference_level,
                                 other_variables= parameter_list$other_variables,
                                 nCounts = parameter_list$nCounts,
                                 min_cells= parameter_list$min_cells,
                                 min_pct_valid_samples = parameter_list$min_pct_valid_samples,
                                 min_pct = parameter_list$min_pct,
                                 assay= parameter_list$assay,
                                 return_full_result = TRUE,
                                 numCores = parameter_list$numCores,
                                 verbose =TRUE)

  print((nebula_result$nebula_result_list[[1]]))
  # format results
  nebula_result_df = nebula_result$degs

  print(nebula_result_df)
  print(head(nebula_result_df))

  # nebula_result
  # nebula_result_df$p_val_adj = p.adjust(
  #   p = nebula_result_df[,which(grepl("p_",colnames(nebula_result_df)))[1]],
  #   method = parameter_list$padjust_method,
  #   n = length(unique(nebula_result_df$gene))
  # )

  # return
  return(nebula_result_df)
}



