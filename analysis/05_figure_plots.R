## 05: Analysis and plots for paper

# libs
library(ggpubr)
library(tidyverse)
library(Seurat)

# output
figure_path = "paper_figures/"
source_path = "paper_sourcedata/"
dir.create(figure_path)
dir.create(source_path)

## project settings
raw_file_path = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/2023-03-almu-snseq-analysis/AgrpPomc/","almu_2023_04_AgrpPomc_raw.rds")
param_file_name = "agrpPomcNTS_params.json"
project_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/2023-03-almu-snseq-analysis/AgrpPomc/processed/"

## extra colors:
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)

# load seurat
AgrpPomcNTS_withIEG = readRDS(paste0(project_path,"AgrpPomcNTS_withIEG.rds"))

# load table
activation_result_short = data.table::fread(file = paste0(project_path,"AgrpPomcNTS_activatedClusters.txt"),data.table = F)
activation_result_short = activation_result_short[activation_result_short$cluster_annotation!="19: low-quality-neurons",]

# remove low qual neurons
AgrpPomcNTS_withIEG = subset(AgrpPomcNTS_withIEG,subset = celltype_annotation != "low-quality-neurons")

### load markers
marker_genes = AgrpPomcNTS_withIEG@misc$marker_genes %>% dplyr::filter(p_val_adj < 0.0001) %>%
  dplyr::mutate(specificity = avg_log2FC * ((pct.1+0.01) / (pct.2+0.01)))

## filter top 20 marker genes and save as additional source data
marker_genes_filt = marker_genes %>% dplyr::select(cluster,gene,specificity,avg_log2FC,pct.1,pct.2,p_val_adj) %>%
  dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = specificity,n = 20) %>%
  dplyr::filter(specificity > 1) %>%
  dplyr::arrange(cluster,desc(specificity))

data.table::fwrite(marker_genes_filt,paste0(source_path,"source_figure5bcadd_markers.txt"),sep="\t")

##########
### Overview UMAP with annotations
##########

nannos = length(unique(AgrpPomcNTS_withIEG@meta.data$celltype_annotation))+1
overview_umap = DimPlot(AgrpPomcNTS_withIEG,group.by = "celltype_annotation",label=TRUE,label.size = 5,repel = TRUE,
                        cols = getOkabeItoPalette(n = nannos),raster = TRUE,raster.dpi = c(1536,1536),pt.size = 1.9)+NoLegend()+NoAxes()+ggtitle(NULL)
overview_umap

ggsave(filename = paste0(figure_path,"overview_umap.pdf"),
       plot = overview_umap, "pdf",dpi=400,width=250,height = 250,units="mm")

## source data -- I collect the source data for figure 5 to avoid replicating thousands of umap coordinates for each plot
source_data_fig5bc =overview_umap$data

##########
### top feature plots
##########

rasterize_px = 1536
seurat_pt_size = 1.4
gene_set_to_plot = c("Rbfox3","Mog","Gfap","Olig2")
p <- FeaturePlot(AgrpPomcNTS_withIEG,features = gene_set_to_plot,keep.scale = "all", combine = FALSE,raster = TRUE,order=TRUE,
                 raster.dpi = c(rasterize_px,rasterize_px),
                 pt.size = seurat_pt_size,reduction = "umap_scvi_native")
for(i in 1:length(p)) {
  legend_plot = p[[i]]
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
  #p[[i]] <- scUtils::rasterize_ggplot(p[[i]],pixel_raster = 2048,pointsize = 1.8)
}
combined_feature_plots = cowplot::plot_grid(plotlist = p,ncol = 4)
combined_feature_plots

## save
ggsave(filename = paste0(figure_path,"feature_umaps.pdf"),
       plot = combined_feature_plots, "pdf",dpi=300,width=600,height = 150,units="mm")

# Convert to a ggplot and print
leg <- ggpubr::get_legend(legend_plot)
leg <- ggpubr::as_ggplot(leg)
leg
ggsave(filename = paste0(figure_path,"feature_umaps_legend.pdf"),
       plot = leg, "pdf",dpi=100,width=20,height = 40,units="mm")

## add features to source data
source_data_fig5bc = cbind(source_data_fig5bc,FetchData(AgrpPomcNTS_withIEG,vars = gene_set_to_plot))
data.table::fwrite(source_data_fig5bc,paste0(source_path,"source_figure5bc_umap.txt"),sep="\t")

##########
### zoom in on neurons
##########

AgrpPomcNTS_withIEG_neurons = subset(AgrpPomcNTS_withIEG,subset = celltype_annotation == "Neurons")
AgrpPomcNTS_withIEG_neurons = RunUMAP(AgrpPomcNTS_withIEG_neurons,reduction = "scvi_native",seed.use = 1234,
                                      dims = 1:ncol(AgrpPomcNTS_withIEG_neurons@reductions$scvi_native@cell.embeddings),reduction.name = "umap_neurons")

#AgrpPomcNTS_withIEG_neurons = readRDS(paste0(project_path,"AgrpPomcNTS_withIEG_neurons.rds"))

# neuron cluster umap
nannos = length(unique(AgrpPomcNTS_withIEG_neurons@meta.data$preliminary_clusters))+1
overview_umapNeurons = DimPlot(AgrpPomcNTS_withIEG_neurons,group.by = "preliminary_clusters",label=TRUE,label.size = 3,repel = F,reduction="umap_neurons",
                               cols = getOkabeItoPalette(n = nannos),raster = TRUE,raster.dpi = c(1536,1536),pt.size = 1.9)+NoLegend()+NoAxes()+ggtitle(NULL)
overview_umapNeurons

ggsave(filename = paste0(figure_path,"overview_neurons_umap.pdf"),
       plot = overview_umapNeurons, "pdf",dpi=400,width=250,height = 250,units="mm")

source_data_fig5de = overview_umapNeurons$data
source_data_fig5de = cbind(source_data_fig5de, cluster_annotation = AgrpPomcNTS_withIEG_neurons@meta.data$cluster_annotation)

##########
### a neuronal Umap showing the expression the expression of Slc17a6 and Slc32a1.
##########

# plot
rasterize_px = 1536
seurat_pt_size = 1.9
gene_set_to_plot2 = c("Slc17a6","Slc32a1")
p <- FeaturePlot(AgrpPomcNTS_withIEG_neurons,features = gene_set_to_plot2,keep.scale = "all", combine = FALSE,raster = TRUE,order=TRUE,raster.dpi = c(rasterize_px,rasterize_px),reduction="umap_neurons",
                 pt.size = seurat_pt_size)
for(i in 1:length(p)) {
  legend_plot = p[[i]]
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
  #p[[i]] <- scUtils::rasterize_ggplot(p[[i]],pixel_raster = 2048,pointsize = 1.8)
}
combined_feature_plots = cowplot::plot_grid(plotlist = p,ncol = 2)
combined_feature_plots

## save
ggsave(filename = paste0(figure_path,"feature_gaba_glut_umaps.pdf"),
       plot = combined_feature_plots, "pdf",dpi=300,width=300,height = 150,units="mm")

# Convert legend to a ggplot and print
leg <- ggpubr::get_legend(legend_plot)
leg <- ggpubr::as_ggplot(leg)
leg
ggsave(filename = paste0(figure_path,"feature_gaba_glut_umaps_legend.pdf"),
       plot = leg, "pdf",dpi=100,width=20,height = 40,units="mm")

source_data_extFig10 = cbind(source_data_fig5de[,c(1,2,3,4)],FetchData(AgrpPomcNTS_withIEG_neurons,vars=gene_set_to_plot2))
data.table::fwrite(source_data_extFig10,paste0(source_path,"source_extFig10ab_umap.txt"),sep="\t")

##########
### top neurons
##########

# highlight top 5 neuron labels and plot ieg score on neurons
AgrpPomcNTS_withIEG_neurons@meta.data$ieg_score[is.na(AgrpPomcNTS_withIEG_neurons@meta.data$ieg_score)] = 0
AgrpPomcNTS_withIEG_neurons@meta.data$ieg_label = NA
AgrpPomcNTS_withIEG_neurons@meta.data$ieg_label[AgrpPomcNTS_withIEG_neurons@meta.data$ieg_score > 4.1] = AgrpPomcNTS_withIEG_neurons@meta.data$cluster_annotation[AgrpPomcNTS_withIEG_neurons@meta.data$ieg_score > 4.1]
Idents(AgrpPomcNTS_withIEG_neurons) = "ieg_label"
ieg_scores_umap = FeaturePlot(AgrpPomcNTS_withIEG_neurons,features = "ieg_score",label=TRUE,label.size = 4,repel = TRUE,reduction="umap_neurons",
                              raster = TRUE,raster.dpi = c(1536,1536),pt.size = 1.9)+NoAxes()+
  scale_color_gradient2(low="darkblue",mid="grey80",high = "darkred")
ieg_scores_umap

ggsave(filename = paste0(figure_path,"ieg_scores_neurons_umap.pdf"),
       plot = ieg_scores_umap, "pdf",dpi=400,width=250,height = 250,units="mm")

source_data_fig5de = cbind(source_data_fig5de,ieg_score = AgrpPomcNTS_withIEG_neurons@meta.data$ieg_score)
source_data_fig5de = cbind(source_data_fig5de, ieg_category = AgrpPomcNTS_withIEG_neurons@meta.data$ieg_category)
source_data_fig5de = cbind(source_data_fig5de, FetchData(AgrpPomcNTS_withIEG_neurons,vars="Th"))
data.table::fwrite(source_data_fig5de,paste0(source_path,"source_figure5de_umap.txt"),sep="\t")

##########
###  Violin plot of the main markers of cluster 82.Tal1/Nell
##########

top_neurons = unique(AgrpPomcNTS_withIEG_neurons@meta.data$ieg_label[AgrpPomcNTS_withIEG_neurons@meta.data$ieg_score > 4.1])
top_neurons

# get top markers
nmarkers = 10
top_n_markers_per_cluster =  marker_genes %>%
  dplyr::filter(cluster %in% "82") %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(specificity > 2 & p_val_adj < 0.01) %>%
  dplyr::filter(!grepl("Rik|Gm|D7",gene)) %>%
  dplyr::filter(gene!="Th") %>%
  dplyr::slice_max(order_by = specificity,n = nmarkers,with_ties = F)

top_n_markers_per_cluster
top_n_markers_per_cluster_vec = unique(c(top_n_markers_per_cluster$gene))

## make violin
Idents(AgrpPomcNTS_withIEG_neurons) = "cluster_annotation"
nell1_vln_plot = VlnPlot(AgrpPomcNTS_withIEG_neurons,features = top_n_markers_per_cluster_vec,idents = "82: Tal1/Nell1",same.y.lims = F,stack = TRUE,cols = getOkabeItoPalette(length(top_n_markers_per_cluster_vec)))+NoLegend()
nell1_vln_plot = nell1_vln_plot+theme(text=element_text(size=20), axis.text.x = element_text(size=15),axis.title.y = element_blank())
nell1_vln_plot#

ggsave(filename = paste0(figure_path,"tal1_82_markers_Violin.pdf"),
       plot = nell1_vln_plot, "pdf",dpi=400,width=250,height = 150,units="mm")

cluster_cells = AgrpPomcNTS_withIEG_neurons@meta.data$Cell_ID[AgrpPomcNTS_withIEG_neurons@meta.data$cluster_annotation == "82: Tal1/Nell1"]
source_data_extFig10d = cbind(Cell_ID = colnames(AgrpPomcNTS_withIEG_neurons) , FetchData(AgrpPomcNTS_withIEG_neurons,vars=top_n_markers_per_cluster_vec)) %>%
  dplyr::filter(Cell_ID %in% cluster_cells)
data.table::fwrite(source_data_extFig10d,paste0(source_path,"source_extFig10d_vln.txt"),sep="\t")

all_markers_per_cluster_82 =  marker_genes %>%
  dplyr::filter(cluster %in% "82") %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(specificity > 2 & p_val_adj < 0.01) %>%
  dplyr::arrange(desc(specificity))

data.table::fwrite(all_markers_per_cluster_82,file = paste0(source_path,"source_extFig10dadd_tal1_82_markers.txt"),sep="\t")

##########
### iegs violin
##########

library(ggh4x)

top_neurons = unique(AgrpPomcNTS_withIEG_neurons@meta.data$ieg_label[AgrpPomcNTS_withIEG_neurons@meta.data$ieg_score > 4.1])
genes_to_show = unique(c(unlist(sapply(activation_result_short$IEGs[activation_result_short$cluster_annotation %in% top_neurons],strsplit,split="/")),"Fos"))
top_neurons_plots =list()
source_data_list = list()
for(n in top_neurons){
  print(n)
  cluster_to_show = n
  significant_genes = unique(unlist(sapply(activation_result_short$IEGs[activation_result_short$cluster_annotation %in% cluster_to_show],strsplit,split="/")))
  
  Idents(AgrpPomcNTS_withIEG_neurons) = "cluster_annotation"
  plist = VlnPlot(object = AgrpPomcNTS_withIEG_neurons,idents = cluster_to_show,features = genes_to_show,split.by = "Condition",ncol = length(genes_to_show),combine = F)
  for(i in 1:length(plist)){
    # save source data as well
    source_data = plist[[i]]$data  
    source_data$Cell_ID = rownames(source_data)
    source_data$gene = colnames(source_data)[1]
    colnames(source_data)[1] = "expression"
    source_data_list[[paste0(n,"_",i)]] = source_data
    # modify last plot
    if(i==length(plist)){
      plegend = plist[[i]]
      plist[[i]]$data$cluster = cluster_to_show
      plist[[i]] = plist[[i]]+ggh4x::facet_grid2(cluster ~ .,scales = "free",
                                                 strip = ggh4x::strip_themed(
                                                   background_y =  elem_list_rect(fill = c("#db6f6b")),
                                                   text_y = elem_list_text(size=c(15)),
                                                   by_layer_y = TRUE,
                                                 )) #facet_grid(cluster ~ .,scales = "free")
      
    }
    # add star if significant
    if(plist[[i]]$labels$title %in% significant_genes){plist[[i]]$labels$title = paste0(plist[[i]]$labels$title,"*")}
    # update layout for all plots
    plist[[i]] = plist[[i]]+theme(text=element_text(size=15), axis.text.x = element_blank(),axis.title.x = element_blank())+NoLegend()
  }
  cowplot::plot_grid(plotlist = plist,nrow = 1)
  top_neurons_plots[[n]] = cowplot::plot_grid(plotlist = plist,nrow = 1)
}

fullPlot = cowplot::plot_grid(plotlist = top_neurons_plots,ncol = 1)
fullPlot

ggsave(filename = paste0(figure_path,"IEG_violin.pdf"),
       plot = fullPlot, "pdf",dpi=300,width=430,height = 300,units="mm")

# Convert legend to a ggplot and print
leg <- ggpubr::get_legend(plegend)
leg <- ggpubr::as_ggplot(leg)
leg
ggsave(filename = paste0(figure_path,"IEG_violin_legend.pdf"),
       plot = leg, "pdf",dpi=100,width=20,height = 40,units="mm")

# reformat source data
source_data_extFig10c = do.call(rbind,source_data_list)
source_data_extFig10c$expression[source_data_extFig10c$expression < 0.001] = 0 # the violin plot has very small non-zero values for the not expressed genes
source_data_extFig10c = source_data_extFig10c %>% tidyr::spread(key = "gene",value="expression")

# and save
data.table::fwrite(source_data_extFig10c,paste0(source_path,"source_extFig10c_IEGvln.txt"),sep="\t")


##########
###  th cluster violin
##########

## th violin
th_clusters = marker_genes$cluster[marker_genes$gene=="Th" & marker_genes$specificity>5]

## th umap
AgrpPomcNTS_withIEG_neurons@meta.data$label = NA
AgrpPomcNTS_withIEG_neurons@meta.data$label[AgrpPomcNTS_withIEG_neurons@meta.data$preliminary_clusters %in% th_clusters] = AgrpPomcNTS_withIEG_neurons@meta.data$cluster_annotation[AgrpPomcNTS_withIEG_neurons@meta.data$preliminary_clusters %in% th_clusters]
Idents(AgrpPomcNTS_withIEG_neurons) = "label"
th_neuron_umap = FeaturePlot(AgrpPomcNTS_withIEG_neurons,features = "Th",label=TRUE,order=TRUE,label.size = 4,repel = TRUE,reduction="umap_neurons",
                             raster = TRUE,raster.dpi = c(1536,1536),pt.size = 1.9)+NoAxes()
th_neuron_umap
ggsave(filename = paste0(figure_path,"th_expression_neurons_umap.pdf"),
       plot = th_neuron_umap, "pdf",dpi=400,width=250,height = 250,units="mm")

# get top markers
nmarkers = 3
top_n_markers_per_cluster =  marker_genes %>%
  dplyr::filter(cluster %in% th_clusters) %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(specificity > 2 & p_val_adj < 0.01) %>%
  dplyr::filter(!grepl("Rik|Gm|D7",gene)) %>%
  dplyr::filter(gene!="Th") %>%
  dplyr::slice_max(order_by = specificity,n = nmarkers,with_ties = F)

top_n_markers_per_cluster_vec = unique(c("Th",top_n_markers_per_cluster$gene,"Calcr","Prlh"))

## make violin

th_clusters_names = unique(AgrpPomcNTS_withIEG_neurons@meta.data$cluster_annotation[AgrpPomcNTS_withIEG_neurons@meta.data$preliminary_clusters %in% th_clusters])
manual_order = c("63: Ebf2/Th","146: Rxfp1/Ebf2","76: Dbh/Phox2b","97: Gfral/Sctr") # or th_clusters
AgrpPomcNTS_withIEG_neurons@meta.data$label2 = NA
AgrpPomcNTS_withIEG_neurons@meta.data$label2[AgrpPomcNTS_withIEG_neurons@meta.data$cluster_annotation %in% manual_order] = AgrpPomcNTS_withIEG_neurons@meta.data$cluster_annotation[AgrpPomcNTS_withIEG_neurons@meta.data$cluster_annotation %in% manual_order]
AgrpPomcNTS_withIEG_neurons@meta.data$label2 = factor(AgrpPomcNTS_withIEG_neurons@meta.data$label2,levels = rev(c(manual_order,NA)))
Idents(AgrpPomcNTS_withIEG_neurons) = "label2"
th_vln_plot = VlnPlot(AgrpPomcNTS_withIEG_neurons,features = top_n_markers_per_cluster_vec,idents = manual_order,same.y.lims = F,stack = TRUE,cols = getOkabeItoPalette(length(top_n_markers_per_cluster_vec)))+NoLegend()
th_vln_plot = th_vln_plot+theme(text=element_text(size=20), axis.text.x = element_text(size=15),axis.title.y = element_blank())
th_vln_plot#

ggsave(filename = paste0(figure_path,"th_calcr_expression_neurons_Violin.pdf"),
       plot = th_vln_plot, "pdf",dpi=400,width=270,height = 150,units="mm")

#source_data_fig5f = th_vln_plot$data %>% tidyr::spread(key = "feature",value="expression")
cluster_cells = AgrpPomcNTS_withIEG_neurons@meta.data$Cell_ID[AgrpPomcNTS_withIEG_neurons@meta.data$cluster_annotation %in% unique(as.character(th_vln_plot$data$ident))]
source_data_fig5f = data.frame(Cell_ID = colnames(AgrpPomcNTS_withIEG_neurons), cluster_annotation = AgrpPomcNTS_withIEG_neurons@meta.data$cluster_annotation)
source_data_fig5f = cbind(source_data_fig5f, FetchData(AgrpPomcNTS_withIEG_neurons,vars=top_n_markers_per_cluster_vec)) %>%
  dplyr::filter(Cell_ID %in% cluster_cells) %>% dplyr::arrange(cluster_annotation)
# and save
data.table::fwrite(source_data_fig5f,paste0(source_path,"source_figure5f_vln.txt"),sep="\t")

## cluster info
th_cluster_info = AgrpPomcNTS_withIEG_neurons@meta.data %>%
  dplyr::filter(preliminary_clusters %in% th_clusters) %>%
  dplyr::group_by(preliminary_clusters) %>%
  dplyr::add_count(name="ncells") %>%
  dplyr::mutate(meanUMI = mean(nCount_RNA),meanFeatures = mean(nFeature_RNA)) %>%
  dplyr::distinct(cluster=preliminary_clusters,cluster_annotation,ncells,meanUMI,meanFeatures,ieg_score)

th_cluster_info

data.table::fwrite(th_cluster_info,file = paste0(source_path,"source_figure5fadd_th_info.txt"),sep="\t")


##########
### iegs violin th clusters
##########

library(ggh4x)

top_neurons = manual_order#unique(AgrpPomcNTS_withIEG_neurons@meta.data$ieg_label[AgrpPomcNTS_withIEG_neurons@meta.data$ieg_score > 4.1])
genes_to_show = unique(c(unlist(sapply(activation_result_short$IEGs[activation_result_short$cluster_annotation %in% top_neurons],strsplit,split="/")),"Fos"))
top_neurons_plots =list()

for(n in top_neurons){
  print(n)
  cluster_to_show = n
  significant_genes = unique(unlist(sapply(activation_result_short$IEGs[activation_result_short$cluster_annotation %in% cluster_to_show],strsplit,split="/")))
  
  Idents(AgrpPomcNTS_withIEG_neurons) = "cluster_annotation"
  plist = VlnPlot(object = AgrpPomcNTS_withIEG_neurons,idents = cluster_to_show,features = genes_to_show,split.by = "Condition",ncol = length(genes_to_show),combine = F)
  for(i in 1:length(plist)){
    if(i==length(plist)){
      plegend = plist[[i]]
      plist[[i]]$data$cluster = cluster_to_show
      plist[[i]] = plist[[i]]+ggh4x::facet_grid2(cluster ~ .,scales = "free",
                                                 strip = ggh4x::strip_themed(
                                                   background_y =  ggh4x::elem_list_rect(fill = c("#db6f6b")),
                                                   text_y = ggh4x::elem_list_text(size=c(15)),
                                                   by_layer_y = TRUE,
                                                 ))
      
      #facet_grid(cluster ~ .,scales = "free")
    }
    if(plist[[i]]$labels$title %in% significant_genes){plist[[i]]$labels$title = paste0(plist[[i]]$labels$title,"*")}
    plist[[i]] = plist[[i]]+theme(text=element_text(size=15), axis.text.x = element_blank(),axis.title.x = element_blank())+NoLegend()
  }
  cowplot::plot_grid(plotlist = plist,nrow = 1)
  top_neurons_plots[[n]] = cowplot::plot_grid(plotlist = plist,nrow = 1)
}

fullPlot = cowplot::plot_grid(plotlist = top_neurons_plots,ncol = 1)
fullPlot

ggsave(filename = paste0(figure_path,"IEG_th_violin.pdf"),
       plot = fullPlot, "pdf",dpi=300,width=255,height = 300,units="mm")

# Convert to a ggplot and print
leg <- ggpubr::get_legend(plegend)
leg <- ggpubr::as_ggplot(leg)
leg
ggsave(filename = paste0(figure_path,"IEG_th_violin_legend.pdf"),
       plot = leg, "pdf",dpi=100,width=20,height = 40,units="mm")



##########
###  save neurons
##########

saveRDS(AgrpPomcNTS_withIEG_neurons,paste0(project_path,"AgrpPomcNTS_withIEG_neurons.rds"))
data.table::fwrite(activation_result_short,file = paste0(figure_path,"AgrpPomcNTS_activatedClusters.txt"),sep="\t")
