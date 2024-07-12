library(SingleCellExperiment)
library(scFlow)
library(scFlowData)
library(ggplot2)
library(RColorBrewer)
library(heatmaply)
library(gplots)
library(viridis)
library(cowplot)

ensembl_fp <- system.file("extdata","ensembl_mappings.tsv",package="scFlowData")
ctd_fp <- system.file("extdata","ctd",package="scFlowData")

# Produce metadata table
df <- data.frame(
  manifest = c("d8", "d14", "d21", "d28", "d35"),
  age = c(8, 14, 21, 28, 35)
)
write.table(df, "sampsheet.tsv", sep = "\t", quote = F,
            col.names = T, row.names = F)

# Produce sample folder path
pf <- data.frame(
  key = c("d8", "d14", "d21", "d28", "d35"),
  filepath = c("~/raw_mat/d8", 
               "~/raw_mat/d14",
               "~/raw_mat/d21",
               "~/raw_mat/d28",
               "~/raw_mat/d35")
  )
write.table(df, "paths.tsv", sep = "\t", quote = F,
            col.names = T, row.names = F)

manifest <- pf
samplesheet <- df
dir_list <- manifest$filepath

# create SingleCellExoeiment objects and QCs with annotaitons
# we did not run emptydrop() due to low counts had been filterd before

outdir <- "~/QC"
j <- 1
for (i in dir_list) {
  mat <- read_sparse_matrix(i)
  metadata <- read_metadata(
    unique_key = manifest$key[j],
    key_colname = "manifest",
    samplesheet_path = "sampsheet.tsv",
    col_classes = list("age" = "factor")
  )
  sce <- generate_sce(mat, metadata)
  sce <- annotate_sce(
    sce,
    min_library_size = 500,
    max_library_size = "adaptive",
    min_features = 200,
    max_features = "adaptive",
    max_mito = 0.05,
    min_ribo = 0,
    max_ribo = 1,
    min_cells = 3,
    drop_mito = TRUE,
    drop_ribo = FALSE,
    ensembl_mapping_file = ensembl_fp
  )
  sce <- filter_sce(sce)
  sce <- find_singlets(sce, "doubletfinder", pK = 0.005,
                       vars_to_regress_out = c("nCount_RNA", "pc_mito"),
                       num.cores = 1)
  sce <- filter_sce(sce)
  report_qc_sce(
    sce = sce,
    report_folder_path = outdir,
    report_file = paste0(manifest$key[j],"_qc_report")
  ) 
  write_sce(sce = sce,
                folder_path = file.path("~/QC/sce_individual",
                                        paste("sce", basename(i), sep = "_")),
                overwrite = TRUE)
  j <- j+1
}

# Read all sce individuals together

dir_sce <- "~/QC/sce_individual"

sce_path <- dir(
  path = dir_sce,
  pattern = "sce_",
  full.names = TRUE
)
sce_pathlist <- list()

for (i in sce_path) {
  sce_pathlist[[i]] <- i
}

sce_list <- lapply(sce_pathlist, read_sce)

# extract barcodes of each cells as a whitelist file for umi_tools
for (i in sce_path) {
  write.table(
    sce_list[[i]]@colData@rownames,
    file = paste0("~/whitelist/",
                  sce_list[[i]]@colData@listData[["manifest"]][1],
                  "_whitelist.txt"),
    quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
  )
}

# merge of post-QC samples
sce <- merge_sce(
  sce_list,
  ensembl_mapping_file = ensembl_fp
)

# generate merged report
sce <- annotate_merged_sce(
  sce,
  plot_vars = c("total_features_by_counts","total_counts","pc_mito","pc_ribo"),
  unique_id_var = "manifest",
  facet_vars = "age",
  outlier_vars = c("total_features_by_counts", "total_counts")
)

report_merged_sce(sce)

write_sce(sce,
          "~/sce_merged",
          write_metadata = T)

# Integration with LIGER

sce <- read_sce("~/sce_merged", read_metadata = T)
sce <- integrate_sce(
  sce,
  method             = "Liger",
  unique_id_var      = "manifest",
  take_gene_union    = FALSE,
  remove_missing     = TRUE,
  num_genes          = 3000,
  combine            = "union",
  keep_unique        = FALSE,
  capitalize         = FALSE,
  use_cols           = TRUE,
  k                  = 30,
  lambda             = 5.0,
  thresh             = 0.0001,
  max_iters          = 100,
  nrep               = 1,
  rand_seed          = 1,
  knn_k              = 20,
  k2                 = 500,
  prune_thresh       = 0.2,
  min_cells          = 2,
  quantiles          = 50,
  nstart             = 10,
  resolution         = 1,
  dims_use           = NULL,
  dist_use           = "CR",
  center             = FALSE,
  small_clust_thresh = 0
)
write_sce(
  sce,
  file.path("~/sce_integrated"),
  write_metadata = TRUE
)

sce@colData@listData[["total_features_by_counts"]] <-
  as.double(sce@colData@listData[["total_features_by_counts"]])

# reducu dimention 

sce <- reduce_dims_sce(sce,
                       pca_dims = 20)

png("~/plots/umap_manifest.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) + 
  scale_colour_manual(
    values = c("#FFA500", "#7CCD7C", "#87CEFA", "#AB82FF", "#EE0000"))
dev.off()

set.seed(321)
sce <- cluster_sce(sce,
                   cluster_method = "leiden",
                   reduction_method = "UMAP_Liger",
                   pca_dims = 20,
                   res = 0.001,
                   k = 100
)

write_sce(
  sce,
  "~/sce_clustered", write_metadata = T
)

sce <- annotate_integrated_sce(
  sce,
  categorical_covariates = c("manifest", "age"),
  input_reduced_dim = "UMAP"
)

report_integrated_sce(sce)

# automatic celltyping

set.seed(123)
sce <- map_celltypes_sce(
  sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

p <- plot_reduced_dim(sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 0.8,
                      alpha = 1)
p <- p + scale_colour_manual(values = rainbow(12)) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())

png("~/plots/umap_celltype.png",
    width = 170, height = 130, units = "mm", res = 600)
p
dev.off()

sce <- annotate_celltype_metrics(
  sce,
  unique_id_var = "manifest",
  cluster_var = "clusters",
  celltype_var = "cluster_celltype",
  facet_vars = c("manifest", "age"),
  metric_vars = c("pc_mito","pc_ribo","total_counts","total_features_by_counts")
)

report_celltype_metrics(sce)
png("~/plots/MALAT1_umap.png",
    width = 100, height = 100, units = "mm", res = 600)

plot_reduced_dim_gene(
  sce,
  reduced_dim = "UMAP_Liger",
  gene = "MALAT1",
  size = 0.3,
  alpha = 1,
  palette = c("grey80", "#CD1076")
)

dev.off()
