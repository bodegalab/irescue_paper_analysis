# libraries required to be installed and loaded
library(data.table)
library(Seurat)
library(SeuratWrappers)

# libraries required to be installed, but not loaded
check_installs <- function(pkgname) {
  if (!pkgname %in% installed.packages()) {
    print(paste('ERROR: package', pkgname, 'is not installed. Please install', pkgname, 'and try again.'))
    stop()
  }
}
check_installs('R.utils')
check_installs('presto')

#
# DATA LOAD
#

# function to integrate TE and Gene expression counts in one Seurat object
read_irescue <- function(gex_dir, te_dir, gex_reader = 'STARsolo',
                         te_gene_column = 1, te_cell_column = 1,
                         gex_min_cells = 3, gex_min_features = 200,
                         project = 'irescue', te_assay_name = 'TE') {
  # import gene expression counts
  if (tolower(gex_reader) %in% c('star','starsolo')) {
    gex.data <- Seurat::ReadSTARsolo(file.path(gex_dir))
  } else if (tolower(gex_reader) %in% c('cellranger','10x')) {
    gex.data <- Seurat::Read10X(file.path(gex_dir))
  } else {
    stop('Invalid value for `gex_reader` parameter.')
  }
  # import TE counts
  te.data <- Seurat::Read10X(file.path(te_dir),
                             gene.column = te_gene_column,
                             cell.column = te_cell_column)
  # make Seurat object from gene counts
  obj <- Seurat::CreateSeuratObject(gex.data,
                                    min.cells = gex_min_cells,
                                    min.features = gex_min_features,
                                    project = project)
  # add TE expression to Seurat object
  te.assay <- Seurat::CreateAssayObject(te.data)
  obj[[te_assay_name]] <- subset(te.assay, colnames(te.assay)[which(colnames(te.assay) %in% colnames(obj))])
  return(obj)
}

# function to run read_irescue() on a sample, and pre-prend the sample name to cell IDs
import_crc_counts <- function(sample, parent_dir) {
  gex_dir <- file.path(parent_dir, sample, 'STAR', paste0(sample, '.Solo.out'), 'Gene/filtered')
  te_dir  <- file.path(parent_dir, sample, 'IRESCUE')
  crc <- read_irescue(gex_dir = gex_dir, te_dir = te_dir, project = sample)
  crc <- RenameCells(crc, new.names = paste(sample, Cells(crc), sep = '_'))
  return(crc)
}

# Download cell-type annotation
anno_url <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144735/suppl/GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz'
anno_file <- basename(anno_url)
if (!file.exists(anno_file)) {
  download.file(anno_url, anno_file)
}

# Import cell annotation
anno <- fread(anno_file, data.table = FALSE)
rownames(anno) <- anno[,1]
anno[,1] <- NULL

# Path to counts data
samples_dir <- file.path('./Data/results/')

# Import all samples
samples <- unique(anno$Sample)
names(samples) <- samples
crc_list <- lapply(samples, import_crc_counts, parent_dir = samples_dir)

# merge cells into one object
crc <- crc_list[[1]]
for (i in 2:length(crc_list)) {
  crc <- merge(crc, crc_list[[i]])
}
rm(crc_list); gc()

#
# DATA ANNOTATION
#

# add annotation to meta data
crc@meta.data <- merge.data.frame(crc@meta.data, anno, by = 0, all.x = TRUE)
rownames(crc@meta.data) <- crc@meta.data$Row.names
crc@meta.data$Row.names <- NULL
crc@meta.data$Sample <- crc@meta.data$orig.ident
crc@meta.data$Patient <- sub("-[BNT]$", "", crc@meta.data$Sample)
crc@meta.data[ grepl("N$", crc@meta.data$Sample), ]$Class <- "Normal"
crc@meta.data[ grepl("B$", crc@meta.data$Sample), ]$Class <- "Border"
crc@meta.data[ grepl("T$", crc@meta.data$Sample), ]$Class <- "Tumor"

# import TE metadata from IRescue output
tetable <- unique(rbindlist(lapply(list.files(samples_dir, pattern = '^features.tsv.gz$', recursive = TRUE, full.names = TRUE), fread, header = FALSE, select = 1:2, col.names = c('subfamily', 'family'))))
tetable[, c('family','class') := tstrsplit(family, '/')]
tetable <- tetable[subfamily %in% rownames(crc[['TE']])]

#
# DATA FILTERING
#

# remove unannotated cells
crc <- subset(crc, cells = rownames(anno))

# filter for epithelial cells
crc <- subset(crc, cells = rownames(crc@meta.data[crc@meta.data$Class!="Border" & crc@meta.data$Cell_type=="Epithelial cells",]))

# remove cells with no TE
crc <- subset(crc, cells = names(which(colSums(GetAssayData(crc, assay = 'TE')) > 0)))

# assign levels to columns
crc@meta.data$Class <- factor(crc@meta.data$Class, c('Normal','Border','Tumor'))
crc@meta.data$Cell_subtype <- factor(crc@meta.data$Cell_subtype,
                                     c('Stem-like/TA','Intermediate',
                                       'Goblet cells','Mature Enterocytes',
                                       'BEST4+ Enterocytes','Tuft cells',
                                       'CMS1','CMS2','CMS3','CMS4'))

#
# DATA ANALYSIS
#

# normalize and scale TE counts
DefaultAssay(crc) <- 'TE'
VariableFeatures(crc) <- rownames(crc)
crc <- NormalizeData(crc, verbose = FALSE)
crc <- ScaleData(crc, verbose = FALSE)

# Run PCA, UMAP and clustering on TE expression
crc <- RunPCA(crc, reduction.name = 'te.pca', verbose = FALSE)
crc <- RunUMAP(crc, reduction = 'te.pca', reduction.name = 'te.umap', dims = 1:9, verbose = FALSE)
crc <- FindNeighbors(crc, dims = 1:20, reduction = 'te.pca', verbose = FALSE)
crc <- FindClusters(crc, resolution = .8, verbose = FALSE, algorithm = 1)
crc[['te.clusters']] <- crc[['seurat_clusters']]

#
# DIFFERENTIAL TE EXPRESSION
#

Idents(crc) <- 'te.clusters'
cluster_markers <- as.data.table(RunPrestoAll(crc, assay = 'TE', only.pos = TRUE, verbose = FALSE), keep.rownames = 'te')

