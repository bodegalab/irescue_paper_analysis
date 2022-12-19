# libraries required to be installed and loaded
library(data.table)
library(Seurat)
library(SeuratDisk) #remotes::install_github("mojaveazure/seurat-disk")

# libraries required to be installed, but not loaded
check_installs <- function(pkgname) {
  if (!pkgname %in% installed.packages()) {
    print(paste('ERROR: package', pkgname, 'is not installed. Please install', pkgname, 'and try again.'))
    stop()
  }
}
check_installs('R.utils')

#
# DATA LOAD
#

f_sim = './Simulations/sims_mex/'
f_star <- './Simulations/tecount/STAR/sim_pbmc8k.Solo.out/Gene/filtered/'
f_irescue <- './Simulations/tecounts/IRESCUE/IRescue_out/'
f_scte <- './Simulations/tecounts/scTE/sim_pbmc8k.h5ad'
f_rmsk <- './Simulations/annotations.bed'

# Import simulated counts
sim <- Read10X(f_sim, gene.column = 1, cell.column = 1)
sim <- CreateSeuratObject(sim, assay = 'TE', project = 'sim')
star_barcodes <- readLines(file.path(f_star,'barcodes.tsv'))
sim <- subset(sim, cells = star_barcodes)

# Import IRescue counts
irescue <- Read10X(f_irescue, gene.column = 1, cell.column = 1)
irescue <- CreateSeuratObject(irescue, assay = 'TE', project = 'irescue')

# Import scTE counts
SeuratDisk::Convert(f_scte, dest = 'h5seurat', overwrite = TRUE)
scte <- SeuratDisk::LoadH5Seurat(sub('h5ad$', 'h5seurat', f_scte))
rmsk <- unique(fread(f_rmsk, select = c(7,8,9), col.names = c('subfamily', 'family', 'class')))
scte <- scte[rownames(scte) %in% rmsk$subfamily]
scte <- RenameAssays(scte, 'RNA' = 'TE')
scte@meta.data$orig.ident <- 'scte'

#
# DATA ANNOTATION
#

# Function to annotate the percentageFeatureSet of a list of features
annotate_metadata <- function(object, features, te.assay = 'TE') {
    require(Seurat)
    for (i in seq(features)) {
        f <- unique(features[[i]][which(features[[i]] %in% rownames(object[[te.assay]]))])
        col <- names(features[i])
        object[[col]] <- Seurat::PercentageFeatureSet(object, features = f, assay = te.assay)
    }
    return(object)
}
sim <- annotate_metadata(sim, features = lapply(split(rmsk[class %in% c('LINE','SINE','LTR')], by = 'class'), '[[', 'subfamily'))
sim <- annotate_metadata(sim, features = lapply(split(rmsk[family %in% c('L1','L2','Alu','MIR','ERV1','ERVL','ERVL-MaLR')], by = 'family'), '[[', 'subfamily'))
irescue <- annotate_metadata(irescue, features = lapply(split(rmsk[class %in% c('LINE','SINE','LTR')], by = 'class'), '[[', 'subfamily'))
irescue <- annotate_metadata(irescue, features = lapply(split(rmsk[family %in% c('L1','L2','Alu','MIR','ERV1','ERVL','ERVL-MaLR')], by = 'family'), '[[', 'subfamily'))
scte <- annotate_metadata(scte, features = lapply(split(rmsk[class %in% c('LINE','SINE','LTR')], by = 'class'), '[[', 'subfamily'))
scte <- annotate_metadata(scte, features = lapply(split(rmsk[family %in% c('L1','L2','Alu','MIR','ERV1','ERVL','ERVL-MaLR')], by = 'family'), '[[', 'subfamily'))

#
# DATA ANALYSIS
#

VariableFeatures(sim) <- rownames(sim)
sim <- NormalizeData(sim, verbose = FALSE)
irescue <- NormalizeData(irescue, verbose = FALSE)
scte <- NormalizeData(scte, verbose = FALSE)
sim <- ScaleData(sim, verbose = FALSE)
irescue <- ScaleData(irescue, verbose = FALSE)
scte <- ScaleData(scte, verbose = FALSE)
sim <- RunPCA(sim, verbose = FALSE)
irescue <- RunPCA(irescue, verbose = FALSE)
scte <- RunPCA(scte, verbose = FALSE)
sim <- FindNeighbors(sim, dims = 1:10, verbose = FALSE)
irescue <- FindNeighbors(irescue, dims = 1:10, verbose = FALSE)
scte <- FindNeighbors(scte, dims = 1:10, verbose = FALSE)

# find clusters
sim <- FindClusters(sim, resolution = .7)

# find the same number of clusters as the simulated ones
irescue <- FindClusters(irescue, resolution = .7)
scte <- FindClusters(scte, resolution = .5)

# prepare data frame of clusters
clusters <- dcast(
    rbind(as.data.table(irescuen@meta.data, keep.rownames = 'cell')[, .(cell, seurat_clusters, orig.ident)],
          as.data.table(simn@meta.data, keep.rownames = 'cell')[, .(cell, seurat_clusters, orig.ident)],
          as.data.table(scten@meta.data, keep.rownames = 'cell')[, .(cell, seurat_clusters, orig.ident)]),
    cell ~ orig.ident,
    value.var = 'seurat_clusters'
)
clusters[, irescue := sub('^0$','c0',irescue)]
clusters[, irescue := sub('^1$','c1',irescue)]
clusters[, irescue := sub('^2$','c3',irescue)]
clusters[, irescue := sub('^3$','c4',irescue)]
clusters[, irescue := sub('^4$','c6',irescue)]
clusters[, irescue := sub('^5$','c5',irescue)]
clusters[, irescue := sub('^6$','c7',irescue)]
clusters[, irescue := sub('^7$','c2',irescue)]
clusters[, irescue := sub('^8$','c8',irescue)]
clusters[, scte := sub('^0$','c0',scte)]
clusters[, scte := sub('^1$','c3',scte)]
clusters[, scte := sub('^2$','c4',scte)]
clusters[, scte := sub('^3$','c5',scte)]
clusters[, scte := sub('^4$','c6',scte)]
clusters[, scte := sub('^5$','c7',scte)]
clusters[, scte := sub('^6$','c1',scte)]
clusters[, scte := sub('^7$','c2',scte)]
clusters[, scte := sub('^8$','c8',scte)]
clusters[, scte := sub('^c','',scte)]

#
# UMI stats
#

umi_table <- './Simulations/tecounts/IRESCUE/tmp/cb_umi_te.bed.gz'
umi <- fread(umi_table, header = FALSE, col.names = c('cell', 'umi', 'te'))
umi[, umi := paste(cell, umi, sep = "_")]
umi_n <- umi[, .N, .(umi,te)]

# uniquely mapped UMIs
uniq_mapped <- umi[, .N, umi][N==1, umi]

# UMIs multimapping in same subfamily
same_sf <- umi_n[! umi %chin% uniq_mapped][, .N, umi][N==1, umi]

# UMIs multimapping in different subfamily
diff_sf <- umi_n[! umi %chin% c(uniq_mapped, same_sf), umi]

umicount <- data.table(
    n = c(length(uniq_mapped),length(same_sf),length(diff_sf)),
    cat = c('Uniquely mapped','Multi-mapped intra-subfamly','Multi-mapped inter-subfamily')
)
umicount[, perc := round(n/sum(n)*100)]
umicount[, mil := n/1e6]
umicount[, cat := factor(cat, c('Uniquely mapped','Multi-mapped intra-subfamly','Multi-mapped inter-subfamily'))]
umicount[, percent := paste0(as.character(perc), "%")]

# EC count & stats
ec <- unique(umi)
ec <- ec[, .(ec_tes = length(te), ec = paste(te, collapse = "|")), .(cell, umi)]
ec[, umis_by_ec_cell := .N, .(ec,cell)]
cells_by_ec <- ec[, .(cells = uniqueN(cell)), .(ec,ec_tes)]
