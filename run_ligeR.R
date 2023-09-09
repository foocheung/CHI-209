library(Signac)
library(Seurat)
library(patchwork)
library(rliger)
library(GenomeInfoDb) 
library("biovizBase")
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(docopt)

# Define docopt usage
doc <- "
Usage:
  script.R [--metadata=METADATA] [--counts=COUNTS] [--fragments=FRAGMENTS] [--tss=TSS] [--pbmc=PBMC]

Options:
  -h --help               Show this help message and exit.
  --metadata=METADATA     Path to metadata file.
  --counts=COUNTS         Path to counts file in h5 format.
  --fragments=FRAGMENTS   Path to fragments file in gz format.
  --tss=TSS               Threshold for TSS enrichment [default: 2].
  --pbmc=PBMC             Path to pbmc RNA file.
"


#Rscript ./run_ligeR.R \
#--metadata /Volumes/CHI/TEMP/230329_VH00286_21_AAAJWCKHV/TEST/outs/fastq_path/AAAJWCKHV/aggr_local/outs/singlecell.csv \
#--counts /Volumes/CHI/TEMP/230329_VH00286_21_AAAJWCKHV/TEST/outs/fastq_path/AAAJWCKHV/aggr_local/outs/filtered_peak_bc_matrix.h5 \
#--fragments /Volumes/CHI/TEMP/230329_VH00286_21_AAAJWCKHV/TEST/outs/fastq_path/AAAJWCKHV/aggr_local/outs/fragments.tsv.gz \
#--tss 2.8  \
#--pbmc ~/OneDrive\ -\ National\ Institutes\ of\ Health/MILO_GALINA_TEST_RUN/CHI_pbmc.rds

# Parse command line arguments
args <- docopt(doc)

metadata <- read.csv(args$'--metadata', header = TRUE, row.names = 1)
counts <- Read10X_h5(filename = args$'--counts')
fragments <- args$'--fragments'
tss_threshold <- as.numeric(args$'--tss')
pbmc_rna <- readRDS(args$'--pbmc')

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = fragments,
  min.cells = 10,
  min.features = 200
)

bm <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(bm) <- annotations

bm <- NucleosomeSignal(object = bm)
bm <- TSSEnrichment(object = bm, fast = FALSE)

bm$pct_reads_in_peaks <- bm$peak_region_fragments / bm$passed_filters * 100
bm$blacklist_ratio <- bm$blacklist_region_fragments / bm$peak_region_fragments

DensityScatter(bm, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

bm$high.tss <- ifelse(bm$TSS.enrichment > tss_threshold, paste0('High (TSS > ', tss_threshold, ')'), paste0('Low (TSS <= ', tss_threshold, ')'))
TSSPlot(bm, group.by = 'high.tss') + NoLegend()

bm$nucleosome_group <- ifelse(bm$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = bm, group.by = 'nucleosome_group')

bm_orig <- bm

bm <- subset(
  x = bm,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 30000 &
    pct_reads_in_peaks > 15 &
    nucleosome_signal < 4 &
    TSS.enrichment > tss_threshold 
)

bm <- RunTFIDF(bm)
bm <- FindTopFeatures(bm, min.cutoff = 'q0')
bm <- RunSVD(bm)
DepthCor(bm)

bm <- RunUMAP(object = bm, reduction = 'lsi', dims = 2:30)
bm <- FindNeighbors(object = bm, reduction = 'lsi', dims = 2:30)
bm <- FindClusters(object = bm, verbose = FALSE, algorithm = 3)
DimPlot(object = bm, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(bm)

bm[['RNA']] <- CreateAssayObject(counts = gene.activities)
bm <- NormalizeData(
  object = bm,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(bm$nCount_RNA)
)

DefaultAssay(bm) <- 'RNA'

FeaturePlot(
  object = bm,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

pbmc_rna <- readRDS(args$'--pbmc')
bm.rna <- pbmc_rna

data <- list(atac = as.matrix(bm@assays$RNA@counts), rna = as.matrix(bm.rna@assays$RNA@counts))
int <- createLiger(data)

int <- rliger::normalize(int)
int <- selectGenes(int)
int <- scaleNotCenter(int)
int <- optimizeALS(int, k = 20)

int <- quantile_norm(int)
int <- louvainCluster(int, resolution = 0.3, verbose = TRUE)
int <- runUMAP(int, distance = 'cosine', n_neighbors = 10, min_dist = 0.3)

plots <- plotByDatasetAndCluster(int, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE, pt.size = 0.5)
plots[[1]] + plots[[2]]
png(paste(tss_threshold, "_liger_umaps.png", sep=""), height = 1300, width = 3500, res = 300)
plots[[1]] + plots[[2]]
dev.off()

int2 <- ligerToSeurat(int, by.dataset = TRUE)
table(int2@active.ident)

int2$liger_clusters <- int2@active.ident
png(paste(tss_threshold, "_liger_split.png", sep="") , height = 1300, width = 3500, res = 300)
DimPlot(int2, group.by = "liger_clusters", label = TRUE, repel = TRUE, split.by = "orig.ident") + NoLegend()
dev.off()
