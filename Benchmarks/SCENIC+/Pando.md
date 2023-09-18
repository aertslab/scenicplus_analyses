Installation

```bash

conda create -n Pando_env

conda activate Pando_env

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install -c bioconda r-signac

conda install -c bioconda bioconductor-bsgenome.hsapiens.ucsc.hg38

conda install -c bioconda bioconductor-ensdb.hsapiens.v86

conda install -c conda-forge r-devtools

conda install -c bioconda bioconductor-biovizbase

conda install -c r r-feather

conda install -c conda-forge r-dosnow

conda install -c r r-doparallel

conda install -c conda-forge r-rmpi

```

```r

devtools::install_github('quadbiolab/Pando')

```

```bash

conda activate Pando_env

```

# Preprocessing using Signac and Seurat

```r

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Matrix)
library(feather)
set.seed(1234)

# load RNA and ATAC data
GEX_counts_fname <- "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/data/count_matrices/scRNA_count_matrix.tsv"
GEX_counts <- read.table(
		file = GEX_counts_fname,
		header = TRUE,
		sep = '\t',
		quote = "",
		row.names = 1)
GEX_counts <- Matrix(data.matrix(GEX_counts), sparse = TRUE)

ACC_counts_fname <- "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/data/count_matrices/sc_fragment_matrix.feather"
ACC_counts <- read_feather(ACC_counts_fname)
ACC_counts <- as.matrix(ACC_counts)
rownames(ACC_counts) <- ACC_counts[, "index"]
ACC_counts <- ACC_counts[, -1]
ACC_counts <- Matrix(ACC_counts, sparse = TRUE)
colnames(ACC_counts) <- unlist(lapply(colnames(ACC_counts), function(x) gsub('___DPLC', '', x)))

fragpath <- "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/data/atac_fragments/fragments.sorted.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
seurat_obj <- CreateSeuratObject(
		counts = GEX_counts,
		assay = "RNA")

# create ATAC assay and add it to the object
seurat_obj[["ATAC"]] <- CreateChromatinAssay(
	  counts = ACC_counts,
	  sep = c(":", "-"),
	  fragments = fragpath,
	  annotation = annotation
	)

# QC stats

DefaultAssay(seurat_obj) <- "ATAC"

seurat_obj <- NucleosomeSignal(seurat_obj)
seurat_obj <- TSSEnrichment(seurat_obj)

plot_dir <- "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/Pando/plots"

pdf(paste0(plot_dir, '/atac_RNA_QC.pdf'), width = 20)
VlnPlot(
  object = seurat_obj,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
dev.off()


DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj)

DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 5)
seurat_obj <- RunTFIDF(seurat_obj)
seurat_obj <- RunSVD(seurat_obj)

out_dir <- "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/Pando/outs"

saveRDS(seurat_obj, file = paste0(out_dir, '/seurat_obj.rds'))

seurat_obj <- readRDS(paste0(out_dir, '/seurat_obj.rds'))

# Subset out special chromosomes, otherwise an error occurs due to: StringToGRanges(rownames(GetAssay(object, assay=peak_assay)))
chromosome_length <- read.table('/staging/leuven/stg_00002/lcb/resources/human/hg38/hg38.chrom.sizes')
chromosome_length <- chromosome_length[unlist(lapply(chromosome_length$V1, function(x) !grepl('_', x))), ]

regions_df <- data.frame(seurat_obj@assays$ATAC@ranges)

region_names_to_keep <- do.call(paste, 
	c(regions_df[regions_df$seqnames %in% unlist(list(chromosome_length$V1)), c('seqnames', 'start', 'end')], sep = "-"))

gene_names <- row.names(seurat_obj@assays$RNA@meta.features)

seurat_obj <- subset(x = seurat_obj, features = c(region_names_to_keep, gene_names))


library(Pando)

data(SCREEN.ccRE.UCSC.hg38)

seurat_obj <- initiate_grn(
		object = seurat_obj,
		regions = SCREEN.ccRE.UCSC.hg38,
		peak_assay = 'ATAC',
		rna_assay = "RNA",
		exclude_exons = TRUE)

library(BSgenome.Hsapiens.UCSC.hg38)
data(motifs)

seurat_obj <- find_motifs(
    			seurat_obj,
    			pfm = motifs,
    			genome = BSgenome.Hsapiens.UCSC.hg38,
			verbose = TRUE)

saveRDS(seurat_obj, file = paste0(out_dir, '/seurat_obj.rds'))


# remove duplicates: https://github.com/quadbiolab/Pando/issues/5#issuecomment-1016210826

regions <- NetworkRegions(seurat_obj)
cand_ranges <- regions@ranges
peak_ranges <- StringToGRanges(rownames(GetAssay(seurat_obj, assay='ATAC')))
peak_overlaps <- findOverlaps(cand_ranges, peak_ranges)

peak_overlaps_df <- data.frame(peak_overlaps)

peak_overlaps_df <- peak_overlaps_df[!duplicated(peak_overlaps_df[['queryHits']]), ]
peak_overlaps_df <- peak_overlaps_df[!duplicated(peak_overlaps_df[['subjectHits']]), ]

new_motif_data <- regions@motifs@data[peak_overlaps_df[['queryHits']], ]
new_regions_peaks <- peak_overlaps_df[['subjectHits']]
new_ranges_data <- regions@ranges[peak_overlaps_df[['queryHits']]]

seurat_obj@grn@regions@motifs@data <- new_motif_data
seurat_obj@grn@regions@peaks <- new_regions_peaks
seurat_obj@grn@regions@ranges <- new_ranges_data

saveRDS(seurat_obj, file = paste0(out_dir, '/seurat_obj.rds'))

seurat_obj <- readRDS(paste0(out_dir, '/seurat_obj.rds'))

library(foreach)
library(doParallel)
registerDoParallel(20)
seurat_obj <- infer_grn(
		    seurat_obj,
		    peak_to_gene_method = 'Signac',
		    method = 'glm',
		    upstream = 150000,
		    downstream = 150000,
                    parallel = FALSE)

```

V2 using same genes and regions as in SCENIC+ analysis

```bash

conda create -n Pando_env

conda activate Pando_env

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install -c bioconda r-signac

conda install -c bioconda bioconductor-bsgenome.hsapiens.ucsc.hg38

conda install -c bioconda bioconductor-ensdb.hsapiens.v86

conda install -c conda-forge r-devtools

conda install -c bioconda bioconductor-biovizbase

conda install -c r r-feather

conda install -c conda-forge r-dosnow

conda install -c r r-doparallel

conda install -c conda-forge r-rmpi

```

```r

devtools::install_github('quadbiolab/Pando')

```

```bash

conda activate Pando_env

```

# Preprocessing using Signac and Seurat

```r

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Matrix)
library(feather)
set.seed(1234)

# load RNA and ATAC data
GEX_counts_fname <- "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/data/count_matrices/scRNA_count_matrix.tsv"
GEX_counts <- read.table(
		file = GEX_counts_fname,
		header = TRUE,
		sep = '\t',
		quote = "",
		row.names = 1)
GEX_counts <- Matrix(data.matrix(GEX_counts), sparse = TRUE)

ACC_counts_fname <- "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/data/count_matrices/sc_fragment_matrix.feather"
ACC_counts <- read_feather(ACC_counts_fname)
ACC_counts <- as.matrix(ACC_counts)
rownames(ACC_counts) <- ACC_counts[, "index"]
ACC_counts <- ACC_counts[, -1]
ACC_counts <- Matrix(ACC_counts, sparse = TRUE)
colnames(ACC_counts) <- unlist(lapply(colnames(ACC_counts), function(x) gsub('___DPLC', '', x)))

fragpath <- "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/data/atac_fragments/fragments.sorted.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
seurat_obj <- CreateSeuratObject(
		counts = GEX_counts,
		assay = "RNA")

# create ATAC assay and add it to the object
seurat_obj[["ATAC"]] <- CreateChromatinAssay(
	  counts = ACC_counts,
	  sep = c(":", "-"),
	  fragments = fragpath,
	  annotation = annotation
	)

# Subset out special chromosomes, otherwise an error occurs due to: StringToGRanges(rownames(GetAssay(object, assay=peak_assay)))
chromosome_length <- read.table('/staging/leuven/stg_00002/lcb/resources/human/hg38/hg38.chrom.sizes')
chromosome_length <- chromosome_length[unlist(lapply(chromosome_length$V1, function(x) !grepl('_', x))), ]

regions_df <- data.frame(seurat_obj@assays$ATAC@ranges)

region_names_to_keep <- do.call(paste, 
	c(regions_df[regions_df$seqnames %in% unlist(list(chromosome_length$V1)), c('seqnames', 'start', 'end')], sep = "-"))

gene_names <- row.names(seurat_obj@assays$RNA@meta.features)

# Use genes and regions also used in SCENIC+ analysis.
scplus_regions = scan('/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/scplus_regions.txt', character())
scplus_regions = unlist(lapply(scplus_regions, function(x) gsub('-', '-', gsub(':', '-', x))))
scplus_genes = scan('/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/scplus_genes.txt', character())

regions_names_to_keep <- region_names_to_keep[region_names_to_keep %in% scplus_regions]
gene_names_to_keep <- gene_names[gene_names %in% scplus_genes]

seurat_obj <- subset(x = seurat_obj, features = c(region_names_to_keep, gene_names_to_keep))

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
#seurat_obj <- FindVariableFeatures(seurat_obj)

DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = NA)
seurat_obj <- RunTFIDF(seurat_obj)
seurat_obj <- RunSVD(seurat_obj, features = regions_names_to_keep)

out_dir <- "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/Pando/outs"

saveRDS(seurat_obj, file = paste0(out_dir, '/seurat_obj_V2.rds'))

seurat_obj <- readRDS(paste0(out_dir, '/seurat_obj_V2.rds'))


library(Pando)

data(SCREEN.ccRE.UCSC.hg38)

seurat_obj <- initiate_grn(
		object = seurat_obj,
		regions = SCREEN.ccRE.UCSC.hg38,
		peak_assay = 'ATAC',
		rna_assay = "RNA",
		exclude_exons = TRUE)

library(BSgenome.Hsapiens.UCSC.hg38)
data(motifs)

seurat_obj <- find_motifs(
    			seurat_obj,
    			pfm = motifs,
    			genome = BSgenome.Hsapiens.UCSC.hg38,
			verbose = TRUE)

saveRDS(seurat_obj, file = paste0(out_dir, '/seurat_obj_V2.rds'))


# remove duplicates: https://github.com/quadbiolab/Pando/issues/5#issuecomment-1016210826

regions <- NetworkRegions(seurat_obj)
cand_ranges <- regions@ranges
peak_ranges <- StringToGRanges(rownames(GetAssay(seurat_obj, assay='ATAC')))
peak_overlaps <- findOverlaps(cand_ranges, peak_ranges)

peak_overlaps_df <- data.frame(peak_overlaps)

peak_overlaps_df <- peak_overlaps_df[!duplicated(peak_overlaps_df[['queryHits']]), ]
peak_overlaps_df <- peak_overlaps_df[!duplicated(peak_overlaps_df[['subjectHits']]), ]

new_motif_data <- regions@motifs@data[peak_overlaps_df[['queryHits']], ]
new_regions_peaks <- peak_overlaps_df[['subjectHits']]
new_ranges_data <- regions@ranges[peak_overlaps_df[['queryHits']]]

seurat_obj@grn@regions@motifs@data <- new_motif_data
seurat_obj@grn@regions@peaks <- new_regions_peaks
seurat_obj@grn@regions@ranges <- new_ranges_data

#saveRDS(seurat_obj, file = paste0(out_dir, '/seurat_obj.rds'))

#seurat_obj <- readRDS(paste0(out_dir, '/seurat_obj.rds'))

library(foreach)
library(doParallel)
registerDoParallel(20)
seurat_obj <- infer_grn(
		    seurat_obj,
		    genes = gene_names_to_keep,
		    peak_to_gene_method = 'Signac',
		    method = 'glm',
		    upstream = 150000,
		    downstream = 150000,
                    parallel = TRUE)
saveRDS(seurat_obj, file = paste0(out_dir, '/seurat_obj_V2.rds'))

seurat_obj <- find_modules(seurat_obj)

modules <- NetworkModules(seurat_obj)

saveRDS(modules, file = paste0(out_dir, '/modules.rds'))

coefs <- coef(seurat_obj)

write.table(modules@meta, '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/Pando/outs/modules_pando.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

write.table(coefs, '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/Pando/outs/coefs_pando.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

```

