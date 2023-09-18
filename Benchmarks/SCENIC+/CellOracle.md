# 1. installation

```bash

conda create -n celloracle_env python=3.6
conda activate celloracle_env

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install numba cython pybedtools jupyter notebook
conda install -c anaconda cmake

pip install git+https://github.com/morris-lab/CellOracle.git

conda install -c conda-forge r-igraph


conda deactivate

conda create --name cicero
conda activate cicero
conda install -c bioconda r-monocle3
conda install -c conda-forge r-devtools

```

```bash

conda activate celloracle_env

```

# 2. Prepare input data

## 2.1 scRNA-seq data preparation: https://morris-lab.github.io/CellOracle.documentation/tutorials/scrnaprocess.html

```python

#import libraries
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.figdir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle/scanpy/plots'

#load data
f_RNA_count_matrix = '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/data/count_matrices/scRNA_count_matrix.tsv'
RNA_mtx = pd.read_csv(f_RNA_count_matrix, sep = '\t', index_col = 0).T

adata = sc.AnnData(RNA_mtx)

#filtering
sc.pp.filter_genes(adata, min_counts=1) #goes from 56,305 to 40,578

#normalization
sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')

#identification of highly variable genes
# Select top 3000 highly-variable genes
filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                              flavor='cell_ranger',
                                              n_top_genes=3000,
                                              log=False)

# Subset the genes
adata = adata[:, filter_result.gene_subset]

# Renormalize after filtering
sc.pp.normalize_per_cell(adata)

#Log transformation

# keep raw cont data before log transformation
adata.raw = adata
adata.layers["raw_count"] = adata.raw.X.copy()


# Log transformation and scaling
sc.pp.log1p(adata)
sc.pp.scale(adata)

#PCA and find neighbors

# PCA
sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True, save = 'PCA_var_ratio.pdf')

sc.pp.neighbors(adata, n_neighbors=4, n_pcs=7)


#Cell clustering
sc.tl.leiden(adata, resolution=0.8)

#Dimensional reduction
sc.tl.umap(adata)

adata.obs['line'] = [bc.split('_')[0] for bc in adata.obs.index]

sc.pl.pca(adata, color = 'line', save = '_line.pdf')
sc.pl.umap(adata, color = 'line', save = '_line.pdf')

#save
adata.write_h5ad('/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle/scanpy/outs/adata.h5ad')

```

## 2.2 Base GRN input data preparation: https://morris-lab.github.io/CellOracle.documentation/tutorials/base_grn.html

```bash

conda activate cicero

```

```r

install.pacakges('feather')

devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")

```

```r

library(cicero)
library(monocle3)
library(feather)

output_folder <- "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle/cicero"
f_scATAC_fragment_matrix <- "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/data/count_matrices/sc_fragment_matrix.feather"

#read and format fragments matrix
indata <- as.matrix(read_feather(f_scATAC_fragment_matrix))
rownames(indata) <- indata[, "index"]
indata <- indata[, -1]

#format cell and peak info
chromosomes <- unlist(lapply(rownames(indata), function(x) strsplit(x, ':')[[1]][1]))
start_end <-  unlist(lapply(rownames(indata), function(x) strsplit(x, ':')[[1]][2]))
start <- unlist(lapply(start_end, function(x) strsplit(x, '-')[[1]][1]))
end <- unlist(lapply(start_end, function(x) strsplit(x, '-')[[1]][2]))

region_names <- paste(chromosomes, start, end, sep = '_')

peakinfo <- data.frame(
		row.names = region_names,
		chr = chromosomes,
		bp1 = as.numeric(start),
		bp2 = as.numeric(end))

row.names(indata) <- row.names(peakinfo)

# Make CDS
input_cds <- new_cell_data_set(
		indata,
		gene_metadata = peakinfo)

input_cds <- monocle3::detect_genes(input_cds)

#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

pdf(paste0(output_folder, "/plots/peak_counts.pdf"))
hist(Matrix::colSums(exprs(input_cds)))
dev.off()

# Data preprocessing
set.seed(2017)

input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI", verbose = TRUE)

saveRDS(input_cds, paste0(output_folder, '/outs/input_cds.Rds'))
input_cds <- readRDS(paste0(output_folder, '/outs/input_cds.Rds'))

# Dimensional reduction with umap
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

saveRDS(cicero_cds, paste0(output_folder, '/outs/cicero_cds.Rds'))

#load reference genome information

chromosome_length <- read.table('/staging/leuven/stg_00002/lcb/resources/human/hg38/hg38.chrom.sizes')
#remove unknown chroms
chromosome_length <- chromosome_length[unlist(lapply(chromosome_length$V1, function(x) !grepl('_', x))), ]
row.names(chromosome_length) <- NULL

#run Cicero
conns <- run_cicero(cicero_cds, chromosome_length)

saveRDS(conns, paste0(output_folder, '/outs/cicero_connections.Rds'))

# Save results for the next step

all_peaks <- row.names(exprs(input_cds))

write.csv(x = all_peaks, file = paste0(output_folder, "/outs/all_peaks.csv"))
write.csv(x = conns, file = paste0(output_folder, "/outs/cicero_connections.csv"))

```

```bash

conda activate celloracle_env

```

```python

output_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle/cicero"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, shutil, importlib, glob
from tqdm import tqdm
from celloracle import motif_analysis as ma
import celloracle as co
co.__version__
# '0.8.5'

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

# Load scATAC-seq peak list.
peaks = pd.read_csv(os.path.join(output_folder, "outs/all_peaks.csv"), index_col=0)
peaks = peaks.x.values
peaks

# Load cicero coaccess score.
cicero_connections = pd.read_csv(os.path.join(output_folder, "outs/cicero_connections.csv"), index_col=0)
cicero_connections.head()

# Get TSS annotation
tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="hg38")
tss_annotated.tail()

# Integrate TSS info and cicero connections
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated, cicero_connections=cicero_connections)
print(integrated.shape)
integrated.head()

# Filter peaks
peak = integrated[integrated.coaccess >= 0.8]
peak = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)
print(peak.shape)
peak.head()

# Save data
peak.to_csv(os.path.join(output_folder, "outs/processed_peak_file.csv"))

```

## 2.3 Transcription factor binding motif scan: https://morris-lab.github.io/CellOracle.documentation/tutorials/motifscan.html

### 2.3.1 using CellOracle motif enrichment tool

```bash

conda activate celloracle_env

```

```python

data_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle/cicero"
output_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle/celloracle"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


import seaborn as sns

import os, sys, shutil, importlib, glob
from tqdm import tqdm

from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object

plt.rcParams['figure.figsize'] = (15,7)
plt.rcParams["savefig.dpi"] = 600

# Rerefence genome data preparation

ref_genome = "hg38"

genome_installation = ma.is_genome_installed(ref_genome=ref_genome)
print(ref_genome, "installation: ", genome_installation)

if not genome_installation:
    import genomepy
    genomepy.install_genome(ref_genome, "UCSC")
else:
    print(ref_genome, "is installed.")

# Load data
# Load annotated peak data.
peaks = pd.read_csv(os.path.join(data_folder, "outs/processed_peak_file.csv"), index_col=0)
peaks.head()

# Define function for quality check
def decompose_chrstr(peak_str):
    """
    Args:
        peak_str (str): peak_str. e.g. 'chr1_3094484_3095479'

    Returns:
        tuple: chromosome name, start position, end position
    """

    *chr_, start, end = peak_str.split("_")
    chr_ = "_".join(chr_)
    return chr_, start, end

from genomepy import Genome

def check_peak_foamat(peaks_df, ref_genome):
    """
    Check peak fomat.
     (1) Check chromosome name.
     (2) Check peak size (length) and remove sort DNAs (<5bp)

    """

    df = peaks_df.copy()

    n_peaks_before = df.shape[0]

    # Decompose peaks and make df
    decomposed = [decompose_chrstr(peak_str) for peak_str in df["peak_id"]]
    df_decomposed = pd.DataFrame(np.array(decomposed))
    df_decomposed.columns = ["chr", "start", "end"]
    df_decomposed["start"] = df_decomposed["start"].astype(np.int)
    df_decomposed["end"] = df_decomposed["end"].astype(np.int)

    # Load genome data
    genome_data = Genome(ref_genome)
    all_chr_list = list(genome_data.keys())


    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])


    # Filter peaks with invalid chromosome name
    n_threshold = 5
    df = df[(lengths >= n_threshold) & df_decomposed.chr.isin(all_chr_list)]

    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])

    # Data counting
    n_invalid_length = len(lengths[lengths < n_threshold])
    n_peaks_invalid_chr = n_peaks_before - df_decomposed.chr.isin(all_chr_list).sum()
    n_peaks_after = df.shape[0]

    #
    print("Peaks before filtering: ", n_peaks_before)
    print("Peaks with invalid chr_name: ", n_peaks_invalid_chr)
    print("Peaks with invalid length: ", n_invalid_length)
    print("Peaks after filtering: ", n_peaks_after)

    return df

peaks = check_peak_foamat(peaks, ref_genome)

"""
Peaks before filtering:  39178
Peaks with invalid chr_name:  0
Peaks with invalid length:  0
Peaks after filtering:  39178
"""

# Instantiate TFinfo object
tfi = ma.TFinfo(peak_data_frame=peaks,
                ref_genome=ref_genome)

%%time
# Scan motifs. !!CAUTION!! This step may take several hours if you have many peaks!
print('Starting')
tfi.scan(fpr=0.02,
         motifs=None,  # If you enter None, default motifs will be loaded.
         verbose=True)
# Save tfinfo object
tfi.to_hdf5(file_path = os.path.join(output_folder, 'outs/test1.celloracle.tfinfo'))

tfi.scanned_df.head()

# Filtering motifs
# Reset filtering
tfi.reset_filtering()

# Do filtering
tfi.filter_motifs_by_score(threshold=10)
# Filtering finished: 10559374 -> 2226333

# Do post filtering process. Convert results into several file format.
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

#get final results

df = tfi.to_dataframe()
df.head()

# Save result as a dataframe
df = tfi.to_dataframe()
df.to_parquet(os.path.join(output_folder, "outs/base_GRN_dataframe.parquet"))

# If you want, you can save the result as a dictionary as follows.
td = tfi.to_dictionary(dictionary_type="targetgene2TFs")
save_as_pickled_object(td, os.path.join(output_folder, "outs/TFinfo_targetgene2TFs.pickled"))

```

### 2.3.2 using cistarget motifs (may take a very long time?)

## 3. GRN model construction and Network analysis

```bash

module load R/4.0.2-foss-2018a-bare

```

```python

data_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle"

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

import celloracle as co
co.__version__
#'0.8.5 

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

co.test_R_libraries_installation()

#load data

adata = sc.read_h5ad(os.path.join(data_folder, 'scanpy/outs/adata.h5ad'))

base_GRN = pd.read_parquet(os.path.join(data_folder, 'celloracle/outs/base_GRN_dataframe.parquet'))

# Instantiate Oracle object
oracle = co.Oracle()

adata.X = adata.layers["raw_count"].copy()

oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="line",
                                   embedding_name="X_umap")

oracle.import_TF_data(TF_info_matrix=base_GRN)

# Perform PCA
oracle.perform_PCA()

# Select important PCs
fig, ax = plt.subplots()
eax.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
fig.savefig(os.path.join(data_folder, 'celloracle/plots/pca_co.pdf'))
print(n_comps)
n_comps = min(n_comps, 50)

#estimate optimal number of nearest neighbors for KNN imputation
n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")

k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")

oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)

oracle.to_hdf5(os.path.join(data_folder, 'celloracle/outs/oracle_obj.celloracle.oracle'))

%%time
# Calculate GRN for each population in "louvain_annot" clustering unit.
# This step may take long time.(~30 minutes)
links = oracle.get_links(cluster_name_for_GRN_unit="line", alpha=10,
                         verbose_level=10, test_mode=False)

# network preprocessing
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

links.plot_degree_distributions(plot_model=True,
                                save=os.path.join(data_folder, 'celloracle/plots/degree_distribution/'))

links.get_score()

links.to_hdf5(file_path=os.path.join(data_folder, "celloracle/outs/links.celloracle.links"))

for cluster in links.cluster:
	links.plot_scores_as_rank(cluster=cluster, n_gene=30, save=os.path.join(data_folder, f'celloracle/plots/ranked_score_{cluster}'))

#format links for benchmarking
df_unfiltered_links = []
for cluster in links.links_dict.keys():
	tmp = links.links_dict[cluster].copy()
	tmp['cluster'] = cluster
	df_unfiltered_links.append(tmp)

df_unfiltered_links = pd.concat(df_unfiltered_links, axis = 0)

df_filtered_links = []
for cluster in links.filtered_links.keys():
	tmp = links.filtered_links[cluster].copy()
	tmp['cluster'] = cluster
	df_filtered_links.append(tmp)

df_filtered_links = pd.concat(df_filtered_links, axis = 0)

df_unfiltered_links.to_csv(os.path.join(data_folder, 'celloracle/outs/target_genes_unfiltered.tsv'), sep = '\t', header = True, index = False)
df_filtered_links.to_csv(os.path.join(data_folder, 'celloracle/outs/target_genes_filtered.tsv'), sep = '\t', header = True, index = False)


```

Rerun using same genes and regions as in SCENIC+ analysis + proper TSS annotation

# 1. installation

```bash

conda create -n celloracle_env python=3.6
conda activate celloracle_env

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install numba cython pybedtools jupyter notebook
conda install -c anaconda cmake

pip install git+https://github.com/morris-lab/CellOracle.git

conda install -c conda-forge r-igraph


conda deactivate

conda create --name cicero
conda activate cicero
conda install -c bioconda r-monocle3
conda install -c conda-forge r-devtools

```

```bash

conda activate celloracle_env

```

# 2. Prepare input data

## 2.1 scRNA-seq data preparation: https://morris-lab.github.io/CellOracle.documentation/tutorials/scrnaprocess.html

```python

#import libraries
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.figdir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle/scanpy/plots'

#load data
f_RNA_count_matrix = '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/data/count_matrices/scRNA_count_matrix.tsv'
RNA_mtx = pd.read_csv(f_RNA_count_matrix, sep = '\t', index_col = 0).T

adata = sc.AnnData(RNA_mtx)

#filtering
#sc.pp.filter_genes(adata, min_counts=1) #goes from 56,305 to 40,578

#normalization
sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')

#identification of highly variable genes
# Select top 3000 highly-variable genes
#filter_result = sc.pp.filter_genes_dispersion(adata.X,
#                                              flavor='cell_ranger',
#                                              n_top_genes=3000,
#                                              log=False)

#load genes used for SCENIC+ analysis
with open('/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/scplus_genes.txt', 'r') as f:
	scplus_genes = f.read().splitlines()


# Subset the genes
adata = adata[:, scplus_genes]

# Renormalize after filtering
sc.pp.normalize_per_cell(adata)

#Log transformation

# keep raw cont data before log transformation
adata.raw = adata
adata.layers["raw_count"] = adata.raw.X.copy()


# Log transformation and scaling
sc.pp.log1p(adata)
sc.pp.scale(adata)

#PCA and find neighbors

# PCA
sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True, save = '_PCA_var_ratio_V2.pdf')

sc.pp.neighbors(adata, n_neighbors=4, n_pcs=7)


#Cell clustering
sc.tl.leiden(adata, resolution=0.8)

#Dimensional reduction
sc.tl.umap(adata)

adata.obs['line'] = [bc.split('_')[0] for bc in adata.obs.index]

sc.pl.pca(adata, color = 'line', save = '_line_V2.pdf')
sc.pl.umap(adata, color = 'line', save = '_line_V2.pdf')

#save
adata.write_h5ad('/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle/scanpy/outs/adata_V2.h5ad')

```

## 2.2 Base GRN input data preparation: https://morris-lab.github.io/CellOracle.documentation/tutorials/base_grn.html

```bash

conda activate cicero

```

```r

install.pacakges('feather')

devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")

```

```r

library(cicero)
library(monocle3)
library(feather)

output_folder <- "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle/cicero"
f_scATAC_fragment_matrix <- "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/data/count_matrices/sc_fragment_matrix.feather"

#read and format fragments matrix
indata <- as.matrix(read_feather(f_scATAC_fragment_matrix))
rownames(indata) <- indata[, "index"]
indata <- indata[, -1]

#format cell and peak info
chromosomes <- unlist(lapply(rownames(indata), function(x) strsplit(x, ':')[[1]][1]))
start_end <-  unlist(lapply(rownames(indata), function(x) strsplit(x, ':')[[1]][2]))
start <- unlist(lapply(start_end, function(x) strsplit(x, '-')[[1]][1]))
end <- unlist(lapply(start_end, function(x) strsplit(x, '-')[[1]][2]))

region_names <- paste(chromosomes, start, end, sep = '_')

peakinfo <- data.frame(
		row.names = region_names,
		chr = chromosomes,
		bp1 = as.numeric(start),
		bp2 = as.numeric(end))

row.names(indata) <- row.names(peakinfo)

# Make CDS
input_cds <- new_cell_data_set(
		indata,
		gene_metadata = peakinfo)

input_cds <- monocle3::detect_genes(input_cds)


# subset for regions used in SCENIC+ analysis
scplus_regions = scan('/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/scplus_regions.txt', character())
scplus_regions = unlist(lapply(scplus_regions, function(x) gsub('-', '_', gsub(':', '_', x))))

input_cds <- input_cds[scplus_regions]

#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

# Data preprocessing
set.seed(2017)

input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI", verbose = TRUE)

saveRDS(input_cds, paste0(output_folder, '/outs/input_cds_V2.Rds'))
input_cds <- readRDS(paste0(output_folder, '/outs/input_cds_V2.Rds'))

# Dimensional reduction with umap
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

saveRDS(cicero_cds, paste0(output_folder, '/outs/cicero_cds_V2.Rds'))

#load reference genome information

chromosome_length <- read.table('/staging/leuven/stg_00002/lcb/resources/human/hg38/hg38.chrom.sizes')
#remove unknown chroms
chromosome_length <- chromosome_length[unlist(lapply(chromosome_length$V1, function(x) !grepl('_', x))), ]
row.names(chromosome_length) <- NULL

#run Cicero
#conns <- run_cicero(cicero_cds, chromosome_length, window = 150000)

distance_parameters <- estimate_distance_parameter(cicero_cds, window = 150000, distance_constraint = 75000, genomic_coords = chromosome_length)
distance_parameter <- mean(unlist(distance_parameters))
cicero_cds <- generate_cicero_models(cicero_cds, distance_parameter,  window = 150000, genomic_coords = chromosome_length)
conns <- assemble_connections(cicero_cds)

saveRDS(conns, paste0(output_folder, '/outs/cicero_connections_V2.Rds'))

# Save results for the next step

all_peaks <- row.names(exprs(input_cds))

write.csv(x = all_peaks, file = paste0(output_folder, "/outs/all_peaks_V2.csv"))
write.csv(x = conns, file = paste0(output_folder, "/outs/cicero_connections_V2.csv"))

```

```bash

conda activate celloracle_env

```

```python

output_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle/cicero"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, shutil, importlib, glob
from tqdm import tqdm
from celloracle import motif_analysis as ma
import celloracle as co
co.__version__
# '0.8.5'

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

# Load scATAC-seq peak list.
peaks = pd.read_csv(os.path.join(output_folder, "outs/all_peaks_V2.csv"), index_col=0)
peaks = peaks.x.values
peaks

# Load cicero coaccess score.
cicero_connections = pd.read_csv(os.path.join(output_folder, "outs/cicero_connections_V2.csv"), index_col=0)
cicero_connections.head()

# Get TSS annotation (manually, the annotation included in celloracle is not compatible with our data)
#tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="hg38")
#tss_annotated.tail()

import pybiomart as pbm
dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://oct2016.archive.ensembl.org/')
annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
annot['Chromosome Name'] = annot['Chromosome Name'].astype('str')
filter = annot['Chromosome Name'].str.contains('CHR|GL|JH|MT')
annot = annot[~filter]
annot['Chromosome Name'] = annot['Chromosome Name'].str.replace(r'(\b\S)', r'chr\1')
annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
annot = annot[annot.Transcript_type == 'protein_coding']
annot['End'] = annot['Start'] + 1
annot['Strand'] = ['+' if s == 1 else '-' for s in annot['Strand']]
annot.drop('Transcript_type', axis = 1)
annot['score'] = 0
annot = annot[['Chromosome', 'Start', 'End', 'Gene', 'score', 'Strand']]
annot.rename({'Chromosome': 'chr', 'Start': 'start', 'End': 'end', 'Gene': 'gene_short_name', 'Strand': 'strand'}, axis = 1, inplace = True)
annot.reset_index(drop = True, inplace = True)
annot.sort_values(['chr', 'start']).to_csv(
	'/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/hg38_tss_annot_ensembl_86.bed', 
	sep = '\t',
	header = False,
	index = False)

tss_annotated = ma.get_tss_info(
	peak_str_list = peaks, 
	ref_genome=None, 
	custom_tss_file_path = '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/hg38_tss_annot_ensembl_86.bed')


# Integrate TSS info and cicero connections
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated, cicero_connections=cicero_connections)
print(integrated.shape)
integrated.head()

integrated_formatted = integrated.copy()
integrated_formatted['peak_id'] = [f"{x.split('_')[0]}:{x.split('_')[1]}-{x.split('_')[2]}" for x in integrated['peak_id']]

integrated_formatted.to_csv(os.path.join(output_folder, 'outs/peak_to_gene_V2.tsv'), sep = '\t', header = True, index = False)

# Filter peaks
peak = integrated[integrated.coaccess >= 0.8]
peak = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)
print(peak.shape)
peak.head()

# Save data
peak.to_csv(os.path.join(output_folder, "outs/processed_peak_file_V2.csv"))

```

## 2.3 Transcription factor binding motif scan: https://morris-lab.github.io/CellOracle.documentation/tutorials/motifscan.html

### 2.3.1 using CellOracle motif enrichment tool

```bash

conda activate celloracle_env

```

```python

data_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle/cicero"
output_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle/celloracle"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns

import os, sys, shutil, importlib, glob
from tqdm import tqdm

from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object

plt.rcParams['figure.figsize'] = (15,7)
plt.rcParams["savefig.dpi"] = 600

# Rerefence genome data preparation

ref_genome = "hg38"

genome_installation = ma.is_genome_installed(ref_genome=ref_genome)
print(ref_genome, "installation: ", genome_installation)

if not genome_installation:
    import genomepy
    genomepy.install_genome(ref_genome, "UCSC")
else:
    print(ref_genome, "is installed.")

# Load data
# Load annotated peak data.
peaks = pd.read_csv(os.path.join(data_folder, "outs/processed_peak_file_V2.csv"), index_col=0)
peaks.head()

# Define function for quality check
def decompose_chrstr(peak_str):
    """
    Args:
        peak_str (str): peak_str. e.g. 'chr1_3094484_3095479'

    Returns:
        tuple: chromosome name, start position, end position
    """

    *chr_, start, end = peak_str.split("_")
    chr_ = "_".join(chr_)
    return chr_, start, end

from genomepy import Genome

def check_peak_foamat(peaks_df, ref_genome):
    """
    Check peak fomat.
     (1) Check chromosome name.
     (2) Check peak size (length) and remove sort DNAs (<5bp)

    """

    df = peaks_df.copy()

    n_peaks_before = df.shape[0]

    # Decompose peaks and make df
    decomposed = [decompose_chrstr(peak_str) for peak_str in df["peak_id"]]
    df_decomposed = pd.DataFrame(np.array(decomposed))
    df_decomposed.columns = ["chr", "start", "end"]
    df_decomposed["start"] = df_decomposed["start"].astype(np.int)
    df_decomposed["end"] = df_decomposed["end"].astype(np.int)

    # Load genome data
    genome_data = Genome(ref_genome)
    all_chr_list = list(genome_data.keys())


    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])


    # Filter peaks with invalid chromosome name
    n_threshold = 5
    df = df[(lengths >= n_threshold) & df_decomposed.chr.isin(all_chr_list)]

    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])

    # Data counting
    n_invalid_length = len(lengths[lengths < n_threshold])
    n_peaks_invalid_chr = n_peaks_before - df_decomposed.chr.isin(all_chr_list).sum()
    n_peaks_after = df.shape[0]

    #
    print("Peaks before filtering: ", n_peaks_before)
    print("Peaks with invalid chr_name: ", n_peaks_invalid_chr)
    print("Peaks with invalid length: ", n_invalid_length)
    print("Peaks after filtering: ", n_peaks_after)

    return df

peaks = check_peak_foamat(peaks, ref_genome)

"""
			   V1     V2
Peaks before filtering:  39178	14947
Peaks with invalid chr_name:  0	0
Peaks with invalid length:  0	0
Peaks after filtering:  39178	14947
"""

# Instantiate TFinfo object
tfi = ma.TFinfo(peak_data_frame=peaks,
                ref_genome=ref_genome)

#%%time
# Scan motifs. !!CAUTION!! This step may take several hours if you have many peaks!
print('Starting')
tfi.scan(fpr=0.02,
         motifs=None,  # If you enter None, default motifs will be loaded.
         verbose=True)
# Save tfinfo object
tfi.to_hdf5(file_path = os.path.join(output_folder, 'outs/test1_V2.celloracle.tfinfo'))

tfi.scanned_df.head()

# Filtering motifs
# Reset filtering
tfi.reset_filtering()

# Do filtering
tfi.filter_motifs_by_score(threshold=10)
# V1: Filtering finished: 10559374 -> 2226333
# V2: Filtering finished: 3294962 -> 711923

# Do post filtering process. Convert results into several file format.
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

#get final results

df = tfi.to_dataframe()
df.head()

# Save result as a dataframe
df = tfi.to_dataframe()
df.to_parquet(os.path.join(output_folder, "outs/base_GRN_dataframe_V2.parquet"))

# If you want, you can save the result as a dictionary as follows.
td = tfi.to_dictionary(dictionary_type="targetgene2TFs")
save_as_pickled_object(td, os.path.join(output_folder, "outs/TFinfo_targetgene2TFs_V2.pickled"))

#generate cistrome like dictionary
df['formated_region_name'] = [f"{x.split('_')[0]}:{x.split('_')[1]}-{x.split('_')[2]}" for x in df['peak_id']]

cistromes = {x: df.loc[df[x] == 1, 'formated_region_name'].tolist() for x in set(df.columns) - set(['peak_id', 'gene_short_name', 'formated_region_name'])}

import pickle
with open(os.path.join(output_folder, 'cistromes.pkl'), 'wb') as f:
	pickle.dump(cistromes, f)

```

### 2.3.2 using cistarget motifs (may take a very long time?)

## 3. GRN model construction and Network analysis

```bash

module load R/4.0.2-foss-2018a-bare

```

```python

data_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/encode_validation/CellOracle"

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

import celloracle as co
co.__version__
#'0.8.5 

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

co.test_R_libraries_installation()

#load data

adata = sc.read_h5ad(os.path.join(data_folder, 'scanpy/outs/adata_V2.h5ad'))

base_GRN = pd.read_parquet(os.path.join(data_folder, 'celloracle/outs/base_GRN_dataframe_V2.parquet'))

# Instantiate Oracle object
oracle = co.Oracle()

adata.X = adata.layers["raw_count"].copy()

oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="line",
                                   embedding_name="X_umap")

oracle.import_TF_data(TF_info_matrix=base_GRN)

# Perform PCA
oracle.perform_PCA()

# Select important PCs
fig, ax = plt.subplots()
ax.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
fig.savefig(os.path.join(data_folder, 'celloracle/plots/pca_co_V2.pdf'))
print(n_comps)
n_comps = min(n_comps, 50)

#estimate optimal number of nearest neighbors for KNN imputation
n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")

k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")

oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)

oracle.to_hdf5(os.path.join(data_folder, 'celloracle/outs/oracle_obj_V2.celloracle.oracle'))

%%time
# Calculate GRN for each population in "louvain_annot" clustering unit.
# This step may take long time.(~30 minutes)
links = oracle.get_links(cluster_name_for_GRN_unit="line", alpha=10,
                         verbose_level=10, test_mode=False)

# network preprocessing
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

links.plot_degree_distributions(plot_model=True,
                                save=os.path.join(data_folder, 'celloracle/plots/degree_distribution_V2/'))

links.get_score()

links.to_hdf5(file_path=os.path.join(data_folder, "celloracle/outs/links_V2.celloracle.links"))

for cluster in links.cluster:
	links.plot_scores_as_rank(cluster=cluster, n_gene=30, save=os.path.join(data_folder, f'celloracle/plots/ranked_score_V2_{cluster}'))

#format links for benchmarking
df_unfiltered_links = []
for cluster in links.links_dict.keys():
	tmp = links.links_dict[cluster].copy()
	tmp['cluster'] = cluster
	df_unfiltered_links.append(tmp)

df_unfiltered_links = pd.concat(df_unfiltered_links, axis = 0)

df_filtered_links = []
for cluster in links.filtered_links.keys():
	tmp = links.filtered_links[cluster].copy()
	tmp['cluster'] = cluster
	df_filtered_links.append(tmp)

df_filtered_links = pd.concat(df_filtered_links, axis = 0)

df_unfiltered_links.to_csv(os.path.join(data_folder, 'celloracle/outs/target_genes_unfiltered_V2.tsv'), sep = '\t', header = True, index = False)
df_filtered_links.to_csv(os.path.join(data_folder, 'celloracle/outs/target_genes_filtered_V2.tsv'), sep = '\t', header = True, index = False)

#filter links with default params
links = co.load_hdf5(file_path=os.path.join(data_folder, "celloracle/outs/links_V2.celloracle.links"))
links.filter_links(p=0.001, weight="coef_abs", threshold_number=10000)
df_filtered_links = []
for cluster in links.filtered_links.keys():
        tmp = links.filtered_links[cluster].copy()
        tmp['cluster'] = cluster
        df_filtered_links.append(tmp)
df_filtered_links = pd.concat(df_filtered_links, axis = 0)
df_filtered_links.to_csv(os.path.join(data_folder, 'celloracle/outs/target_genes_filtered_V2_default_params.tsv'), sep = '\t', header = True, index = False)


```
