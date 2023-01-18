# Code for generating the pannels of figure 2 and accompanying supplementary figures

## Set up python environment & load data

```python

## Set up environment

import os
import dill
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

work_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_granulocyte_sorted_10k/scenicplus_clustered_consensus'
plot_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_pbmc/v3_clustered_consensus/plots'
ranking_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_granulocyte_sorted_10k'


color_dict = {
        'B cells': '#FF5714',
        'CD14+ Monocytes': '#6EEB83',
        'FCGR3A+ Monocytes': '#65B891',
        'Dendritic cells 1': '#575761',
        'Dendritic cells 2': '#1BE7FF',
        'CD4 T cells': '#1D3461',
        'CD8 T cells': '#3D5A6C',
        'NK cells': '#64B6AC'
}

## Load data

scplus_obj = dill.load(
	open(os.path.join(work_dir, 'scplus_obj.pkl'), 'rb')
)

region_ranking = dill.load(
        open(os.path.join(ranking_dir, 'scenicplus/region_ranking.pkl'), 'rb')
)

gene_ranking = dill.load(
        open(os.path.join(ranking_dir, 'scenicplus/gene_ranking.pkl'), 'rb')
)

```


## Supplementary figure S9, pannel a: Filtering eRegulons

```python

## Filter eRegulons

from scenic.preprocessing.filtering import apply_std_filtering_to_eRegulons

apply_std_filtering_to_eRegulons(scplus_obj)

score_eRegulons(scplus_obj,
                ranking = region_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                key_added = 'eRegulon_AUC_filtered',
                enrichment_type= 'region',
                auc_threshold = 0.05,
                normalize = False,
                n_cpu = 5)
score_eRegulons(scplus_obj,
                gene_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                key_added = 'eRegulon_AUC_filtered',
                enrichment_type = 'gene',
                auc_threshold = 0.05,
                normalize= False,
                n_cpu = 5)

from scenicplus.cistromes import TF_cistrome_correlation, generate_pseudobulks

generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'GEX_celltype',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Gene_based'
)
generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'GEX_celltype',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Region_based'
)

TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'GEX_celltype',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Gene_based',
            out_key = 'filtered_gene_based')
TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'GEX_celltype',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Region_based',
            out_key = 'filtered_region_based')

n_targets = [int(x.split('(')[1].replace('r)', '')) for x in scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Cistrome']]
rho = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'].to_list()
adj_pval = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Adjusted_p-value'].to_list()

thresholds = {
        'rho': [-0.8, 0.7],
        'n_targets': 0
}
import seaborn as sns
plot_dir = plotdir
fig, ax = plt.subplots()
sc = ax.scatter(rho, n_targets, c = -np.log10(adj_pval), s = 5)
ax.set_xlabel('Correlation coefficient')
ax.set_ylabel('nr. target regions')
#ax.hlines(y = thresholds['n_targets'], xmin = min(rho), xmax = max(rho), color = 'black', ls = 'dashed', lw = 1)
ax.vlines(x = thresholds['rho'], ymin = 0, ymax = max(n_targets), color = 'black', ls = 'dashed', lw = 1)
ax.text(x = thresholds['rho'][0], y = max(n_targets), s = str(thresholds['rho'][0]))
ax.text(x = thresholds['rho'][1], y = max(n_targets), s = str(thresholds['rho'][1]))
sns.despine(ax = ax)
fig.colorbar(sc, label = '-log10(adjusted_pvalue)', ax = ax)
fig.savefig(os.path.join(plot_dir, 'eRegulon_selection_n_targets_rho.pdf'))

selected_cistromes = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based'].loc[
        np.logical_or(
                scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] > thresholds['rho'][1],
                scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] < thresholds['rho'][0]
        )]['Cistrome'].to_list()
selected_eRegulons = [x.split('_(')[0] for x in selected_cistromes]
selected_eRegulons_gene_sig = [
        x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys()
        if x.split('_(')[0] in selected_eRegulons]
selected_eRegulons_region_sig = [
        x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'].keys()
        if x.split('_(')[0] in selected_eRegulons]

scplus_obj.uns['selected_eRegulon'] = {'Gene_based': selected_eRegulons_gene_sig, 'Region_based': selected_eRegulons_region_sig}

## Set celltype names
scplus_obj.metadata_cell['celltype_name'] = [x.replace('_', ' ') for x in scplus_obj.metadata_cell['GEX_celltype']]

## Save changes
dill.dump(scplus_obj, open(os.path.join(work_dir, 'scplus_obj.pkl'), 'wb'))

```

## Export to loom for visualizations in R

```python
binarize_AUC(scplus_obj, 
             auc_key='eRegulon_AUC_filtered',
             out_key='eRegulon_AUC_filtered_thresholds',
             signature_keys=['Gene_based'],
             n_cpu=5)

def remove_second_sign(x):
	if 'extended' not in x:
		TF, first, second, n = x.split('_')
		return f'{TF}_{first}_{n}'
	else:
		TF, extended, first, second, n = x.split('_')
		return f'{TF}_{extended}_{first}_{n}'

scplus_obj.uns['eRegulon_metadata_filtered']['Gene_signature_name'] = [remove_second_sign(x) for x in scplus_obj.uns['eRegulon_metadata_filtered']['Gene_signature_name']]
scplus_obj.uns['eRegulon_metadata_filtered']['Region_signature_name'] = [remove_second_sign(x) for x in scplus_obj.uns['eRegulon_metadata_filtered']['Region_signature_name']]

del(scplus_obj.dr_cell['GEX_X_pca'])
del(scplus_obj.dr_cell['GEX_rep'])

from scenicplus.loom import *
export_to_loom(scplus_obj, 
       signature_key = 'Gene_based',
       eRegulon_metadata_key = 'eRegulon_metadata_filtered',
       auc_key = 'eRegulon_AUC_filtered',
       auc_thr_key = 'eRegulon_AUC_filtered_thresholds',
       keep_direct_and_extended_if_not_direct = True,
       tree_structure = (),
       title =  'Gene based eGRN',
       nomenclature = 'hg38',
       out_fname=os.path.join(work_dir,'SCENIC+_gene_based.loom'))

binarize_AUC(scplus_obj,
             auc_key='eRegulon_AUC_filtered',
             out_key='eRegulon_AUC_filtered_thresholds',
             signature_keys=['Region_based'],
             n_cpu=5)

export_to_loom(scplus_obj, 
       signature_key = 'Region_based',
       eRegulon_metadata_key = 'eRegulon_metadata_filtered',
       auc_key = 'eRegulon_AUC_filtered',
       auc_thr_key = 'eRegulon_AUC_filtered_thresholds',
       keep_direct_and_extended_if_not_direct = True,
       tree_structure = (),
       title =  'Region based eGRN',
       nomenclature = 'hg38',
       out_fname=os.path.join(work_dir,'SCENIC+_region_based.loom'))

with open(os.path.join(work_dir, 'selected_eRegulons.txt'), 'w') as f:
	for eRegulon_name in scplus_obj.uns['selected_eRegulon']['Gene_based']:
		f.write(eRegulon_name.rsplit('_', 1)[0])
		f.write('\n')

```

## Main Figure 2, pannel A: eRegulons tSNE

```r

work_dir <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_granulocyte_sorted_10k/scenicplus_clustered_consensus'
plot_dir <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_pbmc/v3_clustered_consensus/plots'

# load data
library(SCopeLoomR)
loom <- open_loom(paste0(work_dir, '/SCENIC+_gene_based.loom'))
cell_data <- get_cell_annotation(loom)
embeddings <- get_embeddings(loom)
eRegulons_tSNE <- embeddings$eRegulons_tSNE
colnames(eRegulons_tSNE) <- c('tSNE_1', 'tSNE_2')
cell_data[] <- lapply(cell_data, sub, pattern = " ", replacement = "-")
cell_data[] <- lapply(cell_data, sub, pattern = "-(.*)", replacement = "")

# Color by cell type
cell_plot_data <- cbind(eRegulons_tSNE, cell_data[rownames(eRegulons_tSNE),])

color_dict = c(
        'B_cells' = '#FF5714',
        'CD14+_Monocytes' = '#6EEB83',
        'FCGR3A+_Monocytes' = '#65B891',
        'Dendritic_cells_1' = '#575761',
        'Dendritic_cells_2' = '#1BE7FF',
        'CD4_T_cells' = '#1D3461',
        'CD8_T_cells' = '#3D5A6C',
        'NK_cells' = '#64B6AC')

#plot
library(ggplot2)
source('/staging/leuven/stg_00002/lcb/cbravo/software/plotting_aux.R')
plot <- ggplot(cell_plot_data, aes(x=tSNE_1, y=tSNE_2, colour=ACC_celltype))
plot <- plot + geom_point(size = 0.2)
plot <- plot + theme_classic()
plot <- plot + theme(legend.position = "none")
plot <- plot + labs(x = NULL, y = NULL)
plot <- plot + scale_color_manual(values = color_dict)
plot <- plot + guides(x = "none", y = "none")  
ggsave(filename = paste0(plot_dir, '/R_eRegulons_tSNE.png'), device='png', bg = "transparent", width=5.24, height=4.5)

pdf(paste0(plot_dir, '/R_eRegulons_tSNE.pdf'), width = 5.24, height = 4.5)
LabelClusters(plot, 'ACC_celltype', split.by ='ACC_celltype', box=FALSE, repel=TRUE)  + scale_color_manual(values = color_dict)
dev.off()

```

## Main Figure 2, Pannel b: n regions per gene

```python

import os
import dill
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

work_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_granulocyte_sorted_10k/scenicplus_clustered_consensus'
plot_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_pbmc/v3_clustered_consensus/plots'

color_dict = {
        'B cells': '#FF5714',
        'CD14+ Monocytes': '#6EEB83',
        'FCGR3A+ Monocytes': '#65B891',
        'Dendritic cells 1': '#575761',
        'Dendritic cells 2': '#1BE7FF',
        'CD4 T cells': '#1D3461',
        'CD8 T cells': '#3D5A6C',
        'NK cells': '#64B6AC'
}

## Load data

scplus_obj = dill.load(
        open(os.path.join(work_dir, 'scplus_obj.pkl'), 'rb')
)

from tqdm import tqdm
from random import sample
scplus_obj.uns['search_space'] = scplus_obj.uns['search_space'].explode('Distance')

scplus_obj.uns['search_space']['abs_distance'] = abs(scplus_obj.uns['search_space']['Distance'])

from scipy.stats import rankdata
def get_rank_region_gene(region, gene):
        search_space_region = scplus_obj.uns['search_space'].query('Name == @region')
        ranking = rankdata(search_space_region['abs_distance'], method = 'dense')
        return ranking[np.where(search_space_region['Gene'] == gene)[0][0]]

from scenicplus.utils import *
selected_regions = list(set(flatten_list([scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'][x] for x in   [x for x in scplus_obj.uns['selected_eRegulon']['Region_based'] if '_+_' in x]])))
selected_genes = list(set(flatten_list([scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'][x] for x in   [x for x in scplus_obj.uns['selected_eRegulon']['Gene_based'] if '_+_' in x]])))

n_regions_per_gene = scplus_obj.uns['eRegulon_metadata_filtered'].query('Gene in @selected_genes').query('Region in @selected_regions')[['Gene', 'Region']].drop_duplicates().groupby('Gene').count()['Region'].to_numpy()
nth_gene_region = []
for region, gene in tqdm(list(scplus_obj.uns['eRegulon_metadata_filtered'].query('Gene in @selected_genes').query('Region in @selected_regions')[['Region', 'Gene']].to_numpy()), total = 161180):
        nth_gene_region.append(get_rank_region_gene(region, gene))

import seaborn as sns
fig, axs = plt.subplots(ncols = 2, figsize = (3.75, 2.5), sharey = True)
ax = axs[0]
ax.hist(
        n_regions_per_gene,
        color = 'darkgray',
        edgecolor = 'black',
        bins = np.arange(0, 30, 1))
sns.despine(ax = ax)
ax.set_title('Number of regions linked to gene.', fontsize = 6)
ax.set_ylabel('Frequency', fontsize = 6)
ax.set_xticks(np.arange(1, 30, 5))
ax.spines['left'].set_position('zero')
ax.set_yticks([500, 1000, 1500, 2000])
ax.tick_params(axis='both', which='major', labelsize=5)
ax.set_xlim(left=0)
ax = axs[1]
ax.hist(nth_gene_region, bins = np.arange(max(nth_gene_region)) - 0.5, color = 'darkgray', edgecolor = 'black')
ax.set_xticks(np.arange(1, 30, 5))
sns.despine(ax = ax)
ax.set_title('nth. gene per region', fontsize = 6)
#ax.set_ylabel('Frequency', fontsize = 5)
ax.set_xlim(left=0)
ax.tick_params(axis='both', which='major', labelsize=5)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'n_regions_per_gene_nth_regions_per_gene.pdf'))

```


## Main Figure 2, Pannel d: DotPlot

```r

work_dir <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_granulocyte_sorted_10k/scenicplus_clustered_consensus'
plot_dir <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_pbmc/v3_clustered_consensus/plots'

# load data
library(SCopeLoomR)
loom <- open_loom(paste0(work_dir, '/SCENIC+_gene_based.loom'))
cells_AUC_gene <- AUCell::getAUC(get_regulons_AUC(loom, column.attr.name='RegulonsAUC'))
rownames(cells_AUC_gene) <- gsub('_extended', '', rownames(cells_AUC_gene) )
cell_data <- get_cell_annotation(loom)
dgem <- get_dgem(loom)

loom <- open_loom(paste0(work_dir, '/SCENIC+_region_based.loom'))
cells_AUC_region <- AUCell::getAUC(get_regulons_AUC(loom, column.attr.name='RegulonsAUC'))
rownames(cells_AUC_region) <- gsub('_extended', '', rownames(cells_AUC_region))

cells_AUC_gene <- cells_AUC_gene[rownames(cells_AUC_region), which(colnames(cells_AUC_gene) %in% colnames(cells_AUC_region))]
cells_AUC_region <- cells_AUC_region[rownames(cells_AUC_gene), colnames(cells_AUC_gene)]

selected_eRegulon_names = scan(paste0(work_dir, '/selected_eRegulons.txt'), what="", sep="\n")

# Get eRegulons
loom <- open_loom(paste0(work_dir, '/SCENIC+_gene_based.loom'))
regulons_gene <- get_regulons(loom, column.attr.name='Regulons')
regulon_gene_list <- list()
for (row in rownames(regulons_gene)){
  regulon_gene_list[[row]] <- colnames(regulons_gene)[which(regulons_gene[row,] == 1)]
}
regulon_gene_list <- regulon_gene_list[selected_eRegulon_names]
names(regulon_gene_list) <- gsub('_extended', '', names(regulon_gene_list))

regulon_gene_list <- tapply(unlist(regulon_gene_list , use.names = FALSE), rep(names(regulon_gene_list), lengths(regulon_gene_list)), FUN = c)

positive <- names(regulon_gene_list)[grep('_+', names(regulon_gene_list), fixed=TRUE)]
negative <- names(regulon_gene_list)[grep('_-', names(regulon_gene_list), fixed=TRUE)]

# calculate RSS
library(SCENIC)
rss_values <- calcRSS(cells_AUC_gene, cell_data$ACC_celltype)
rss_values <- sweep(rss_values,2,colSums(rss_values),`/`)*100
rss_values <- rss_values[,sort(colnames(rss_values))]

# Prepare data for plotting

expression_list <- list()
for (x in unique(cell_data$ACC_celltype)){
  print(x)
  expression_list[[x]] <- as.data.frame(t(log(rowSums(dgem[,grep(x, cell_data$ACC_celltype, fixed = TRUE)])/sum(rowSums(dgem[,grep(x, cell_data$ACC_celltype, fixed = TRUE)]))*10^6+1)))
}
exp_mat <- t(data.table::rbindlist(expression_list))
colnames(exp_mat) <- names(expression_list)
exp_mat <- exp_mat

sel_rel <- c(positive, negative)
order_list <- list()
for (i in 1:ncol(rss_values)){
  order_list[[i]] <- vector()
}
for (x in sel_rel){
  i <- which.max(rss_values[x,])
  order_list[[i]] <- c(order_list[[i]], x)
}
sel_rel <- unlist(order_list)

rss_values <- t(apply(rss_values, 1, function(x)(x-min(x))/(max(x)-min(x))))
rel_list <- sel_rel
rel_data <- data.frame()
for (rel in rel_list){
  for (name in colnames(rss_values)){
      tf <- strsplit(rel, split = "_")[[1]][1]
      row_to_add <- t(c(rel, name, exp_mat[tf, name], rss_values[rel,name]))
      rel_data <- rbind(rel_data, row_to_add)
   }
}
colnames(rel_data) <- c('Regulon', 'Cell_type', 'Expression', 'RSS')
sel_rel_color <- rep('grey',length(sel_rel))
sel_rel_color[which(rev(sel_rel) %in% positive)] <- 'forestgreen'
sel_rel_color[which(rev(sel_rel) %in% negative)] <- 'red'

#plot
library(ggplot2)
rel_data$Regulon <- factor(rel_data$Regulon, levels=rev(sel_rel))
rel_data$Cell_type <- gsub('_', ' ', rel_data$Cell_type)
rel_data$Cell_type <- factor(rel_data$Cell_type, levels= gsub('_', ' ', colnames(rss_values)))
rel_data$RSS <- as.numeric(rel_data$RSS)
rel_data$Expression <- as.numeric(rel_data$Expression)
rel_data$Expression_scale <- ave(rel_data$Expression, rel_data$Regulon, FUN = scale)
rel_data$Repressor_activator <- unlist(lapply(grepl("+", rel_data$Regulon, fixed=TRUE), function(x) if (x) 'Activators' else 'Repressors'))
g <- ggplot(data = rel_data, mapping = aes_string(y = 'Cell_type', x = 'Regulon')) +
    facet_grid(rows = vars(Repressor_activator), scales = "free_y", space = 'free') + 
    geom_tile(mapping = aes_string(fill = 'Expression_scale')) + 
    geom_point(mapping = aes_string(size = 'RSS'), colour="black",pch=21, fill = 'black') +
    scale_radius(range = c(0.1, 3), limits=c(0, 1)) +
    scale_fill_distiller(palette = "RdYlBu") +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip() +
    theme(legend.position = "none")
ggsave(paste0(plot_dir, '/R_Dotplot_RSS_no_legend.pdf'), g, width=3.5, height=13)

#plot
library(ggplot2)
rel_data$Regulon <- factor(rel_data$Regulon, levels=rev(sel_rel))
rel_data$Cell_type <- gsub('_', ' ', rel_data$Cell_type)
rel_data$Cell_type <- factor(rel_data$Cell_type, levels= gsub('_', ' ', colnames(rss_values)))
rel_data$RSS <- as.numeric(rel_data$RSS)
rel_data$Expression <- as.numeric(rel_data$Expression)
rel_data$Expression_scale <- ave(rel_data$Expression, rel_data$Regulon, FUN = scale)
rel_data$Repressor_activator <- unlist(lapply(grepl("+", rel_data$Regulon, fixed=TRUE), function(x) if (x) 'Activators' else 'Repressors'))
g <- ggplot(data = rel_data, mapping = aes_string(y = 'Cell_type', x = 'Regulon')) +
    facet_grid(rows = vars(Repressor_activator), scales = "free_y", space = 'free') +
    geom_tile(mapping = aes_string(fill = 'Expression_scale')) +
    geom_point(mapping = aes_string(size = 'RSS'), colour="black",pch=21, fill = 'black') +
    scale_radius(range = c(0.1, 3), limits=c(0, 1)) +
    scale_fill_distiller(palette = "RdYlBu") +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip()
ggsave(paste0(plot_dir, '/R_Dotplot_RSS.pdf'), g, width=3.5, height=13)

#plot
library(ggplot2)
rel_data$Regulon <- factor(rel_data$Regulon, levels=rev(sel_rel))
rel_data$Cell_type <- gsub('_', ' ', rel_data$Cell_type)
rel_data$Cell_type <- factor(rel_data$Cell_type, levels= gsub('_', ' ', colnames(rss_values)))
rel_data$RSS <- as.numeric(rel_data$RSS)
rel_data$Expression <- as.numeric(rel_data$Expression)
rel_data$Expression_scale <- ave(rel_data$Expression, rel_data$Regulon, FUN = scale)
rel_data$Repressor_activator <- unlist(lapply(grepl("+", rel_data$Regulon, fixed=TRUE), function(x) if (x) 'Activators' else 'Repressors'))
g <- ggplot(data = rel_data[rel_data$Repressor_activator == 'Activators', ], mapping = aes_string(y = 'Cell_type', x = 'Regulon')) +
    geom_tile(mapping = aes_string(fill = 'Expression_scale')) +
    geom_point(mapping = aes_string(size = 'RSS'), colour="black",pch=21, fill = 'black') +
    scale_radius(range = c(0.1, 3), limits=c(0, 1)) +
    scale_fill_distiller(palette = "RdYlBu") +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip() +
    theme(legend.position = "none")
ggsave(paste0(plot_dir, '/R_Dotplot_RSS_no_legend_activators.pdf'), g, width=3.5, height=13)

#plot
library(ggplot2)
rel_data$Regulon <- factor(rel_data$Regulon, levels=rev(sel_rel))
rel_data$Cell_type <- gsub('_', ' ', rel_data$Cell_type)
rel_data$Cell_type <- factor(rel_data$Cell_type, levels= gsub('_', ' ', colnames(rss_values)))
rel_data$RSS <- as.numeric(rel_data$RSS)
rel_data$Expression <- as.numeric(rel_data$Expression)
rel_data$Expression_scale <- ave(rel_data$Expression, rel_data$Regulon, FUN = scale)
rel_data$Repressor_activator <- unlist(lapply(grepl("+", rel_data$Regulon, fixed=TRUE), function(x) if (x) 'Activators' else 'Repressors'))
g <- ggplot(data = rel_data[rel_data$Repressor_activator == 'Activators', ], mapping = aes_string(y = 'Cell_type', x = 'Regulon')) +
    geom_tile(mapping = aes_string(fill = 'Expression_scale')) +
    geom_point(mapping = aes_string(size = 'RSS'), colour="black",pch=21, fill = 'black') +
    scale_radius(range = c(0.1, 3), limits=c(0, 1)) +
    scale_fill_distiller(palette = "RdYlBu") +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip()
ggsave(paste0(plot_dir, '/R_Dotplot_RSS_activators.pdf'), g, width=3.5, height=13)

```

## Main Figure 2, pannel e: overlap of regions

```python

import os
import dill
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

work_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_granulocyte_sorted_10k/scenicplus_clustered_consensus'
plot_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_pbmc/v3_clustered_consensus/plots'

color_dict = {
        'B cells': '#FF5714',
        'CD14+ Monocytes': '#6EEB83',
        'FCGR3A+ Monocytes': '#65B891',
        'Dendritic cells 1': '#575761',
        'Dendritic cells 2': '#1BE7FF',
        'CD4 T cells': '#1D3461',
        'CD8 T cells': '#3D5A6C',
        'NK cells': '#64B6AC'
}

## Load data

scplus_obj = dill.load(
        open(os.path.join(work_dir, 'scplus_obj.pkl'), 'rb')
)

celltype_to_group = {
        'CD4_T_cells': 'T_cells',
        'B_cells': 'B_cells',
        'NK_cells': 'NK_cells',
        'CD8_T_cells': 'T_cells',
        'CD14+_Monocytes': 'Myeloid',
        'Dendritic_cells_1': 'Myeloid',
        'FCGR3A+_Monocytes': 'Myeloid',
        'Dendritic_cells_2': 'Myeloid'
}

scplus_obj.metadata_cell['cell_type_group'] = [celltype_to_group[ct] for ct in scplus_obj.metadata_cell['GEX_celltype']]

from scenicplus.RSS import *

regulon_specificity_scores(
        scplus_obj,
        variable = 'cell_type_group',
        auc_key = 'eRegulon_AUC_filtered',
        signature_keys = ['Region_based'],
        selected_regulons = [x for x in scplus_obj.uns['selected_eRegulon']['Region_based'] if '-' not in x],
        out_key_suffix = 'filtered')

flat_list = lambda t: [item for sublist in t for item in sublist]
selected_markers = list(set(flat_list([scplus_obj.uns['RSS']['cell_type_groupfiltered'].loc[celltype].sort_values(ascending = False).head(10).index.to_list() for celltype in scplus_obj.uns['RSS']['cell_type_groupfiltered'].index])))

from scenicplus.plotting.correlation_plot import *

region_intersetc_data, Z = jaccard_heatmap(
        scplus_obj,
        method = 'intersect',
        gene_or_region_based = 'Region_based',
        save = os.path.join(plot_dir, 'hm_overlap_regions_manuel_selection_w_cbar.pdf'),
        use_plotly = False,
        selected_regulons = selected_markers,
        signature_key = 'eRegulon_signatures_filtered',
        figsize = (10, 10), return_data = True, vmax = 0.5, cmap = 'plasma')

fig, ax = plt.subplots(figsize = (10, 10))
sns.heatmap(region_intersetc_data, ax = ax, cmap = 'plasma', vmax = 0.5, cbar=False, xticklabels = False)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'hm_overlap_regions_manuel_selection_wo_cbar.pdf'))

selected_markers = list(
	set(flat_list([scplus_obj.uns['RSS']['cell_type_groupfiltered'].loc[celltype].sort_values(ascending = False).head(5).index.to_list() for celltype in scplus_obj.uns['RSS']['cell_type_groupfiltered'].index])))
region_intersetc_data, Z = jaccard_heatmap(
        scplus_obj,
        method = 'intersect',
        gene_or_region_based = 'Region_based',
        save = os.path.join(plot_dir, 'hm_overlap_regions_manuel_selection_w_cbar_top_5.pdf'),
        use_plotly = False,
        selected_regulons = selected_markers,
        signature_key = 'eRegulon_signatures_filtered',
        figsize = (10, 10), return_data = True, vmax = 0.4, cmap = 'plasma')

fig, ax = plt.subplots(figsize = (10, 10))
sns.heatmap(region_intersetc_data, ax = ax, cmap = 'plasma', vmax = 0.4, cbar=False, xticklabels = False)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'hm_overlap_regions_manuel_selection_wo_cbar_top_5.pdf'))


clustered_jaccard_df,p_adj_df, Z = fisher_exact_test_heatmap(
	scplus_obj = scplus_obj,
	gene_or_region_based = 'Region_based',
	signature_key = 'eRegulon_signatures_filtered',
	selected_regulons = selected_markers,
	save = os.path.join(plot_dir, 'hm_overlap_regions_manuel_selection_w_cbar_top_5_fisher.pdf'),
	use_plotly = False,
	figsize = (10, 10), return_data = True, vmax = 3, cmap = 'plasma')


```

## Main Figure 2, Pannel f: ChIP-seq coverage

```python

eGRNs_of_interest = ['PAX5_+_(205r)', 'POU2AF1_+_(170r)', 'EBF1_+_(274r)', 'POU2F2_+_(489r)']

regions_POU2AF1 = set(scplus_obj.uns['eRegulon_signatures_filtered']['Region_based']['POU2AF1_+_(170r)'])
regions_POU2F2 = set(scplus_obj.uns['eRegulon_signatures_filtered']['Region_based']['POU2F2_+_(489r)'])
regions_EBF = set(scplus_obj.uns['eRegulon_signatures_filtered']['Region_based']['EBF1_+_(274r)'])
regions_PAX = set(scplus_obj.uns['eRegulon_signatures_filtered']['Region_based']['PAX5_+_(205r)'])

import pyranges as pr

summits_PAX5 = pr.read_bed('/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/DPCL/data/ChIP-seq/summits/GM12878_PAX5_ENCSR000BHD.bed')
summits_EBF1 = pr.read_bed('/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/DPCL/data/ChIP-seq/summits/GM12878_EBF1_ENCSR000BGU.bed')
summits_POU2F2 = pr.read_bed('/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_pbmc/v3_clustered_consensus/data/ENCFF934JFA.bed.gz')

import pyBigWig
bw_PAX5 = pyBigWig.open('/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/DPCL/data/ChIP-seq/bigwig/GM12878_PAX5_ENCFF702MTT.bigWig')
bw_EBF1 = pyBigWig.open('/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/DPCL/data/ChIP-seq/bigwig/GM12878_EBF1_ENCFF107LDM.bigWig')
bw_POU2F2 = pyBigWig.open('/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_pbmc/v3_clustered_consensus/data/ENCFF803HIP.bigWig')

def region_name_to_coord(r):
        return r.replace('-', ':').split(':')

def calc_coverage(regions, n_bins, bw, offset = 500):
        coverage = []
        for region in regions:
                chrom, start, end = region_name_to_coord(region)
                coverage.append(bw.stats(chrom, int(start) - offset, int(end) + offset, nBins = n_bins))
        return pd.DataFrame(coverage)

def calc_coverage_keep_region_names(regions, n_bins, bw, offset = 500):
	coverage = []
	for region in tqdm(regions, total = len(regions)):
		chrom, start, end = region_name_to_coord(region)
		coverage.append(bw.stats(chrom, int(start) - offset, int(end) + offset, nBins = n_bins))
	return pd.DataFrame(coverage, index = regions)

from scenicplus.utils import region_names_to_coordinates
from scenicplus.utils import coord_to_region_names
from pycistarget.utils import target_to_query
def intersect_w_summit(regions, summits):
        def _extend(r):
                chrom, start, end = r.replace('-', ':').split(':')
                return f'{chrom}:{int(start) - 249}-{int(end) + 250}'
        pr_regions = pr.PyRanges(region_names_to_coordinates(regions))
        target_query = target_to_query(summits, pr_regions).drop_duplicates(subset = 'Query', keep = 'first')
        overlapping_regions = [_extend(r) for r in target_query['Target']]
        non_overlapping_regions = set(regions) - set(target_query['Query'])
        return [*overlapping_regions, *non_overlapping_regions]

n_bins = 50

region_sets_of_interest = {'PAX5': regions_PAX, 'EBF1': regions_EBF, 'POU2F2': regions_POU2F2}
bw_of_interest = {'PAX5':bw_PAX5, 'EBF1':bw_EBF1 , 'POU2F2':bw_POU2F2}
summits_of_interest = {'PAX5': summits_PAX5, 'EBF1': summits_EBF1, 'POU2F2': summits_POU2F2}
eRegulons_of_interest = set(['PAX5', 'EBF1', 'POU2F2'])

coverages = {}
coverages_w_region_names = {}
from itertools import combinations
for ChIP_seq_track in eRegulons_of_interest:
	print(f'ChIP-seq: {ChIP_seq_track}')
	for n_region_sets in range(1, len(eRegulons_of_interest) + 1):
		for combination in combinations(eRegulons_of_interest, n_region_sets):
			others = eRegulons_of_interest - set(combination)
			regions = region_sets_of_interest[combination[0]]
			if n_region_sets > 1:
				for partner in combination[1:]:
					regions = regions & region_sets_of_interest[partner]
			if len(others) > 0:
				for other in others:
					regions = regions - region_sets_of_interest[other]
			print(f'c_{ChIP_seq_track}_r_{"_n_".join(combination)}_d_{"_".join(others)}: {len(regions)}')
			if len(regions) > 0:
				coverages[f'c_{ChIP_seq_track}_r_{"_n_".join(combination)}_d_{"_".join(others)}'] = calc_coverage(
					intersect_w_summit(regions, summits_of_interest[ChIP_seq_track]), n_bins, bw_of_interest[ChIP_seq_track])
			
color_dict = {
	'EBF1': 'Blue',
	'PAX5': 'Green',
	'POU2AF1': 'Red'
}

#set global min max for scaling
min_max_factor = {}
for ChIP_seq_track in eRegulons_of_interest:
        maxs = []
        mins = []
        for x in coverages.keys():
                if x.startswith(f'c_{ChIP_seq_track}') and len(x.split('_d_')[1]) > 0:
                        print(f'{ChIP_seq_track} -- {x}')
                        df = pd.DataFrame(coverages[x])
                        mins.append(df.mean(0).min())
                        maxs.append(df.mean(0).max())
        min_max_factor[ChIP_seq_track] = [min(mins), max(maxs)]

fig, axs = plt.subplots(ncols = 3, nrows = 3, figsize = (3.75, 5.75))
for j, region_set_1 in enumerate(eRegulons_of_interest):
        axs[0, j].set_title(region_set_1)
        for i, region_set_2 in enumerate(eRegulons_of_interest):
                axs[i, 0].set_ylabel(region_set_2)
                if j <= i: #this will only plot one triangle
                        for ChIP_seq_track in eRegulons_of_interest:
                                ChIP_min, ChIP_max = min_max_factor[ChIP_seq_track]
                                data_name = [
                                        x for x in coverages.keys()
                                        if x.startswith(f'c_{ChIP_seq_track}')
                                        and region_set_1 in x.split('_r_')[1].split('_d_')[0].split('_n_')
                                        and region_set_2 in x.split('_r_')[1].split('_d_')[0].split('_n_')][0]
                                print(f'{i}, {j}: {data_name} region 1: {region_set_1} region 2: {region_set_2}')
                                coverage = coverages[data_name]
                                axs[i, j].plot(
                                        np.arange(n_bins),
                                        (coverage.mean(0) - ChIP_min) / (ChIP_max - ChIP_min), #scale between global 1 and 0
                                        lw = 1,
                                        label = ChIP_seq_track,
                                        color = color_dict[ChIP_seq_track]
                                )
                                sns.despine(ax = axs[i, j])
                                axs[i,j].set_ylim(0,1)
                                axs[i, j].set_yticks([0, 0.25, 0.5, 0.75,  1])
                                axs[i,j].set_xticks(np.arange(0, 51, 16.66666))
                                axs[i, j].set_xticklabels(np.round((np.arange(0, 51, 16.66666) * 30 - 249.5)).astype(int), fontsize = 5)
                                axs[i, j].set_axisbelow(True)
                                axs[i, j].grid()
                else:
                        sns.despine(ax = axs[i, j], top = True, bottom = True, left = True, right = True)
                        axs[i, j].set_xticks([])
                        axs[i, j].set_yticks([])
handles, labels = axs[0,0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center right', frameon = False)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'aggreation_chip_combinations_POU2AF1_EBF1_PAX5.pdf'))

```

## Supplementary Figure S9, pannel m: Accessibility in B Cells


```python

B_cell_bw = pyBigWig.open('/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_granulocyte_sorted_10k/pycistopic/consensus_peak_calling/pseudobulk_bw_files/B_cells.bw')

coverages = {}
for n_region_sets in range(1, len(eRegulons_of_interest) + 1):
	for combination in combinations(eRegulons_of_interest, n_region_sets):
		others = eRegulons_of_interest - set(combination)
		regions = region_sets_of_interest[combination[0]]
		if n_region_sets > 1:
			for partner in combination[1:]:
				regions = regions & region_sets_of_interest[partner]
		if len(others) > 0:
			for other in others:
				regions = regions - region_sets_of_interest[other]
		print(f'c_{ChIP_seq_track}_r_{"_n_".join(combination)}_d_{"_".join(others)}: {len(regions)}')
		if len(regions) > 5:
			coverages[f'c_B_cell_r_{"_n_".join(combination)}_d_{"_".join(others)}'] = calc_coverage(regions, n_bins, B_cell_bw)

fig, ax = plt.subplots(figsize = (3.5, 3.5))
for region_combination in coverages.keys():
	region_combination_name = region_combination.split('_r_')[1]
	region_combination_coverage = coverages[region_combination]
	ax.plot(
		np.arange(n_bins), 
		region_combination_coverage.mean(0), 
		label = region_combination_name.split('_d_')[0],
		color = color_dict[region_combination_name],
		lw = 1)
a = list(pd.concat([coverages[x] for x in coverages.keys() if '_n_' in x], axis = 0)[24])
b = list(pd.concat([coverages[x] for x in coverages.keys() if '_n_' not in x], axis = 0)[24])
pval = round(mannwhitneyu(a, b, nan_policy = 'omit', alternative = 'greater').pvalue, 4)
ax.text(1, 5.2 , f'p = {pval}', fontsize = 8)
ax.set_ylabel('B cells Accessibility')
ax.set_xticks(np.arange(0, 51, 16.66666))
ax.set_xticklabels((np.arange(0, 51, 16.66666) * 30 - 249.5).astype(int))
sns.despine(ax = ax)
plt.legend(frameon = False,  ncol = 1, loc = 'best', columnspacing = 0.5, fontsize = 5)
plt.grid(which = 'both')
ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'ACC_B_cell.pdf'))

fig, ax = plt.subplots(figsize = (3.5, 3.5))
ax.plot(
	np.arange(n_bins),
	pd.concat([coverages[x] for x in coverages.keys() if '_n_' in x]).mean(0),
	label = '> 1',
	color = 'blue',
	lw = 1)
ax.plot(
        np.arange(n_bins),
        pd.concat([coverages[x] for x in coverages.keys() if '_n_' not in x]).mean(0),
        label = '< 1',
        color = 'red',
        lw = 1)
ax.set_ylabel('B cells Accessibility')
ax.set_xticks(np.arange(0, 51, 16.66666))
ax.set_xticklabels((np.arange(0, 51, 16.66666) * 30 - 249.5).astype(int))
sns.despine(ax = ax)
plt.legend(frameon = False,  ncol = 1, loc = 'best', columnspacing = 0.5, fontsize = 5)
plt.grid(which = 'both')
ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'ACC_B_cell_1_v_2.pdf'))

```

## Supplementary Figure S9, pannel L: Cell type specificity

```python

expr_thr = 0
B_cells = scplus_obj.metadata_cell.query('GEX_celltype == "B_cells"').index.to_list()
X_B_cells = scplus_obj.to_df('EXP').loc[B_cells]

genes_PAX5 = set(scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based']['PAX5_+_(132g)'])
genes_POU2AF1 = set(scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based']['POU2AF1_+_(123g)'])
genes_EBF1 = set(scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based']['EBF1_+_(161g)'])

gene_sets_of_interest = {'PAX5': genes_PAX5, 'EBF1': genes_EBF1, 'POU2AF1': genes_POU2AF1}
eRegulons_of_interest = set(['PAX5', 'EBF1', 'POU2AF1'])

from itertools import combinations
gene_sets = {}
for n_gene_sets in range(1, len(eRegulons_of_interest) + 1):
        for combination in combinations(eRegulons_of_interest, n_gene_sets):
                others = eRegulons_of_interest - set(combination)
                genes = gene_sets_of_interest[combination[0]]
                if n_gene_sets > 1:
                        for partner in combination[1:]:
                                genes = genes & gene_sets_of_interest[partner]
                if len(others) > 0:
                        for other in others:
                                genes = genes - gene_sets_of_interest[other]
                	print(f'{"_n_".join(combination)}_d_{"_".join(others)}: {len(genes)}')
                	if len(genes) > 5:
                        	gene_sets[f'{"_n_".join(combination)}_d_{"_".join(others)}'] = genes
		else:
			print(f'{"_n_".join(combination)}: {len(genes)}')
			gene_sets[f'{"_n_".join(combination)}'] = genes


import tspex
X = scplus_obj.to_df('EXP')
ct_indeces = {ct: scplus_obj.metadata_cell.query('GEX_celltype == @ct').index for ct in set(scplus_obj.metadata_cell['GEX_celltype'])}
X_bulk = {}
for ct in ct_indeces:
        X_bulk[ct] = np.log1p((X.loc[ct_indeces[ct]].T / X.loc[ct_indeces[ct]].T.sum() *1e6).mean(1))
X_bulk=pd.DataFrame(X_bulk)

flatten_list = lambda t: [item for sublist in t for item in sublist]
scale = lambda X: [(x - min(X)) / (max(X) - min(X)) for x in X]
spm = tspex.TissueSpecificity(X_bulk.loc[flatten_list(gene_sets.values())], 'tau', log = False, transform = False).tissue_specificity
spm.update(pd.Series(scale(spm.values), index = spm.index))

tau_gene_sets = {combo: spm.loc[list(gene_sets[combo])] for combo in gene_sets.keys()}
data = pd.DataFrame(
        [(np.repeat(combo, len(tau_gene_sets[combo])), tau_gene_sets[combo]) for combo in gene_sets.keys()],
        columns = ['Gene_set', 'tau']).explode(['Gene_set', 'tau']).reset_index(drop = True)
data['label'] = [x.split('_d_')[0].replace('_n_', '\n') for x in data['Gene_set']]
data['label'] = data['label'].astype('category')
data['tau'] = data['tau'].astype(float)

color_dict['EBF1_n_PAX5_n_POU2AF1'] = '#4CB963'
color_dict_genes = {x.split('_d_')[0].replace('_n_', '\n'): color_dict[x] for x in color_dict.keys()}

data['number_of_factors'] = ['<2' if '\n' not in x else '>=2' for x in data['label']]

a = data.query('number_of_factors == "<2"')['tau'].tolist()
b = data.query('number_of_factors == ">=2"')['tau'].tolist()
pval = round(mannwhitneyu(b, a, nan_policy = 'omit', alternative = 'greater').pvalue, 4)

fig, ax = plt.subplots(figsize = (3.5, 3.5))
sns.boxplot(
	ax = ax, 
	y = 'label', 
	x = 'tau', 
	data = data, 
	order = [x.split('_d_')[0].replace('_n_', '\n') for x in gene_sets.keys()], 
	palette = color_dict_genes)
ax.text(0.9, 3.5, f'p = {pval}', fontsize = 5)
sns.despine(ax = ax)
ax.set_ylabel(None)
plt.grid()
ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'tau.pdf'))

fig, ax = plt.subplots(figsize = (3.75, 3.5))
sns.boxplot(
        ax = ax,
	y = 'number_of_factors',
	x = 'tau',
	data = data,
	palette = {'<2': 'red', '>=2': 'blue'})
sns.despine(ax = ax)
ax.set_ylabel(None)
plt.grid()
ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'tau_1_2_.pdf'))

```

## Export Region-to-gene links for visualization inR

```python

import os
import dill

work_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_granulocyte_sorted_10k/scenicplus_clustered_consensus'

## Load data

scplus_obj = dill.load(
        open(os.path.join(work_dir, 'scplus_obj.pkl'), 'rb')
)

with open(os.path.join(work_dir, 'scplus_cell_names.txt'), 'w') as f:
	for cell_barcode in scplus_obj.cell_names:
		f.write(cell_barcode.split('-10x_pbmc')[0])
		f.write('\n')

with open(os.path.join(work_dir, 'scplus_region_names.txt'), 'w') as f:
	for region_name in scplus_obj.region_names:
		f.write(region_name.replace(':', '-'))
		f.write('\n')

with open(os.path.join(work_dir, 'scplus_gene_names.txt'), 'w') as f:
	for gene_name in scplus_obj.gene_names:
		f.write(gene_name)
		f.write('\n')

scplus_obj.metadata_cell.to_csv(os.path.join(work_dir, 'scplus_metadata_cell.tsv'), sep = '\t')

save_path = work_dir
path_bedToBigBed = '/staging/leuven/stg_00002/lcb/sdewin/PhD/De_Winter_hNTorg/COMBINED_ANALYSIS/r2g/'
biomart_host = 'http://sep2019.archive.ensembl.org'
assembly = 'hg38'
species = 'hsapiens'
from scenicplus.enhancer_to_gene import export_to_UCSC_interact 
r2g_data = export_to_UCSC_interact(scplus_obj,
		    species,
		    os.path.join(save_path,'r2g.rho.bed'),
		    path_bedToBigBed=path_bedToBigBed,
		    bigbed_outfile=os.path.join(save_path,'r2g.rho.bb'),
		    region_to_gene_key='region_to_gene',
		    pbm_host=biomart_host,
		    assembly=assembly,
		    ucsc_track_name='R2G',
		    ucsc_description='SCENIC+ region to gene links',
		    cmap_neg='Reds',
		    cmap_pos='Greens',
		    key_for_color='rho',
		    scale_by_gene=False,
		    subset_for_eRegulons_regions=True,
		    eRegulons_key='eRegulons')
r2g_data = export_to_UCSC_interact(scplus_obj,
		    species,
		    os.path.join(save_path,'r2g.importance.bed'),
		    path_bedToBigBed=path_bedToBigBed,
		    bigbed_outfile=os.path.join(save_path,'r2g.importance.bb'),
		    region_to_gene_key='region_to_gene',
		    pbm_host=biomart_host,
		    assembly=assembly,
		    ucsc_track_name='R2G',
		    ucsc_description='SCENIC+ region to gene links',
		    cmap_neg='Reds',
		    cmap_pos='Greens',
		    key_for_color='importance',
		    scale_by_gene=True,
		    subset_for_eRegulons_regions=True,
		    eRegulons_key='eRegulons')

from scenicplus.utils import export_eRegulons
regions = export_eRegulons(scplus_obj,
	os.path.join(save_path,'eRegulons.bed'),
	assembly,
	eRegulon_signature_key = 'eRegulon_signatures_filtered',
	bigbed_outfile = os.path.join(save_path,'eRegulons.bb'),
	path_bedToBigBed=path_bedToBigBed)

```

## Main Figure 2, pannel g and Supplementary Figure S9, pannel k: Genome Browser plots

```r

# create signac object
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)

counts <- Read10X_h5('/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_combined/data/pbmc_granulocyte_sorted_10k/outs/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5')
fragpath <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_combined/data/pbmc_granulocyte_sorted_10k/outs/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz'
work_dir <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_granulocyte_sorted_10k/scenicplus_clustered_consensus'
plot_dir <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_pbmc/v3_clustered_consensus/plots'


annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

# load data
library(SCopeLoomR)
loom <- open_loom(paste0(work_dir, '/SCENIC+_gene_based.loom'))
cell_data <- get_cell_annotation(loom)

cells_to_keep <- gsub('-10x_pbmc', '', rownames(cell_data))

pbmc <- subset(pbmc, cells = cells_to_keep)

r2g_links <- read.table(
	paste0(work_dir, '/r2g.rho.bed'), sep = '\t', skip = 1,
	col.names = c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'value', 'exp',
       		      'color', 'sourceChrom', 'sourceStart', 'sourceEnd', 'sourceName',
       	              'sourceStrand', 'targetChrom', 'targetStart', 'targetEnd', 'targetName',
       	              'targetStrand'))

cell_metadata <- read.table(paste0(work_dir, '/scplus_metadata_cell.tsv'), sep = '\t', header = T)
rownames(cell_metadata) = cell_metadata$ACC_barcode
pbmc <- AddMetaData(pbmc, cell_metadata)

DefaultAssay(pbmc) <- 'ATAC'

library(GenomicRanges)
r2g_links <- transform(r2g_links, start_ = pmin(sourceStart, targetEnd))
r2g_links <- transform(r2g_links, end_ = pmax(sourceStart, targetEnd))
r2g_links <- r2g_links[c('chrom', 'start_', 'end_', 'targetName', 'value', 'chrom', 'chromStart', 'chromEnd')]
colnames(r2g_links) <- c('seqnames', 'start', 'end', 'gene', 'score', 'chrom', 'chromStart', 'chromEnd')

eRegulon_regions <- read.table(paste0(work_dir, '/eRegulons.bed'), sep = '\t')
colnames(eRegulon_regions) <- c('seqnames', 'start', 'end', 'name')


color_dict = c(
        'B_cells' = '#FF5714',
        'CD14+_Monocytes' = '#6EEB83',
        'FCGR3A+_Monocytes' = '#65B891',
        'Dendritic_cells_1' = '#575761',
        'Dendritic_cells_2' = '#1BE7FF',
        'CD4_T_cells' = '#1D3461',
        'CD8_T_cells' = '#3D5A6C',
        'NK_cells' = '#64B6AC')

selected_eRegulons <- c('POU2AF1_+', 'EBF1_+', 'PAX5_+')

color_dict_eRegulons = c(
	'POU2AF1_+' = 'Red',  
	'EBF1_+' = 'Green',  
	'PAX5_+' = 'Blue',
	'EBF1_+, POU2AF1_+' = 'goldenrod1',
	'EBF1_+, PAX5_+' =  'cyan4',
	'PAX5_+, POU2AF1_+' = 'purple4',
	'EBF1_+, PAX5_+, POU2AF1_+' = 'black'
)

bw_PAX5 <- '/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/DPCL/data/ChIP-seq/bigwig/GM12878_PAX5_ENCFF702MTT.bigWig'
bw_EBF1 <- '/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/DPCL/data/ChIP-seq/bigwig/GM12878_EBF1_ENCFF107LDM.bigWig'
bw_POU2F2 <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_pbmc/v3_clustered_consensus/data/ENCFF803HIP.bigWig'

library(dplyr)
eRegulon_regions <- subset(eRegulon_regions, name %in% selected_eRegulons)
eRegulon_regions <- eRegulon_regions[order(eRegulon_regions$name), ]
eRegulon_regions <- eRegulon_regions %>%
                        group_by(seqnames, start, end) %>%
                        mutate(name = paste0(name, collapse = ", "))
gr_eRegulon_regions <- makeGRangesFromDataFrame(eRegulon_regions, keep.extra.columns = TRUE)
r2g_links_gene <- subset(r2g_links, gene == "BLNK")
r2g_links_gr <- makeGRangesFromDataFrame(subset(r2g_links_gene, chromStart %in% eRegulon_regions$start)[c('seqnames', 'start', 'end', 'gene', 'score')], keep.extra.columns = TRUE)
Links(pbmc) <- r2g_links_gr

gene_to_region_of_interest = c(
	'BLNK' = 'chr10-96226082-96316945',
	'PLEKHG1' = 'chr6-150562217-150730918',
	'CD22' ='chr19-35316521-35485142'
)

library(ggplot2)
library(RColorBrewer)

for(gene in names(gene_to_region_of_interest)) {
	print(gene)
	region = gene_to_region_of_interest[[gene]]
	eRegulon_regions <- eRegulon_regions[eRegulon_regions[['name']] %in% selected_eRegulons, ]
	eRegulon_regions <- eRegulon_regions[order(eRegulon_regions$name), ]
	eRegulon_regions <- eRegulon_regions %>%
				group_by(seqnames, start, end) %>%
				mutate(name = paste0(name, collapse = ", "))
	gr_eRegulon_regions <- makeGRangesFromDataFrame(eRegulon_regions, keep.extra.columns = TRUE)
	r2g_links_gene <- r2g_links[r2g_links[['gene']] == gene, ]
	r2g_links_gr <- makeGRangesFromDataFrame(r2g_links_gene[r2g_links_gene[['chromStart']] %in% eRegulon_regions$start,][c('seqnames', 'start', 'end', 'gene', 'score')], keep.extra.columns = TRUE)
	Links(pbmc) <- r2g_links_gr
	peak_plot <- PeakPlot(
		object = pbmc,
		region = region,
		peaks = gr_eRegulon_regions,
		group.by = 'name'
	) & scale_fill_manual(values=color_dict_eRegulons)
	cov_plot <- CoveragePlot(
	  object = pbmc,
	  region = region,
	  group.by = 'GEX_celltype',
	  annotation = FALSE,
	  peaks = FALSE,
	  links=FALSE,
	  bigwig = list('PAX5 ChIP-seq' = bw_PAX5, 'EBF1 ChIP-seq' = bw_EBF1, 'POU2F2 ChIP-seq' = bw_POU2F2),
	  bigwig.scale = 'separate'
	)  & scale_fill_manual(values=color_dict)
	cov_plot_wo_bw <- CoveragePlot(
          object = pbmc,
          region = region,
          group.by = 'GEX_celltype',
          annotation = FALSE,
          peaks = FALSE,
          links=FALSE
        )  & scale_fill_manual(values=color_dict)
	expr_plot <- ExpressionPlot(
		object = pbmc,
		features = c(gene, 'PAX5', 'POU2AF1', 'EBF1'),
		assay = 'RNA',
		group.by = 'GEX_celltype'
	) & scale_fill_manual(values=color_dict)
	# Annotation
	gene_plot <- AnnotationPlot(
	  object = pbmc,
	  region = region
	)
	# Links
	link_plot <- LinkPlot(
	  object = pbmc,
	  region = region,
	  min.cutoff = -1) + scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Blues"), limits=c(0, NA))
	library(patchwork)
	p <- CombineTracks(
	  plotlist = list(cov_plot, peak_plot, link_plot, gene_plot),
	  heights = c(5, 1, 1, 1),
	  widths = c(10,4)
	)
	pdf(paste0(plot_dir, paste0('/R_R2G_', gene, '.pdf')))
	print(p)
	dev.off()	
	p <- CombineTracks(
		plotlist = list(cov_plot_wo_bw, peak_plot, link_plot, gene_plot),
		expression.plot = expr_plot,
		heights = c(5, 1, 1, 1),
		widths = c(10,4))
	pdf(paste0(plot_dir, paste0('/R_R2G_', gene, '_w_expr.pdf')))
	print(p)
	dev.off()
}

saveRDS(pbmc, paste0(work_dir, 'SIGNAC_obj.RDS'))

```

## Supplementary Figure S9, pannel c: Comparison with HiC

```python

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import os

wdir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/figures_reviewer_questions/HiC_pbmc'
data_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/HIC_PBMC'

merged_data_hic = pd.read_csv(
	os.path.join(data_dir, 'merged_data_hic.tsv'), sep = '\t')

merged_data_scplus = pd.read_csv(
	os.path.join(data_dir, 'merged_data_scplus.tsv'), sep = '\t')

import dill
scplus_obj = dill.load(
        open('/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_granulocyte_sorted_10k/scenicplus_clustered_consensus/scplus_obj.pkl', 'rb'))

from scenicplus.utils import *
selected_regions = list(set(flatten_list([scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'][x] for x in   [x for x in scplus_obj.uns['selected_eRegulon']['Region_based'] if '_+_' in x]])))
selected_genes = list(set(flatten_list([scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'][x] for x in   [x for x in scplus_obj.uns['selected_eRegulon']['Gene_based'] if '_+_' in x]])))

n_regions_per_gene = scplus_obj.uns['eRegulon_metadata_filtered'].query('Gene in @selected_genes').query('Region in @selected_regions')[['Gene', 'Region']].drop_duplicates().groupby('Gene').count()['Region'].to_numpy()


#histogram of nr of regions per gene and rank
fig, axs = plt.subplots(nrows = 2, ncols = 2, figsize = (8,8), sharex = 'col', sharey = 'col')
axs[0, 0].hist(
	merged_data_scplus.groupby('Gene')['Region'].apply(lambda x: len(set(x))),
	bins = np.arange(merged_data_scplus.groupby('Gene')['Region'].apply(lambda x: len(set(x))).max()), color = 'gray', edgecolor = 'black')
axs[1, 0].hist(
	merged_data_hic.groupby('Gene')['Region'].apply(lambda x: len(set(x))),
	bins = np.arange(merged_data_hic.groupby('Gene')['Region'].apply(lambda x: len(set(x))).max()), color = 'gray', edgecolor = 'black')
rank_max_importance_splus = merged_data_scplus.groupby('Region')[['rank_gene', 'R2G_importance']].apply(lambda x: np.array(x['rank_gene'])[np.argmax(x['R2G_importance'])])
axs[0, 1].hist(
	rank_max_importance_splus,
	bins = np.arange(rank_max_importance_splus.max()), color = 'gray', edgecolor = 'black')
rank_max_score = merged_data_hic.groupby('Region')[['rank_gene', 'Score']].apply(lambda x: np.array(x['rank_gene'])[np.argmax(x['Score'])])
axs[1, 1].hist(
	rank_max_score,
	bins = np.arange(rank_max_score.max()), color = 'gray', edgecolor = 'black')
for ax in axs.ravel():
	sns.despine(ax = ax)
axs[0, 0].set_ylabel('Frequency')
axs[1, 0].set_ylabel('Frequency')
axs[0, 0].set_title('SCENIC+: # regions per gene')
axs[0, 1].set_title('SCENIC+: nth. gene per region')
axs[1, 0].set_title('Hi-C: # regions per gene')
axs[1, 1].set_title('Hi-C: nth. gene per region')
fig.tight_layout()
fig.savefig(os.path.join(wdir, 'hist_n_regions_per_gene_and_rank.pdf'))
fig.savefig(os.path.join(wdir, 'hist_n_regions_per_gene_and_rank.png'))

hic_correlation_data = pd.read_csv(
	os.path.join(data_dir, 'correlation_w_hic_B_cell_marker_genes.tsv'), sep = '\t')

order = ['GBM_rnd', 'GBM', 'rho_rnd', 'rho']
hic_correlation_data['method'] = pd.Categorical(hic_correlation_data['method'], categories = order)
color_dict = {
	'GBM_rnd': '#D0CFEC', 
	'GBM': '#3454D1', 
	'rho_rnd': '#F8C0C8', 
	'rho': '#CC5A71'}
sns.set_palette([color_dict[t] for t in order])


from scipy.stats import mannwhitneyu
pairs = [('GBM_rnd', 'GBM'), ('rho_rnd', 'rho')]
p_values = [
	mannwhitneyu(hic_correlation_data.loc[hic_correlation_data['method'] == t1, 'correlation'], hic_correlation_data.loc[hic_correlation_data['method'] == t2, 'correlation']).pvalue
	for t1, t2 in pairs]

p_values_formated = [f'p={p:.1e}' for p in p_values]
from statannotations.Annotator import Annotator

fig, ax = plt.subplots()
sns.boxplot(
	data = hic_correlation_data,
	x = 'method', y = 'correlation')
ax.hlines(y = 0, xmin = -0.5, xmax = 4, ls = 'dashed', color = 'gray')
ax.set_ylim((-1, 1))
annotator = Annotator(ax, pairs, data = hic_correlation_data, x = 'method', y = 'correlation')
annotator.set_custom_annotations(p_values_formated)
annotator.annotate()
sns.despine(ax = ax)
fig.tight_layout()
fig.savefig(os.path.join(wdir, 'bxplot_corr_hic.pdf'))
fig.savefig(os.path.join(wdir, 'bxplot_corr_hic.png'))

sum(merged_data_scplus.groupby('Gene')['Region'].apply(lambda x: len(set(x))).values < 10) #5928
len(set(merged_data_scplus['Gene'])) #6605

sum(merged_data_hic.groupby('Gene')['Region'].apply(lambda x: len(set(x))).values < 10) #4047
len(set(merged_data_hic['Gene'])) #6825

sum(rank_max_importance_splus == 1) #11533
len(rank_max_importance_splus) #23470

sum(rank_max_score == 1) #17903
len(rank_max_score) #23450

```

## Supplementary Figure S9, pannel e: TFs found by SCENIC+, Signac and ArchR

```python

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

outdir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/figures_reviewer_questions/pbmc_signac_archr'

import os
import pandas as pd
import numpy as np

import dill
scplus_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_granulocyte_sorted_10k/scenicplus_clustered_consensus"
scplus_obj = dill.load(
        open(os.path.join(scplus_folder, 'scplus_obj.pkl'), 'rb'))

data_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/extra_methods/pbmc/extract_data"

archr_signac_outdir = '/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/archr_signac_comparison/pbmc/'

tool_to_data = {
	'SCENIC+': scplus_obj.uns['eRegulon_metadata_filtered'].copy(),
	'Signac DEM': pd.read_table(os.path.join(archr_signac_outdir, 'find_motifs_dars_matrix.tsv')),
	'Signac ChromVar': pd.read_table(os.path.join(archr_signac_outdir, 'chromVAR_cistromes.tsv')),
	'ArchR motif DEM': pd.read_table(os.path.join(archr_signac_outdir, 'ArchR/ArchR_motifs_DARs_matrix.tsv')),
	'ArchR track DEM': pd.read_table(os.path.join(archr_signac_outdir, 'ArchR/ArchR_ENCODE_DARs_matrix.tsv')),
	'ArchR ChromVar': pd.read_table(os.path.join(archr_signac_outdir, 'ArchR/ArchR_chromVAR_cistromes_matrix.tsv'))}

from tqdm import tqdm
def split_dimers_into_two_columns(df: pd.DataFrame):
	for column in tqdm(df.columns, total = len(df.columns)):
		if '::' in column:
			col1, col2 = column.split('::')
			if col1 in df.columns:
				df[col1] = df[col1].to_numpy() + df[column].to_numpy()
			else:
				df[col1] = df[column].to_numpy()
			if col2 in df.columns:
				df[col2] = df[col2].to_numpy() + df[column].to_numpy()
			else:
				df[col2] = df[column].to_numpy()
			df.drop(column, axis = 1, inplace = True)

for tool in set(tool_to_data.keys()) - set(['SCENIC+']):
	print(tool)
	split_dimers_into_two_columns(tool_to_data[tool])


tool_to_data['ArchR track DEM'].columns =  [x.split('_')[0].split('-')[0] for x in tool_to_data['ArchR track DEM'].columns]

tool_to_data['ArchR track DEM'] = tool_to_data['ArchR track DEM'].groupby(lambda x: x, axis = 1).sum()

tool_to_data['SCENIC+']['Consensus_name'] = [x.split('_(')[0] for x in tool_to_data['SCENIC+']['Region_signature_name']]

tool_to_data['SCENIC+'].rename({'Region': 'region', 'Gene': 'gene'}, axis = 1, inplace = True)

TFs_all_tools = list(set(tool_to_data['SCENIC+']['TF']))

TFs_all_tools = []
for tool in set(tool_to_data.keys()) - set(['SCENIC+']):
	TFs_all_tools.extend(list(tool_to_data[tool].columns[tool_to_data[tool].sum() >= 10]))
TFs_all_tools = set(TFs_all_tools)

TFs_for_tools = pd.DataFrame(
	index = TFs_all_tools,
	columns = tool_to_data.keys()).fillna(0)

TFs_for_tools.loc[tool_to_data['SCENIC+']['TF']] = 1 

for tool in set(tool_to_data.keys() - set(['SCENIC+'])):
	TFs_for_tools.loc[list(tool_to_data[tool].columns[tool_to_data[tool].sum() >= 10]), tool] = 1

TFs_for_tools = TFs_for_tools.loc[list(set(TFs_for_tools.index) & set(scplus_obj.gene_names))]

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.clustermap(
	TFs_for_tools.T,
	cmap = ['Red', 'Green'],
	yticklabels = True,
	xticklabels = False,
	figsize = (13, 4),
	dendrogram_ratio = 0.05,
	**dict(cbar_kws=dict(ticks=[0, 1], orientation='horizontal', shrink = 0.01)))
x0, _y0, _w, _h = g.cbar_pos
g.ax_heatmap.set_yticklabels(labels = g.ax_heatmap.get_yticklabels(), va = 'center', rotation = 0)
g.ax_heatmap.set_xticklabels(labels = g.ax_heatmap.get_xticklabels(), ha = 'right', rotation = 45, rotation_mode="anchor")
g.ax_cbar.set_position([0.02, 0.98, 0.01, 0.02])
g.ax_cbar.set_title('TF found')
#g.ax_cbar.set_ticks([0.25,0.75])
#g.ax_cbar.set_ticklabels(['0', '1'])
g.savefig(os.path.join(outdir, "scplus_TFs_found_by_methods.pdf"))
g.savefig(os.path.join(outdir, "scplus_TFs_found_by_methods.png"), dpi = 300)


TFs_for_tools_w_cistromes = TFs_for_tools.copy()

cistrome_TFs = list(set([x.split('_')[0] for x in scplus_obj.uns['Cistromes']['Unfiltered'].keys()]))

new_TFs = list(set(cistrome_TFs) - set(TFs_for_tools_w_cistromes.index))

TFs_for_tools_w_cistromes = pd.concat([pd.DataFrame(index = new_TFs, columns = TFs_for_tools_w_cistromes.columns), TFs_for_tools_w_cistromes]).fillna(0)

TFs_for_tools_w_cistromes.loc[cistrome_TFs, 'SCENIC+ Cistromes'] = 1
TFs_for_tools_w_cistromes = TFs_for_tools_w_cistromes.fillna(0)

g = sns.clustermap(
	TFs_for_tools_w_cistromes.T,
	cmap = ['Red', 'Green'],
	yticklabels = True,
	xticklabels = False,
	figsize = (13, 4),
	dendrogram_ratio = 0.05,
	**dict(cbar_kws=dict(ticks=[0, 1], orientation='horizontal', shrink = 0.01)))
x0, _y0, _w, _h = g.cbar_pos
g.ax_heatmap.set_yticklabels(labels = g.ax_heatmap.get_yticklabels(), va = 'center', rotation = 0)
g.ax_heatmap.set_xticklabels(labels = g.ax_heatmap.get_xticklabels(), ha = 'right', rotation = 45, rotation_mode="anchor")
g.ax_cbar.set_position([0.02, 0.98, 0.01, 0.02])
g.ax_cbar.set_title('TF found')
#g.ax_cbar.set_ticks([0.25,0.75])
#g.ax_cbar.set_ticklabels(['0', '1'])
g.savefig(os.path.join(outdir, "scplus_TFs_found_by_methods_w_cistromes.pdf"))
g.savefig(os.path.join(outdir, "scplus_TFs_found_by_methods_w_cistromes.png"), dpi = 300)

```

## Supplementary Figure S9, pannel f: TFs found by SCENIC+, CellOracle, FigR and SCENIC


```python

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

outdir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/figures_reviewer_questions/pbmc_other_eGRN_methods'
import os
os.makedirs(outdir)

import os
import pandas as pd
import numpy as np

import dill
scplus_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/pbmc_granulocyte_sorted_10k/scenicplus_clustered_consensus"
scplus_obj = dill.load(
        open(os.path.join(scplus_folder, 'scplus_obj.pkl'), 'rb'))

data_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/extra_methods/pbmc/extract_data"

tool_to_regulon_metadata = {
        'SCENIC': pd.read_table(os.path.join(data_folder, 'SCENIC_regulons.tsv')),
        'CellOracle': pd.read_table(os.path.join(data_folder, 'celloracle_eRegulons.tsv')),
        'FigR': pd.read_table(os.path.join(data_folder, 'figr_eRegulons.tsv')),
        'Pando': pd.read_table(os.path.join(data_folder, 'pando_eRegulons.tsv')),
        #'SCENIC+': scplus_obj.uns['eRegulon_metadata_filtered'].copy()}

tool_to_regulon_metadata['SCENIC']['TF'] = [x.split('_')[0] for x in tool_to_regulon_metadata['SCENIC']['TF']]

tool_to_regulon_metadata['SCENIC']['Consensus_name'] = tool_to_regulon_metadata['SCENIC']['TF']

tool_to_regulon_metadata['CellOracle']['Consensus_name'] = [
	f"{TF}_+" if score >= 0 else f"{TF}_-" for TF, score in tool_to_regulon_metadata['CellOracle'][['TF', 'coef_mean']].to_numpy()]

tool_to_regulon_metadata['FigR']['Consensus_name'] = [
	f"{TF}_+" if score >= 0 else f"{TF}_-" for TF, score in tool_to_regulon_metadata['FigR'][['TF', 'Corr']].to_numpy()]

tool_to_regulon_metadata['Pando']['Consensus_name'] = [
	 f"{TF}_+" if score >= 0 else f"{TF}_-" for TF, score in tool_to_regulon_metadata['Pando'][['TF', 'estimate']].to_numpy()]

tool_to_regulon_metadata['SCENIC+']['Consensus_name'] = [x.split('_(')[0] for x in tool_to_regulon_metadata['SCENIC+']['Region_signature_name']]

tool_to_regulon_metadata['SCENIC+'].rename({'Region': 'region', 'Gene': 'gene'}, axis = 1, inplace = True)

min_target_genes = 10
for tool in tool_to_regulon_metadata.keys():
	eRegulon_and_counts = tool_to_regulon_metadata[tool].groupby('Consensus_name')['gene'].apply(lambda x: len(set(x)))
	eRegulons_to_keep = eRegulon_and_counts.index[eRegulon_and_counts >= min_target_genes]
	tool_to_regulon_metadata[tool] = tool_to_regulon_metadata[tool].query('Consensus_name in @eRegulons_to_keep')

TFs_all_tools = []
for tool in tool_to_regulon_metadata.keys():
	TFs_all_tools.extend(tool_to_regulon_metadata[tool]['TF'])
TFs_all_tools = set(TFs_all_tools)

TFs_for_tools = pd.DataFrame(
	index = TFs_all_tools,
	columns = tool_to_regulon_metadata.keys()).fillna(0)
for tool in tool_to_regulon_metadata.keys():
	TFs_found_by_tool = list(set(tool_to_regulon_metadata[tool]['TF']))
	TFs_for_tools.loc[TFs_found_by_tool, tool] = 1

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.clustermap(
	TFs_for_tools.T[[x.split('_')[0] for x in scplus_obj.uns['selected_eRegulon']['Gene_based']]],
	cmap = ['Red', 'Green'],
	yticklabels = True,
	xticklabels = True,
	figsize = (13, 4),
	dendrogram_ratio = 0.05,
	**dict(cbar_kws=dict(ticks=[0, 1], orientation='horizontal', shrink = 0.01)))
x0, _y0, _w, _h = g.cbar_pos
g.ax_heatmap.set_yticklabels(labels = g.ax_heatmap.get_yticklabels(), va = 'center', rotation = 0)
g.ax_heatmap.set_xticklabels(labels = g.ax_heatmap.get_xticklabels(), ha = 'right', rotation = 45, rotation_mode="anchor")
g.ax_cbar.set_position([0.02, 0.98, 0.01, 0.02])
g.ax_cbar.set_title('TF found')
#g.ax_cbar.set_ticks([0.25,0.75])
#g.ax_cbar.set_ticklabels(['0', '1'])
g.savefig(os.path.join(outdir, "scplus_TFs_found_by_methods.pdf"))
g.savefig(os.path.join(outdir, "scplus_TFs_found_by_methods.png"), dpi = 300)


g = sns.clustermap(
	TFs_for_tools.T,
	cmap = ['Red', 'Green'],
	yticklabels = True,
	xticklabels = False,
	figsize = (13, 4),
	dendrogram_ratio = 0.05,
	**dict(cbar_kws=dict(ticks=[0, 1], orientation='horizontal', shrink = 0.01)))
x0, _y0, _w, _h = g.cbar_pos
g.ax_heatmap.set_yticklabels(labels = g.ax_heatmap.get_yticklabels(), va = 'center', rotation = 0)
g.ax_heatmap.set_xticklabels(labels = g.ax_heatmap.get_xticklabels(), ha = 'right', rotation = 45, rotation_mode="anchor")
g.ax_cbar.set_position([0.02, 0.98, 0.01, 0.02])
g.ax_cbar.set_title('TF found')
#g.ax_cbar.set_ticks([0.25,0.75])
#g.ax_cbar.set_ticklabels(['0', '1'])
g.savefig(os.path.join(outdir, "TFs_found_by_methods.pdf"))
g.savefig(os.path.join(outdir, "TFs_found_by_methods.png"), dpi = 300)

```

## Supplementary Fig S9, pannel g: GO enrichment on TFs found by Signac, SCENIC+ and ArchR

```python


import pandas as pd

out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/figures_reviewer_questions/pbmc_signac_archr'

import os

data_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/extra_methods/pbmc/pbmc_signac_archr/'
geo_result_files = {
	'SCENIC+': pd.read_csv(os.path.join(data_dir, 'SCPLUS_targets.csv')),
	'Signac': pd.read_csv(os.path.join(data_dir, 'Signac_targets.csv')),
	'ArchR':  pd.read_csv(os.path.join(data_dir, 'ArchR_targets.csv'))}

top_n = 10
top_on = 'negative_log10_of_adjusted_p_value'

geo_result = []
for tool in geo_result_files.keys():
	tmp = geo_result_files[tool]
	tmp = tmp.nlargest(top_n, top_on)
	tmp['Tool'] = tool
	geo_result.append(tmp)

geo_result = pd.concat(geo_result)

import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import seaborn as sns

fig, axs = plt.subplots(ncols = 3, figsize = (3 * 4, 5))
for tool,ax in zip(set(geo_result['Tool']), axs.ravel()):
	data = geo_result.loc[geo_result['Tool'] == tool]
	ax.scatter(
		y = np.arange(len(data))[::-1], x = data['negative_log10_of_adjusted_p_value'], color = 'black')
	texts = [ax.text(x, y, t) for x, y, t in zip(data['negative_log10_of_adjusted_p_value'], np.arange(len(data))[::-1], data['term_name'])]
	#adjust_text(texts)
	ax.set_title(f'{tool}')
	ax.set_xlabel('log10(Adjusted p value)')
	sns.despine(ax = ax)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'geo_scatter.pdf'))
fig.savefig(os.path.join(out_dir, 'geo_scatter.png'))

```

## Supplementary Fig S9, pannel h: Supplementary Fig S9, pannel g: GO enrichment on TFs found by Pando, SCENIC+ and CellOracle

```python

import pandas as pd

out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/figures_reviewer_questions/pbmc_other_eGRN_methods'

import os

data_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/extra_methods/pbmc/extract_data'
geo_result_files = {
	'SCENIC+': pd.read_csv(os.path.join(data_dir, 'gProfiler_scenic_plus.csv')),
	'CellOracle': pd.read_csv(os.path.join(data_dir, 'gProfiler_celloracle.csv')),
	'Pando':  pd.read_csv(os.path.join(data_dir, 'gProfiler_pando.csv')),
	'SCENIC': pd.read_csv(os.path.join(data_dir, 'gProfiler_scenic.csv'))}

top_n = 10
top_on = 'negative_log10_of_adjusted_p_value'

geo_result = []
for tool in geo_result_files.keys():
	tmp = geo_result_files[tool]
	tmp = tmp.nlargest(top_n, top_on)
	tmp['Tool'] = tool
	geo_result.append(tmp)

geo_result = pd.concat(geo_result)

import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import seaborn as sns

fig, axs = plt.subplots(ncols = 3, figsize = (3 * 4, 5), sharex = False)
for tool,ax in zip(set(geo_result['Tool']), axs.ravel()):
	data = geo_result.loc[geo_result['Tool'] == tool]
	ax.scatter(
		y = np.arange(len(data))[::-1], x = data['negative_log10_of_adjusted_p_value'], color = 'black')
	texts = [ax.text(x, y, t) for x, y, t in zip(data['negative_log10_of_adjusted_p_value'], np.arange(len(data))[::-1], data['term_name'])]
	#adjust_text(texts)
	ax.set_title(f'{tool}')
	ax.set_xlabel('log10(Adjusted p value)')
	sns.despine(ax = ax)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'geo_scatter.pdf'))
fig.savefig(os.path.join(out_dir, 'geo_scatter.png'))

```

## Supplementary Fig S9, pannel i: Accuracy of target regions identified by SCENIC+, Signac and ArchR


```python

data_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/extra_methods/pbmc/pbmc_signac_archr'

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

outdir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/figures_reviewer_questions/pbmc_signac_archr'
out_dir = outdir

import pandas as pd
import os
data_f1 = pd.read_csv(os.path.join(data_dir, 'data_f1.tsv'), sep = '\t')
data_pr = pd.read_csv(os.path.join(data_dir, 'data_pr.tsv'), sep = '\t')
data_re = pd.read_csv(os.path.join(data_dir, 'data_re.tsv'), sep = '\t')

tool_to_color = {
	'SCENIC+ Cistromes': '#f70404',
	'Signac ChromVar': '#3772FF',
	'Signac DEM': '#30BCED',
	'ArchR motif DEM': '#44AF69',
	'ArchR ChromVar': '#008148',
	'ArchR track DEM': '#C8D5B9'}

def plot(params, ax):
	sns.violinplot(
		inner = None,
		ax = ax,
		**params)
	sns.stripplot(
		color="#2b2b2b", size = 2,
		ax = ax,
		**params)

def set_y_log(ax):
	ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
	ymin, ymax = ax.get_ylim()
	tick_range = np.arange(np.floor(ymin), ymax)
	ax.yaxis.set_ticks(tick_range)
	ax.yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)

import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from itertools import combinations
from scenicplus.utils import p_adjust_bh
from statannotations.Annotator import Annotator
from matplotlib import ticker as mticker
import numpy as np

fig, axs = plt.subplots(ncols = 3, nrows = 1, figsize = (12, 4))
order = data_f1.groupby('Tool')['F1'].median().sort_values(ascending = False).index
data_f1['Tool'] = pd.Categorical(data_f1['Tool'], categories = order)
sns.set_palette(sns.color_palette([tool_to_color[tool] for tool in order]))
plot(params = {'data': data_f1, 'x': 'Tool','y': 'F1'}, ax = axs[0])
axs[0].set_xticklabels(order, rotation=45, ha='right')
set_y_log(axs[0])
sns.despine(ax = axs[0])
pairs = []
for i in range(0, len(order) - 1, 1):
	pairs.append((order[i], order[i + 1]))
pairs = pairs[::-1]
p_values = [mannwhitneyu(data_f1.loc[data_f1['Tool'] == t1, 'F1'], data_f1.loc[data_f1['Tool'] == t2, 'F1']).pvalue for t1, t2 in pairs]
p_values = p_adjust_bh(p_values)
p_values_formatted = [f'p={p:.1e}' for p in p_values]
annotator = Annotator(axs[0], pairs, **{'data': data_f1, 'x': 'Tool','y': 'F1'})
annotator.set_custom_annotations(p_values_formatted)
annotator.annotate()
order = data_pr.groupby('Tool')['Precision'].median().sort_values(ascending = False).index
data_pr['Tool'] = pd.Categorical(data_pr['Tool'], categories = order)
sns.set_palette(sns.color_palette([tool_to_color[tool] for tool in order]))
plot(params = {'data': data_pr, 'x': 'Tool','y': 'Precision'}, ax = axs[1])
axs[1].set_xticklabels(order, rotation=45, ha='right')
set_y_log(axs[1])
sns.despine(ax = axs[1])
pairs = []
for i in range(0, len(order) - 1, 1):
	pairs.append((order[i], order[i + 1]))
pairs = pairs[::-1]
p_values = [mannwhitneyu(data_pr.loc[data_pr['Tool'] == t1, 'Precision'], data_pr.loc[data_pr['Tool'] == t2, 'Precision']).pvalue for t1, t2 in pairs]
p_values = p_adjust_bh(p_values)
p_values_formatted = [f'p={p:.1e}' for p in p_values]
annotator = Annotator(axs[1], pairs, **{'data': data_pr, 'x': 'Tool','y': 'Precision'})
annotator.set_custom_annotations(p_values_formatted)
annotator.annotate()
order = data_re.groupby('Tool')['Recall'].median().sort_values(ascending = False).index
data_re['Tool'] = pd.Categorical(data_re['Tool'], categories = order)
sns.set_palette(sns.color_palette([tool_to_color[tool] for tool in order]))
plot(params = {'data': data_re,'x': 'Tool','y': 'Recall'}, ax = axs[2])
axs[2].set_xticklabels(order, rotation=45, ha='right')
set_y_log(axs[2])
sns.despine(ax = axs[2])
pairs = []
for i in range(0, len(order) - 1, 1):
	pairs.append((order[i], order[i + 1]))
pairs = pairs[::-1]
p_values = [mannwhitneyu(data_re.loc[data_re['Tool'] == t1, 'Recall'], data_re.loc[data_re['Tool'] == t2, 'Recall']).pvalue for t1, t2 in pairs]
p_values = p_adjust_bh(p_values)
p_values_formatted = [f'p={p:.1e}' for p in p_values]
annotator = Annotator(axs[2], pairs, **{'data': data_re, 'x': 'Tool','y': 'Recall'})
annotator.set_custom_annotations(p_values_formatted)
annotator.annotate()
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'pr_re_f1.pdf'))
fig.savefig(os.path.join(out_dir, 'pr_re_f1.png'))

```

## Supplementary Fig S9, pannel j: Accuracy of target regions identified by SCENIC+, CellOracle and Pando

```python

tool_to_n_regions = {}
for tool in ['SCENIC+', 'Pando', 'CellOracle']:
	tool_to_n_regions[tool] = tool_to_regulon_metadata[tool].groupby('Consensus_name')['region'].apply(lambda x: len(set(x))).values

tool_to_color = {
          'SCENIC+': '#f70404',
          'GRaNIE': '#f9a204',
          'Pando': '#1c94fb',
          'CellOracle': '#228b21',
	  'SCENIC': '#ff1494',
	  'FigR': '#a020f0'}

from scipy.stats import mannwhitneyu

SCENIC = np.log10(list(tool_to_n_regions.values())[0])
Pando = np.log10(list(tool_to_n_regions.values())[1])
CellOracle = np.log10(list(tool_to_n_regions.values())[2])

data = pd.concat([
	pd.DataFrame({'SCENIC+': SCENIC}).melt().rename({'variable': 'tool', 'value': 'log(n_regions)'}, axis = 1),
	pd.DataFrame({'Pando': Pando}).melt().rename({'variable': 'tool', 'value': 'log(n_regions)'}, axis = 1),
	pd.DataFrame({'CellOracle': CellOracle}).melt().rename({'variable': 'tool', 'value': 'log(n_regions)'}, axis = 1)
])

from statannotations.Annotator import Annotator
pvalues = [
	mannwhitneyu(SCENIC, Pando).pvalue,
	mannwhitneyu(Pando, CellOracle).pvalue,
	mannwhitneyu(SCENIC, CellOracle).pvalue]
formatted_pvalues = [f'p={pvalue:.1e}' for pvalue in pvalues]
pairs = [('SCENIC+', 'Pando'), ('Pando', 'CellOracle'), ('SCENIC+', 'CellOracle')]

plotting_parameters = {
	'data': data,
	'x': 'tool',
	'y': 'log(n_regions)'}

from matplotlib import ticker as mticker
fig, ax = plt.subplots()
sns.violinplot(
	inner = None,
#	hue = 'tool',
	palette = tool_to_color,
	ax = ax,
	**plotting_parameters)
sns.stripplot(
	color="#2b2b2b", size = 4,
	ax =ax, alpha = 0.4,
	**plotting_parameters)
#plt.legend([],[], frameon=False)
annotator = Annotator(ax, pairs, **plotting_parameters)
annotator.set_custom_annotations(formatted_pvalues)
annotator.annotate()
ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
ymin, ymax = ax.get_ylim()
tick_range = np.arange(np.floor(ymin), ymax)
ax.yaxis.set_ticks(tick_range)
ax.yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)
ax.set_xticklabels(list(tool_to_n_regions.keys()))
ax.set_ylabel('Nr. of target regions per eRegulon')
fig.tight_layout()
sns.despine(ax = ax)
fig.savefig(os.path.join(outdir, 'n_target_regions_vln.pdf'))
fig.savefig(os.path.join(outdir, 'n_target_regions_vln.png'))
ax_n_targets = ax


data_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/extra_methods/pbmc/ChIP_seq_comparison/'
F1 = pd.read_csv(os.path.join(data_dir, 'F1.tsv'), sep = '\t')
PREC =  pd.read_csv(os.path.join(data_dir, 'PREC.tsv'), sep = '\t')
REC =  pd.read_csv(os.path.join(data_dir, 'REC.tsv'), sep = '\t')

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker as mticker
from adjustText import adjust_text
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
from statannotations.Annotator import Annotator
from itertools import combinations
from scipy.stats import mannwhitneyu
from scenicplus.utils import p_adjust_bh

def plot(ax, data, x, y, TFs_to_label):
	data = data.copy()
	tools_to_median = data.groupby(x)[y].median()
	tool_order = tools_to_median.sort_values(ascending = False).index
	data[y] = np.log10(data[y])
	data[x] = pd.Categorical(data[x], categories = tool_order)
	pairs = list(combinations(tool_order, 2))
	p_values = [
		mannwhitneyu(data.loc[data[x] == tool_a, y], data.loc[data[x] == tool_b, y]).pvalue
		for tool_a, tool_b in pairs]
	p_values = p_adjust_bh(p_values)
	#pairs = [pair for pair, pval in zip(pairs, p_values) if pval < 0.05]
	#p_values = [pval for pval in p_values if pval < 0.05]
	p_values_formatted = [f'p={p:.1e}' for p in p_values]
	sns.set_palette(sns.color_palette([tool_to_color[tool] for tool in tool_order]))
	sns.violinplot(
		inner = None, ax = ax,
		data = data, x = x, y = y)
	sns.swarmplot(
		color="#2b2b2b", size = 4, ax = ax,
		data = data, x = x, y = y)
	y_TF_to_label = {TF: 
		[
			data.loc[np.logical_and(data[x] == tool, data['TF'] == TF), y].values[0]
			if TF in data.loc[data[x] == tool, 'TF'].to_numpy() else None
			for tool in tool_order]
		for TF in TFs_to_label}
	labels = []
	xs = []
	ys = []
	for tool_id in range(len(tool_order)):
		for TF in y_TF_to_label.keys():
			if not y_TF_to_label[TF][tool_id] == None:
				#print(TF)
				xy = ax.collections[tool_id + len(tool_order)].get_offsets()
				xy_TF = xy[
					np.where(
				[np.round(_xy[1], 5) == np.round(y_TF_to_label[TF][tool_id], 5) for _xy in xy])[0][0]]
				labels.append(ax.text(xy_TF[0], xy_TF[1], TF, fontweight = 'bold', fontname = "Arial", fontsize = 11))
				xs.append(xy_TF[0])
				ys.append(xy_TF[1])
	adjust_text(labels, 
		ax = ax, 
		add_objects = ax.collections[0:3], 
		x = xs, y = ys, 
		arrowprops=dict(arrowstyle='->', lw = 2, color = 'black'),
		autoalign = 'x')
		#expand_points = (2.5, 1.5))
	#ax.set_yscale('log')
	ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
	ymin, ymax = ax.get_ylim()
	tick_range = np.arange(np.floor(ymin), ymax)
	ax.yaxis.set_ticks(tick_range)
	ax.yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)
	#ax.set_ylim(-7, 0)
	annotator = Annotator(ax, pairs, data = data, x = x, y = y)
	annotator.set_custom_annotations(p_values_formatted)
	annotator.annotate()
	sns.despine(ax = ax)

TFs_to_label = ['PAX5', 'CEBPB', 'BCL11A', 'EBF1']
out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/figures_reviewer_questions/pbmc_other_eGRN_methods'

fig, axs = plt.subplots(ncols = 3, nrows = 1, figsize = (15, 5))
plot(axs[0], F1, 'tool', 'F1', TFs_to_label)
plot(axs[1], PREC, 'tool', 'Precision', TFs_to_label)
plot(axs[2], REC, 'tool', 'Recall', TFs_to_label)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'f1_prec_rec.pdf'))
fig.savefig(os.path.join(out_dir, 'f1_prec_rec.png'))

```



