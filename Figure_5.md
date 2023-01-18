# Code for generating the pannels of figure 5 and accompanying supplementary figures

## Setup python environment and load data

```python

import dill
import os
import matplotlib
from scenicplus.utils import *

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

f_scplus_obj = '/staging/leuven/stg_00002/lcb/saibar/Projects/PanCancer/1_runs_separateCancers/MMlines_1_sc/scenicplus_output_v10_mmRegions/scplus_obj.pkl'
plot_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2/plots'
out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2'
work_dir = '/staging/leuven/stg_00002/lcb/saibar/Projects/PanCancer/1_runs_separateCancers/MMlines_1_sc/scenicplus_output_v10_mmRegions'

scplus_obj = dill.load(open(f_scplus_obj, 'rb'))

region_ranking = dill.load(open(os.path.join(work_dir, 'region_ranking.pkl'), 'rb'))
gene_ranking = dill.load(open(os.path.join(work_dir, 'gene_ranking.pkl'), 'rb'))

## Filter eRegulons
from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons

apply_std_filtering_to_eRegulons(scplus_obj)

```

## Supplementary Figure S21, pannel a: eRegulon selection

```python

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
        variable = 'MMline',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Gene_based'
)
generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'MMline',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Region_based'
)

TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'MMline',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Gene_based',
            out_key = 'filtered_gene_based')
TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'MMline',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Region_based',
            out_key = 'filtered_region_based')

n_targets = [int(x.split('(')[1].replace('r)', '')) for x in scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Cistrome']]
rho = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'].to_list()
adj_pval = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Adjusted_p-value'].to_list()
thresholds = {
        'rho': [-0.75, 0.65],
        'n_targets': 0
}
import seaborn as sns
fig, ax = plt.subplots(figsize = (5.5, 4.5))
sc = ax.scatter(rho, n_targets, c = -np.log10(adj_pval), s = 5)
ax.set_xlabel('Correlation coefficient')
ax.set_ylabel('nr. target regions')
#ax.hlines(y = thresholds['n_targets'], xmin = min(rho), xmax = max(rho), color = 'black', ls = 'dashed', lw = 1)
ax.vlines(x = thresholds['rho'], ymin = 0, ymax = max(n_targets), color = 'black', ls = 'dashed', lw = 1)
ax.text(x = thresholds['rho'][0], y = max(n_targets), s = str(thresholds['rho'][0]))
ax.text(x = thresholds['rho'][1], y = max(n_targets), s = str(thresholds['rho'][1]))
sns.despine(ax = ax)
fig.colorbar(sc, label = '-log10(adjusted_pvalue)', ax = ax)
fig.savefig(os.path.join(plot_dir, 'eRegulon_selection_n_targets_rho_stringent.pdf'))

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
scplus_obj.uns['selected_eRegulons'] = {'Gene_based': selected_eRegulons_gene_sig, 'Region_based': selected_eRegulons_region_sig}

```

## Export Loom file for visualization in R

```python

from scenicplus.dimensionality_reduction import run_eRegulons_pca
run_eRegulons_pca(
	scplus_obj,
	auc_key = 'eRegulon_AUC_filtered',
	reduction_name = 'eRegulons_PCA',
	selected_regulons = flatten_list([scplus_obj.uns['selected_eRegulons'][k] for k in scplus_obj.uns['selected_eRegulons'].keys()]))

## Save changes
dill.dump(scplus_obj, open(os.path.join(out_dir, 'scplus_obj.pkl'), 'wb'))

binarize_AUC(scplus_obj,
             auc_key='eRegulon_AUC_filtered',
             out_key='eRegulon_AUC_filtered_thresholds',
             signature_keys=['Gene_based'],
             n_cpu=1)

def remove_second_sign(x):
        if 'extended' not in x:
                TF, first, second, n = x.split('_')
                return f'{TF}_{first}_{n}'
        else:
                TF, extended, first, second, n = x.split('_')
                return f'{TF}_{extended}_{first}_{n}'

scplus_obj.uns['eRegulon_metadata_filtered']['Gene_signature_name'] = [remove_second_sign(x) for x in scplus_obj.uns['eRegulon_metadata_filtered']['Gene_signature_name']]
scplus_obj.uns['eRegulon_metadata_filtered']['Region_signature_name'] = [remove_second_sign(x) for x in scplus_obj.uns['eRegulon_metadata_filtered']['Region_signature_name']]


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
       out_fname=os.path.join(out_dir,'SCENIC+_gene_based.loom'))

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
       out_fname=os.path.join(out_dir,'SCENIC+_region_based.loom'))

with open(os.path.join(out_dir, 'selected_eRegulons.txt'), 'w') as f:
        for eRegulon_name in scplus_obj.uns['selected_eRegulons']['Gene_based']:
                f.write(eRegulon_name.rsplit('_', 1)[0])
                f.write('\n')

```

## Main Figure 5, pannel A: PCA plot

```r

plot_dir <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2/plots'
out_dir <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2'

# load data
library(SCopeLoomR)
loom <- open_loom(paste0(out_dir, '/SCENIC+_gene_based.loom'))
cell_data <- get_cell_annotation(loom)
embeddings <- get_embeddings(loom)
eRegulons_PCA <- embeddings$eRegulons_PCA
colnames(eRegulons_PCA) <- c('PCA_1', 'PCA_2')
cell_data[] <- lapply(cell_data, sub, pattern = " ", replacement = "-")
cell_data[] <- lapply(cell_data, sub, pattern = "-(.*)", replacement = "")

cell_plot_data <- cbind(eRegulons_PCA, cell_data[rownames(eRegulons_PCA),])

color_dict = c(
    'MM047' = '#a6747b',
    'MM029' = '#a64b57',
    'MM099' = '#a63241',
    'MM001' = '#d6c496',
    'MM011' = '#d6b560',
    'MM031' = '#d6ac40',
    'MM074' = '#8c7862',
    'MM087' = '#8c673f',
    'MM057' = '#8c5d2a'
)


library(ggplot2)
source('/staging/leuven/stg_00002/lcb/cbravo/software/plotting_aux.R')
plot <- ggplot(cell_plot_data, aes(x=PCA_1, y=PCA_2, colour=MMline))
plot <- plot + geom_point(size = 0.4)
plot <- plot + theme_classic()
plot <- plot + theme(legend.position = "none")
plot <- plot + labs(x = NULL, y = NULL)
plot <- plot + scale_color_manual(values = color_dict)
plot <- plot + guides(x = "none", y = "none")
ggsave(filename = paste0(plot_dir, '/R_eRegulons_PCA.png'), device='png', bg = "transparent", width=5.25, height=4.25)

pdf(paste0(plot_dir, '/R_eRegulons_PCA.pdf'), width = 5.25, height = 4.25)
LabelClusters(plot, 'MMline', split.by ='MMline', box=FALSE, repel=TRUE)  + scale_color_manual(values = color_dict)
dev.off()

```


## Main Figure 5, pannel b: heatmap-dotplot

```r

plot_dir <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2/plots'
work_dir <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2'

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
rss_values <- calcRSS(cells_AUC_gene, cell_data$MMline)
rss_values <- sweep(rss_values,2,colSums(rss_values),`/`)*100
rss_values <- rss_values[,sort(colnames(rss_values))]

# Prepare data for plotting

expression_list <- list()
for (x in unique(cell_data$MMline)){
  print(x)
  expression_list[[x]] <- as.data.frame(t(log(rowSums(dgem[,grep(x, cell_data$MMline, fixed = TRUE)])/sum(rowSums(dgem[,grep(x, cell_data$MMline, fixed = TRUE)]))*10^6+1)))
}
exp_mat <- t(data.table::rbindlist(expression_list))
colnames(exp_mat) <- names(expression_list)
exp_mat <- exp_mat

line_order <- c('MM001', 'MM011', 'MM031', 'MM074', 'MM087', 'MM057', 'MM047', 'MM029', 'MM099')

sel_rel <- c(positive, negative)
order_list <- list()
for (i in 1:ncol(rss_values)){
  order_list[[i]] <- vector()
}
for (x in sel_rel){
  i <- which.max(rss_values[x,line_order])
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
sel_rel_color[which(sel_rel %in% positive)] <- 'forestgreen'
sel_rel_color[which(sel_rel %in% negative)] <- 'red'

#plot
library(ggplot2)
rel_data$Regulon <- factor(rel_data$Regulon, levels=sel_rel)
rel_data$Cell_type <- gsub('_', ' ', rel_data$Cell_type)
rel_data$Cell_type <- factor(rel_data$Cell_type, levels= line_order)
rel_data$RSS <- as.numeric(rel_data$RSS)
rel_data$Expression <- as.numeric(rel_data$Expression)
g <- ggplot(data = rel_data, mapping = aes_string(x = 'Cell_type', y = 'Regulon')) +
    geom_point(mapping = aes_string(size = 'RSS', fill = 'Expression'), colour="black",pch=21) +
    scale_radius(range = c(0.1, 3), limits=c(0, 1)) +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip() + theme(axis.text.x = element_text(colour = sel_rel_color)) +
    theme(legend.position = "none")
ggsave(paste0(plot_dir, '/R_Dotplot_RSS_no_legend.pdf'), g, width=11.5, height=4.5)


library(ggplot2)
rel_data$Regulon <- factor(rel_data$Regulon, levels=sel_rel)
rel_data$Cell_type <- gsub('_', ' ', rel_data$Cell_type)
rel_data$Cell_type <- factor(rel_data$Cell_type, levels= line_order)
rel_data$RSS <- as.numeric(rel_data$RSS)
rel_data$Expression <- as.numeric(rel_data$Expression)
g <- ggplot(data = rel_data, mapping = aes_string(x = 'Cell_type', y = 'Regulon')) +
    geom_point(mapping = aes_string(size = 'RSS', fill = 'Expression'), colour="black",pch=21) +
    scale_radius(range = c(0.1, 3), limits=c(0, 1)) +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip() + theme(axis.text.x = element_text(colour = sel_rel_color))
ggsave(paste0(plot_dir, '/R_Dotplot_RSS.pdf'), g, width=11.5, height=4.5)

rel_data$Regulon <- factor(rel_data$Regulon, levels=rev(sel_rel))
rel_data$Cell_type <- gsub('_', ' ', rel_data$Cell_type)
rel_data$Cell_type <- factor(rel_data$Cell_type, levels= line_order)
rel_data$RSS <- as.numeric(rel_data$RSS)
rel_data$Expression <- as.numeric(rel_data$Expression)
rel_data$Expression_scale <- ave(rel_data$Expression, rel_data$Regulon, FUN = scale)
rel_data$Repressor_activator <- unlist(lapply(grepl("+", rel_data$Regulon, fixed=TRUE), function(x) if (x) 'Activators' else 'Repressors'))
g <- ggplot(data = rel_data, mapping = aes_string(x = 'Cell_type', y = 'Regulon')) +
    facet_grid(cols = vars(Repressor_activator), scales = "free_x", space = 'free') +
    geom_tile(mapping = aes_string(fill = 'Expression_scale')) +
    geom_point(mapping = aes_string(size = 'RSS'), colour="black",pch=21, fill = 'black') +
    scale_radius(range = c(0.1, 3), limits=c(0, 1)) +
    scale_fill_distiller(palette = "RdYlBu") +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip() +
    theme(legend.position = "none")
ggsave(paste0(plot_dir, '/R_Dotplot_RSS_hm_no_legend.pdf'), g, width=11.5, height=4.5)

rel_data$Regulon <- factor(rel_data$Regulon, levels=rev(sel_rel))
rel_data$Cell_type <- gsub('_', ' ', rel_data$Cell_type)
rel_data$Cell_type <- factor(rel_data$Cell_type, levels= line_order)
rel_data$RSS <- as.numeric(rel_data$RSS)
rel_data$Expression <- as.numeric(rel_data$Expression)
rel_data$Expression_scale <- ave(rel_data$Expression, rel_data$Regulon, FUN = scale)
rel_data$Repressor_activator <- unlist(lapply(grepl("+", rel_data$Regulon, fixed=TRUE), function(x) if (x) 'Activators' else 'Repressors'))
g <- ggplot(data = rel_data, mapping = aes_string(x = 'Cell_type', y = 'Regulon')) +
    facet_grid(cols = vars(Repressor_activator), scales = "free_x", space = 'free') +
    geom_tile(mapping = aes_string(fill = 'Expression_scale')) +
    geom_point(mapping = aes_string(size = 'RSS'), colour="black",pch=21, fill = 'black') +
    scale_radius(range = c(0.1, 3), limits=c(0, 1)) +
    scale_fill_distiller(palette = "RdYlBu") +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip()
ggsave(paste0(plot_dir, '/R_Dotplot_RSS_hm.pdf'), g, width=11.5, height=4.5)

rel_data$Regulon <- factor(rel_data$Regulon, levels=rev(sel_rel))
rel_data$Cell_type <- gsub('_', ' ', rel_data$Cell_type)
rel_data$Cell_type <- factor(rel_data$Cell_type, levels= line_order)
rel_data$RSS <- as.numeric(rel_data$RSS)
rel_data$Expression <- as.numeric(rel_data$Expression)
rel_data$Expression_scale <- ave(rel_data$Expression, rel_data$Regulon, FUN = scale)
rel_data$Repressor_activator <- unlist(lapply(grepl("+", rel_data$Regulon, fixed=TRUE), function(x) if (x) 'Activators' else 'Repressors'))
g <- ggplot(data = rel_data[rel_data$Repressor_activator == 'Activators', ], mapping = aes_string(x = 'Cell_type', y = 'Regulon')) +
    geom_tile(mapping = aes_string(fill = 'Expression_scale')) +
    geom_point(mapping = aes_string(size = 'RSS'), colour="black",pch=21, fill = 'black') +
    scale_radius(range = c(0.1, 3), limits=c(0, 1)) +
    scale_fill_distiller(palette = "RdYlBu") +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip() +
    theme(legend.position = "none")
ggsave(paste0(plot_dir, '/R_Dotplot_RSS_hm_no_legend_activators.pdf'), g, width=11.5, height=4.5)

```

## Supplementary Figure S21, pannel d: eRegulon AUC values on PCA

```python

import dill
import os
import matplotlib
from scenicplus.utils import *
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

plot_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2/plots'
out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2'

from scenicplus.dimensionality_reduction import plot_AUC_given_ax

scplus_obj = dill.load(open(os.path.join(out_dir, 'scplus_obj.pkl'), 'rb'))

eRegulons_to_plot = ['MITF_+_(993g)', 'SOX10_+_(863g)', 'JUN_+_(941g)', 'ZEB1_+_(572g)']
fig, axs = plt.subplots(nrows = 2, ncols = 2, figsize = (5.5, 4.5))
for eRegulon, ax in zip(eRegulons_to_plot, axs.ravel()):
	plot_AUC_given_ax(
		scplus_obj = scplus_obj,
		reduction_name = 'eRegulons_PCA',
		feature = eRegulon,
		ax = ax,
		auc_key = 'eRegulon_AUC_filtered',
		dot_size = 0.4)
	sns.despine(ax = ax, top = True, bottom = True, left = True, right = True)
	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_xlabel(None)
	ax.set_ylabel(None)
	ax.set_title(eRegulon, fontsize = 6)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'eRegulons_AUC_on_PCA.pdf'))

```

## Supplementary Figure S21, pannel e: overlap of target regions of TFAP2A, RUNX3, MITF and SOX10

```r

plot_dir <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2/plots'
work_dir <- '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2'

# load data
library(SCopeLoomR)
# Get eRegulons
loom <- open_loom(paste0(work_dir, '/SCENIC+_region_based.loom'))
regulons_region <- get_regulons(loom, column.attr.name='Regulons')
regulon_region_list <- list()
for (row in rownames(regulons_region)){
  regulon_region_list[[row]] <- colnames(regulons_region)[which(regulons_region[row,] == 1)]
}
names(regulon_region_list) <- gsub('_extended', '', names(regulon_region_list))

regulon_region_list <- tapply(unlist(regulon_region_list , use.names = FALSE), rep(names(regulon_region_list), lengths(regulon_region_list)), FUN = c)

regulons_of_interest <- c('MITF_+', 'TFAP2A_+', 'RUNX3_+', 'SOX10_+')

library(ggplot2)
g <- ggVennDiagram(regulon_region_list[regulons_of_interest], label = 'count')
ggsave(paste0(plot_dir, '/R_VENN_MTRS.pdf'), g, width=4.5, height=4.5)

library(UpSetR)

g <- upset(fromList(regulon_region_list[regulons_of_interest]), order.by = 'freq')

pdf(paste0(plot_dir, '/R_upset_MTRS.pdf'), width=4.5, height=4.5)
g
dev.off()

```

## Supplementary Figure S21, pannel f: ChIP-seq coverages

```python

import dill
import os
import matplotlib
from scenicplus.utils import *
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

plot_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2/plots'
out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2'

from scenicplus.dimensionality_reduction import plot_AUC_given_ax

scplus_obj = dill.load(open(os.path.join(out_dir, 'scplus_obj.pkl'), 'rb'))


bw_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/melanoma_mixed_lines/hg38/ChIP-seq'

f_MITF_bw = os.path.join(bw_dir, 'MITF_ChIP.bw')
f_TFAP2A_bw = os.path.join(bw_dir, 'TFAP2A_ChIP.bw')
f_SOX10_bw = os.path.join(bw_dir, 'SOX10_ChIP.bw')

import pyBigWig

MITF_bw = pyBigWig.open(f_MITF_bw)
TFAP2A_bw = pyBigWig.open(f_TFAP2A_bw)
SOX10_bw = pyBigWig.open(f_SOX10_bw)

def region_name_to_coord(r):
        return r.replace('-', ':').split(':')

from tqdm import tqdm
def calc_coverage(regions, n_bins, bw, offset = 500):
        coverage = []
        for region in tqdm(regions, total = len(regions)):
                chrom, start, end = region_name_to_coord(region)
                coverage.append(bw.stats(chrom, int(start) - offset, int(end) + offset, nBins = n_bins))
        return pd.DataFrame(coverage)

n_bins = 50
region_sets_of_interest = {
	'SOX10': set(scplus_obj.uns['eRegulon_signatures_filtered']['Region_based']['SOX10_+_(3192r)']), 
	'TFAP2A':  set(scplus_obj.uns['eRegulon_signatures_filtered']['Region_based']['TFAP2A_+_(758r)']),
	'MITF': set(scplus_obj.uns['eRegulon_signatures_filtered']['Region_based']['MITF_+_(1252r)'])}
bw_of_interest = {
	'SOX10': SOX10_bw,
	'MITF': MITF_bw,
	'TFAP2A': TFAP2A_bw}
eRegulons_of_interest = set(['SOX10', 'MITF', 'TFAP2A'])

from itertools import combinations

coverages = {}
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
                                coverages[f'c_{ChIP_seq_track}_r_{"_n_".join(combination)}_d_{"_".join(others)}'] = calc_coverage(regions, n_bins, bw_of_interest[ChIP_seq_track])

color_dict = {
	'SOX10': '#A4303F',
	'MITF': '#477998',
	'TFAP2A': '#6A994E'
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

fig, axs = plt.subplots(ncols = 3, nrows = 3, figsize = (3.25, 4.5))
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
				axs[i, j].set_yticks([0,1])
				axs[i,j].set_xticks(np.arange(0, 51, 16.66666))
				axs[i, j].set_xticklabels(np.round((np.arange(0, 51, 16.66666) * 30 - 249.5)).astype(int), fontsize = 5)
		else:
			sns.despine(ax = axs[i, j], top = True, bottom = True, left = True, right = True)
			axs[i, j].set_xticks([])
			axs[i, j].set_yticks([])
handles, labels = axs[0,0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center right', frameon = False)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'aggreation_chip_combinations_TFAP2A_MITF_SOX10.pdf'))

fig, ax = plt.subplots()
for data_name in [x for x in coverages.keys() if x.endswith('_d_')]:
	ChIP_seq_track = data_name.split('_r_')[0].replace('c_', '')
	ChIP_min, ChIP_max = min_max_factor[ChIP_seq_track]
	coverage = coverages[data_name]
	ax.plot(
		np.arange(n_bins),
		(coverage.mean(0) - ChIP_min) / (ChIP_max - ChIP_min),
		lw = 1,
		label = ChIP_seq_track,
		color = color_dict[ChIP_seq_track])
	sns.despine(ax = ax)
	ax.set_xticks(np.arange(0, 51, 16.66666))
	ax.set_xticklabels(np.round((np.arange(0, 51, 16.66666) * 30 - 249.5)).astype(int), fontsize = 5)
fig.legend(frameon = False)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'aggreation_chip_TFAP2A_MITF_SOX10.pdf'))

```

## Supplementary Figure S21, pannel g: SOX10 KD vulcano plot

```python

import dill
import os
import matplotlib
from scenicplus.utils import *
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

plot_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2/plots'
out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2'

scplus_obj = dill.load(open(os.path.join(out_dir, 'scplus_obj.pkl'), 'rb'))

RNA_KD_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/melanoma_mixed_lines/hg38/KD_data/v2/analysis/RNA'

from random import sample
from scipy import stats
from scipy.stats import ttest_ind
TF_KD = 'SOX10'
TF_logfc = pd.read_csv(os.path.join(RNA_KD_dir, f'{TF_KD}_v_NTC_log2fc.tsv'), sep = '\t')
TF_logfc = TF_logfc.dropna()
logfcs = {}
gene_set_lens = []
for TF in set(scplus_obj.uns['eRegulon_metadata_filtered']['TF']):
    genset_names_pos = [k for k in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys() if TF in k and '+' in k and 'extended' not in k]
    genset_names_neg = [k for k in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys() if TF in k and '-' in k and 'extended' not in k]
    if len(genset_names_pos) > 0:
        genset_name_pos = genset_names_pos[0]
    else:
        genset_name_pos = None

    if len(genset_names_neg) > 0:
        genset_name_neg = genset_names_neg[0]
    else:
        genset_name_neg = None

    #genset_name_neg = [k for k in scplus_obj.uns['eRegulon_signatures']['Gene_based'].keys() if TF in k and '-' in k and 'extended' not in k][0]
    gene_set_pos = scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'][genset_name_pos] if genset_name_pos is not None else None
    gene_set_neg = scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'][genset_name_neg] if genset_name_neg is not None else None
    gene_set_pos = list(set(gene_set_pos) & set(TF_logfc.index)) if gene_set_pos is not None else None
    gene_set_neg = list(set(gene_set_neg) & set(TF_logfc.index)) if gene_set_neg is not None else None
    if genset_name_pos is not None:
        logfcs[genset_name_pos] = TF_logfc.loc[gene_set_pos]['stat'].to_numpy()
        gene_set_lens.append(len(list(set(gene_set_pos) & set(TF_logfc.index))))
    if genset_name_neg is not None:
        logfcs[genset_name_neg] = TF_logfc.loc[gene_set_neg]['stat'].to_numpy()
        gene_set_lens.append(len(list(set(gene_set_neg) & set(TF_logfc.index))))
random_gene_set = sample(list(TF_logfc.index), k = round(np.mean(gene_set_lens)))
pvals = {k: ttest_ind(logfcs[k], TF_logfc.loc[random_gene_set]['stat'].to_numpy()).pvalue for k in logfcs.keys()}
TF_change = {k: TF_logfc.loc[k.split('_')[0]]['log2FoldChange'] for k in logfcs.keys() if k.split('_')[0] in TF_logfc.index}
TF_to_use = set(pvals.keys()) & set(TF_change.keys())
TF_to_use = [TF for TF in TF_to_use if '+' in TF]

scale = lambda X: [(x - min(X)) / (max(X) - min(X)) for x in X]
from adjustText import adjust_text
fig, ax = plt.subplots(figsize = (3.25, 4.5))
ax.scatter(
        [np.mean(logfcs[TF]) for TF in TF_to_use],
        [-np.log10(pvals[TF]) for TF in TF_to_use],
        c = stats.zscore([TF_change[TF] for TF in TF_to_use]),
        cmap = 'bwr',
        vmin = -0.5,
        vmax = 0.5,
	s = 3
)
texts = []
for x, y, s in zip(
    [np.mean(logfcs[TF]) for TF in TF_to_use],
    [-np.log10(pvals[TF]) for TF in TF_to_use],
    TF_to_use):
    if  (y > 4 and x < 0) or (y > 3 and x > 0):
        texts.append(ax.text(x, y, s, fontsize = 5))
ax.set_xlabel('avg. logFoldChange of target genes')
ax.set_ylabel('-log10Pval')
ax.set_xlim((-1, 1.2))
adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
sns.despine(ax = ax)
fig.savefig(os.path.join(plot_dir,'vulcano_sox10_kd.pdf'))

```

## Main Figure 5, pannel d: Perturbation simulations over iterations


```python

import dill
import os
import matplotlib
from scenicplus.utils import *
import seaborn as sns
from scenicplus.simulation import *

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

plot_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2/plots'
out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2'

simulation_out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/melanoma_mixed_lines/hg38_v10_clust_consensus/outs'
scplus_obj = dill.load(open(os.path.join(simulation_out_dir, 'scplus_obj_w_perturbation.pkl'), 'rb'))

# genes over iter
flatten_list = lambda t: [item for sublist in t for item in sublist]
DEGs = list(set(flatten_list([list(scplus_obj.uns['DEGs']['lineState'][k].index) for k in scplus_obj.uns['DEGs']['lineState'].keys()])))
genes_to_use = list(set([*DEGs, * [x.split('_')[0] for x in scplus_obj.uns['selected_eRegulons']['Gene_based']]]))

regressors = train_gene_expression_models(
        scplus_obj,
        eRegulon_metadata_key = 'eRegulon_metadata',
        genes = genes_to_use,
        eRegulons_to_use = None)

perturbed_matrices = simulate_perturbation(
        scplus_obj,
        perturbation = {'SOX10': 0},
        regressors = regressors,
        keep_intermediate = True,
        n_iter = 10)

mtx_over_time = perturbed_matrices

factors_of_interest = [
        'TFAP2A', 'MITF', 'IRF4',
        'JUN', 'FOSL1', 'FOSL2', 'NFE2L3'
]

color_dict = {
        'TFAP2A': '#CC7E85',
        'MITF': '#FCAA67',
        'IRF4': '#F03A47',
        'JUN': '#26FFE6',
        'FOSL1': '#037971',
        'FOSL2': '#A9FFF7',
        'NFE2L3': '#00B295',
}

fig, ax = plt.subplots(figsize = (3.5, 4.0))
for factor in factors_of_interest:
        expr = [
        np.log1p(
                mtx_over_time[str(i)].loc[['MM001' in bc for bc in mtx_over_time[str(i)].index]].T/(mtx_over_time[str(i)].loc[['MM001' in bc for bc in mtx_over_time[str(i)].index]].sum(1) * 1e6)).mean(1)[factor] for i in range(12)]
        fc = [(expr[i] + 0.0000000001) / (expr[0] + 0.0000000001) for i in range(12)]
        ax.plot(np.log2(fc), label = factor, color = color_dict[factor])
ax.legend(loc = 'lower center', ncol = 3, fancybox = False, shadow = False, frameon = False)
ax.set_xlabel('iteration nr.')
ax.set_ylabel('Predicted log2FC')
ax.set_ylim((-0.6, 0.7))
ax.hlines(y = 0, xmin = 0, xmax = 12, lw = 1)
sns.despine(ax = ax)
ax.set_yticks([-0.6, -0.3, 0, 0.3, 0.6])
ax.set_xticks([0, 4, 8, 12])
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'perturb_iter.pdf'))

```

## Main Figure 5, pannel f: Perturbation arrows and pannel g: heatmap

```python

import dill
import os
import matplotlib
from scenicplus.utils import *
import seaborn as sns
from scenicplus.simulation import *

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

plot_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2/plots'
out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2'

scplus_obj = dill.load(open(os.path.join(out_dir, 'scplus_obj.pkl'), 'rb'))

from pycisTopic.diff_features import *
from pycisTopic.signature_enrichment import *
def make_rankings(X, seed = 123):
        rng = np.random.default_rng(seed=seed)
        # Function to make rankings per array
        def rank_scores_and_assign_random_ranking_in_range_for_ties(
                scores_with_ties_for_motif_or_track_numpy: np.ndarray
        ) -> np.ndarray:
                #
                # Create random permutation so tied scores will have a different ranking each time.
                random_permutations_to_break_ties_numpy = rng.permutation(
                scores_with_ties_for_motif_or_track_numpy.shape[0]
                )
                ranking_with_broken_ties_for_motif_or_track_numpy = random_permutations_to_break_ties_numpy[
                (-scores_with_ties_for_motif_or_track_numpy)[
                    random_permutations_to_break_ties_numpy].argsort()
                ].argsort().astype(imputed_acc_obj_ranking_db_dtype)

                return ranking_with_broken_ties_for_motif_or_track_numpy
        imputed_acc_ranking = CistopicImputedFeatures(
                np.zeros((len(X.columns), len(X.index))),
                            X.columns,
                            X.index,
                            'Ranking')
        imputed_acc_obj_ranking_db_dtype = 'uint32'
        mtx = X.T.to_numpy()
        for col_idx in range(len(imputed_acc_ranking.cell_names)):
                imputed_acc_ranking.mtx[:, col_idx] = rank_scores_and_assign_random_ranking_in_range_for_ties(
                        mtx[:, col_idx].toarray().flatten() if sparse.issparse(mtx) else mtx[:, col_idx].flatten())
        return imputed_acc_ranking

from tqdm import tqdm
TFs_of_interest = list(set([x.split('_')[0] for x in scplus_obj.uns['selected_eRegulons']['Gene_based']]))

import logging
logging.basicConfig(level=logging.CRITICAL)
import warnings
from scenicplus.eregulon_enrichment import score_eRegulons
from scenicplus.dimensionality_reduction import run_eRegulons_pca
for TF in tqdm(TFs_of_interest, total = len(TFs_of_interest)):
        perturbed_matrices = simulate_perturbation(
                scplus_obj,
                perturbation = {TF: 0},
                regressors = regressors,
                keep_intermediate = True,
                n_iter = 3)
        perturbed_rankings = {k: make_rankings(perturbed_matrices[k]) for k in perturbed_matrices.keys()}
        for k in perturbed_rankings.keys():
                with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        score_eRegulons(
                                scplus_obj,
                                ranking = perturbed_rankings[k],
                                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                                key_added = f'{TF}_KD_sim_eRegulon_AUC_iter_{k}',
                                enrichment_type = 'gene',
                                n_cpu = 5)
                run_eRegulons_pca(
                        scplus_obj,
                        auc_key = f'{TF}_KD_sim_eRegulon_AUC_iter_{k}',
                        reduction_name = f'{TF}_KD_sim_eRegulons_PCA_iter_{k}',
                        selected_regulons = scplus_obj.uns['selected_eRegulons']['Gene_based'],
                        signature_keys = ['Gene_based'])
dill.dump(scplus_obj, open(os.path.join(out_dir, 'scplus_obj_w_perturbation_no_auto_reg.pkl'), 'wb'))

shifts_PC0 = {}
shifts_PC1 = {}

for TF in TFs_of_interest:
        delta_embedding = _project_perturbation_in_embedding(
                scplus_obj,
                original_matrix = scplus_obj.uns[f'{TF}_KD_sim_eRegulon_AUC_iter_4']['Gene_based'],
                perturbed_matrix = scplus_obj.uns[f'{TF}_KD_sim_eRegulon_AUC_iter_0']['Gene_based'],
                reduction_name = f'{TF}_KD_sim_eRegulons_PCA_iter_0')
        mean_shift = pd.DataFrame(delta_embedding).groupby(scplus_obj.metadata_cell['MMline'].to_numpy()).mean()
        shifts_PC0[TF] = mean_shift[0]
        shifts_PC1[TF] = mean_shift[1]

line_color_dict = {
    'MM047': '#a6747b',
    'MM029': '#a64b57',
    'MM099': '#a63241',
    'MM001': '#d6c496',
    'MM011': '#d6b560',
    'MM031': '#d6ac40',
    'MM074': '#8c7862',
    'MM087': '#8c673f',
    'MM057': '#8c5d2a'
}

shift_df = pd.DataFrame(shifts_PC0).T
factors_to_plot = [
        *shift_df.max(1).sort_values(ascending = False).head(10).index,
        *reversed(shift_df.min(1).sort_values(ascending = True).head(10).index)]
line_order = ['MM001', 'MM011', 'MM031', 'MM087', 'MM074', 'MM057', 'MM047', 'MM029', 'MM099']
import seaborn as sns
fig, ax = plt.subplots(figsize = (7, 4.25))
sns.heatmap(
        shift_df.loc[factors_to_plot, line_order].T,
        yticklabels=True,vmin = -0.4, vmax = 0.4, ax = ax, cmap = 'bwr')
for ytick in ax.get_yticklabels():
        ytick.set_color(line_color_dict[ytick.get_text()])
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'hm_PC0_shift.pdf'))

scale = lambda X: [(x - min(X)) / (max(X) - min(X)) for x in X]
n_grid_cols = 25

TF = 'SOX10'
fig, ax = plt.subplots(figsize = (5.24, 4.25))
embedding = scplus_obj.dr_cell[f'{TF}_KD_sim_eRegulons_PCA_iter_0'].to_numpy()
delta_embedding = _project_perturbation_in_embedding(
	scplus_obj,
	original_matrix = scplus_obj.uns[f'{TF}_KD_sim_eRegulon_AUC_iter_4']['Gene_based'],
	perturbed_matrix = scplus_obj.uns[f'{TF}_KD_sim_eRegulon_AUC_iter_0']['Gene_based'],
	reduction_name = f'{TF}_KD_sim_eRegulons_PCA_iter_0')
grid_xy, uv, mask = _calculate_grid_arrows(embedding, delta_embedding, 0.005, n_grid_cols, n_grid_cols, 25, 5)
plot_metadata_given_ax(
	scplus_obj=scplus_obj,
	reduction_name=f'{TF}_KD_sim_eRegulons_PCA_iter_0',
	ax = ax,
	variable = 'MMline',
	color_dictionary = {'MMline':line_color_dict},
	show_label = False,
	show_legend = False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel(None)
ax.set_ylabel(None)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_title(TF)
distances = np.sqrt((uv**2).sum(1))
norm = matplotlib.colors.Normalize(vmin=0.15, vmax=0.5, clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.Greys)
uv[np.logical_or(~mask, np.array(scale(distances)) < 0.15)] = np.nan
ax.streamplot(
	grid_xy.reshape(n_grid_cols,n_grid_cols, 2)[:, :, 0],
	grid_xy.reshape(n_grid_cols,n_grid_cols, 2)[:, :, 1],
	uv.reshape(n_grid_cols,n_grid_cols, 2)[:, :, 0],
	uv.reshape(n_grid_cols,n_grid_cols, 2)[:, :, 1], density = 3, color = np.array(scale(distances)).reshape(n_grid_cols,n_grid_cols), cmap = 'Greys', zorder = 10, norm = norm,
	linewidth = 0.5)
sns.despine(ax = ax, top = True, bottom = True, left = True, right = True)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, f'{TF}_KD_sim_PCA_arrows.pdf'))

TF = 'ZEB1'
fig, ax = plt.subplots(figsize = (5.24, 4.25))
embedding = scplus_obj.dr_cell[f'{TF}_KD_sim_eRegulons_PCA_iter_0'].to_numpy()
delta_embedding = _project_perturbation_in_embedding(
        scplus_obj,
        original_matrix = scplus_obj.uns[f'{TF}_KD_sim_eRegulon_AUC_iter_4']['Gene_based'],
        perturbed_matrix = scplus_obj.uns[f'{TF}_KD_sim_eRegulon_AUC_iter_0']['Gene_based'],
        reduction_name = f'{TF}_KD_sim_eRegulons_PCA_iter_0')
grid_xy, uv, mask = _calculate_grid_arrows(embedding, delta_embedding, 0.005, n_grid_cols, n_grid_cols, 25, 5)
plot_metadata_given_ax(
        scplus_obj=scplus_obj,
        reduction_name=f'{TF}_KD_sim_eRegulons_PCA_iter_0',
        ax = ax,
        variable = 'MMline',
        color_dictionary = {'MMline':line_color_dict},
        show_label = False,
        show_legend = False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel(None)
ax.set_ylabel(None)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_title(TF)
distances = np.sqrt((uv**2).sum(1))
norm = matplotlib.colors.Normalize(vmin=0.15, vmax=0.5, clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.Greys)
uv[np.logical_or(~mask, np.array(scale(distances)) < 0.15)] = np.nan
ax.streamplot(
        grid_xy.reshape(n_grid_cols,n_grid_cols, 2)[:, :, 0],
        grid_xy.reshape(n_grid_cols,n_grid_cols, 2)[:, :, 1],
        uv.reshape(n_grid_cols,n_grid_cols, 2)[:, :, 0],
        uv.reshape(n_grid_cols,n_grid_cols, 2)[:, :, 1], density = 3, color = np.array(scale(distances)).reshape(n_grid_cols,n_grid_cols), cmap = 'Greys', zorder = 10, norm = norm,
        linewidth = 0.5)
sns.despine(ax = ax, top = True, bottom = True, left = True, right = True)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, f'{TF}_KD_sim_PCA_arrows.pdf'))

```

## Main Figure 5, pannel e: validation of perturbation

```python

## Validation

perturbed_matrices = simulate_perturbation(
        scplus_obj,
        perturbation = {'SOX10': 0},
        regressors = regressors,
        keep_intermediate = True,
        n_iter = 4)

KD_dir = "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/melanoma_mixed_lines/hg38/KD_data/v2/analysis/RNA"
lines_w_SOX10KD = ['MM001', 'MM011', 'MM031', 'MM032', 'MM057', 'MM074', 'MM087']

KD_data = pd.concat(
        [
                pd.read_csv(
                        os.path.join(KD_dir, f'{line}_SOX10_v_NTC_log2fc.tsv'), sep = '\t'
                ).rename(
                        {'log2FoldChange': f'{line}_lfc'}, axis = 1)[[f'{line}_lfc']]
                for line in set(lines_w_SOX10KD) & set(scplus_obj.metadata_cell['MMline'])], axis = 1).fillna(0)

genes_in_both = list(set(regressors.keys()) & set(KD_data.index))

KD_data = KD_data.loc[genes_in_both]

avg_SOX10_sim_kd_ctrl = perturbed_matrices['0'][genes_in_both].groupby(scplus_obj.metadata_cell['MMline'].to_numpy()).mean().T
avg_SOX10_sim_kd = perturbed_matrices['4'][genes_in_both].groupby(scplus_obj.metadata_cell['MMline'].to_numpy()).mean().T
avg_SOX10_sim_kd_ctrl.columns = [n + '_ctrl' for n in avg_SOX10_sim_kd_ctrl.columns]
avg_SOX10_sim_kd.columns = [n + '_kd_sim' for n in avg_SOX10_sim_kd.columns]
SOX10_sim_cts_mtx = pd.concat([avg_SOX10_sim_kd_ctrl, avg_SOX10_sim_kd], axis = 1)

SOX10_sim_cts_mtx[(SOX10_sim_cts_mtx.to_numpy() < 0)] = 0

avg_SOX10_sim_kd[(avg_SOX10_sim_kd.to_numpy() < 0)] = 0


SOX10_sim_lfc = pd.DataFrame(index = avg_SOX10_sim_kd.index, columns = [x.split('_')[0] + '_lfc_sim' for x in avg_SOX10_sim_kd.columns], data =  np.log2((avg_SOX10_sim_kd.to_numpy() + 1e-6) / (avg_SOX10_sim_kd_ctrl.to_numpy() + 1e-6)))

MEL_sim_real = pd.melt(
                pd.concat(
                        [
                SOX10_sim_lfc.loc[list(set(scplus_obj.uns['DEGs']['lineState']['MEL'].index) & set(genes_in_both))],
                KD_data.loc[list(set(scplus_obj.uns['DEGs']['lineState']['MEL'].index) & set(genes_in_both))]],
                axis = 1))
MEL_sim_real['Gene'] = 'MEL'
MES_sim_real = pd.melt(
                pd.concat(
                        [
                SOX10_sim_lfc.loc[list(set(scplus_obj.uns['DEGs']['lineState']['MES'].index) & set(genes_in_both))],
                KD_data.loc[list(set(scplus_obj.uns['DEGs']['lineState']['MES'].index) & set(genes_in_both))]],
                axis = 1))
MES_sim_real['Gene'] = 'MES'
data = pd.concat([MEL_sim_real, MES_sim_real]).reset_index(drop = True)

data.columns = ['line', 'logfc', 'gene']

data = data.loc[~(data['line'] == 'MM099_lfc_sim')]
data = data.loc[~(data['line'] == 'MM029_lfc_sim')]
data = data.loc[~(data['line'] == 'MM047_lfc_sim')]


data['real_sim'] = [f'sim_{y}' if 'sim' in x else f'real_{y}' for x, y in zip(data['line'], data['gene'])]
data['line'] = [x.split('_')[0] for x in data['line']]

fig, ax = plt.subplots(figsize = (3.25, 4.0))
sns.boxplot(data =  data, hue = 'real_sim', y = 'line', x = 'logfc', ax = ax,
            order =sorted(list(set(data['line']))),
	    palette = {'sim_MEL': '#773344', 'real_MEL': '#E3655B', 'sim_MES': '#023436', 'real_MES': '#1B998B'},
            fliersize = 2,
            flierprops = {'markerfacecolor': 'gray', 'markeredgecolor': 'none', 'zorder': 0})
plt.legend([],[], frameon=False)
ax.set_xlim((-3, 3))
for ytick in ax.get_yticklabels():
        ytick.set_color(line_color_dict[ytick.get_text().split('_')[0]])
ax.grid()
ax.set_axisbelow(True)
sns.despine(ax = ax)
fig.tight_layout()
fig.legend(loc='lower left')
fig.savefig(os.path.join(plot_dir, 'boxplot_MEL_MES_SOX10KD_sim.pdf'))

```


## Suplementary Figure S21, pannel b: Mel and Mes gene sets


```python

import pandas as pd
gene_sets = pd.read_csv('/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/melanoma_mixed_lines/WIDMER_2012_gene_sets', sep = '\t')

import dill
import os
import matplotlib
from scenicplus.utils import *
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

plot_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2/plots'
out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2'

scplus_obj = dill.load(open(os.path.join(out_dir, 'scplus_obj.pkl'), 'rb'))

work_dir = '/staging/leuven/stg_00002/lcb/saibar/Projects/PanCancer/1_runs_separateCancers/MMlines_1_sc/scenicplus_output_v10_mmRegions'
gene_ranking = dill.load(open(os.path.join(work_dir, 'gene_ranking.pkl'), 'rb'))

from scenicplus.eregulon_enrichment import *
signatures = {
	'Melanocytic': list(set(gene_sets.query('type == "Proliferative"')['symbol'])), 
	'Mesenchymal': list(set(gene_sets.query('type == "Invasive"')['symbol']))}

sign_auc = signature_enrichment(
	gene_ranking,
	signatures,
	enrichment_type='gene')

data = sign_auc.loc[scplus_obj.dr_cell['eRegulons_PCA'].sort_values('PC_0').index]

scale = lambda X: [(x - min(X)) / (max(X) - min(X)) for x in X]

fig, ax = plt.subplots(figsize = (5.5, 4.5))
ax.plot(
	np.arange(len(scale(data['Melanocytic'].rolling(window = 10).mean().dropna()))),
	scale(data['Melanocytic'].rolling(window = 10).mean().dropna()),
	color = '#6F3C48',
	label = 'Melanocytic')
ax.plot(
	np.arange(len(scale(data['Mesenchymal'].rolling(window = 10).mean().dropna()))),
        scale(data['Mesenchymal'].rolling(window = 10).mean().dropna()),
        color = '#092E30',
        label = 'Mesenchymal')
sns.despine(ax = ax)
fig.legend()
ax.set_xlabel('Cells sorted by PCA1')
ax.set_ylabel('Scaled AUC value')
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'PC0_MEL_MES.pdf'))


sign_auc_scaled = sign_auc.copy()
sign_auc_scaled['Melanocytic'] = scale(sign_auc_scaled['Melanocytic'])
sign_auc_scaled['Mesenchymal'] = scale(sign_auc_scaled['Mesenchymal']) 
order = sign_auc_scaled.max(1).sort_values().index

Red_values = sign_auc_scaled.loc[order, 'Melanocytic'].to_numpy()
Blue_values = sign_auc_scaled.loc[order, 'Mesenchymal'].to_numpy()


fig, ax = plt.subplots(figsize = (5.75, 4.5))
ax.scatter(
	x = scplus_obj.dr_cell['eRegulons_PCA'].loc[order, 'PC_0'].to_numpy(),
	y = scplus_obj.dr_cell['eRegulons_PCA'].loc[order, 'PC_1'].to_numpy(),
	c = np.array([Red_values, np.repeat(0, len(Red_values)), Blue_values]).T,
	s = 1)
sns.despine(ax = ax, top = True, bottom = True, left = True, right = True)
ax.set_xticks([])
ax.set_yticks([])
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'PCA_MEL_MES_AUC.pdf'))
fig.savefig(os.path.join(plot_dir, 'PCA_MEL_MES_AUC.png'), dpi = 300)

```

## Supplementary Figure S21, pannel h: TFs found by SCENIC+, GRaNIE, SCENIC, CellOracle and FigR


```python

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/figures_reviewer_questions/mm_other_eGRN_mehtods'

import os

data_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/extra_methods/mm_lines/extract_data"

import dill
scplus_folder = "/staging/leuven/stg_00002/lcb/sdewin/PhD/Multiomics_pipeline/analysis/melanoma_mixed_lines/hg38_v10_clust_consensus/outs"
scplus_obj = dill.load(
        open(os.path.join(scplus_folder, 'scplus_obj.pkl'), 'rb'))

import pandas as pd
#load data
tool_to_regulon_metadata = {
	'SCENIC': pd.read_table(os.path.join(data_folder, 'SCENIC_regulons.tsv')),
	'CellOracle': pd.read_table(os.path.join(data_folder, 'celloracle_eRegulons.tsv')),
	'FigR': pd.read_table(os.path.join(data_folder, 'figr_eRegulons.tsv')),
	'Pando': pd.read_table(os.path.join(data_folder, 'pando_eRegulons.tsv')),
	'GRaNIE': pd.read_table(os.path.join(data_folder, 'GRaNIE_eRegulons.tsv')),
	'SCENIC+': scplus_obj.uns['eRegulon_metadata_filtered'].copy()}
#seperate positive and negative eRegulons
tool_to_regulon_metadata['SCENIC'].rename({'TF': 'Consensus_name'}, axis = 1, inplace = True)
tool_to_regulon_metadata['CellOracle']['Consensus_name'] = [
	f"{TF}_+" if score >= 0 else f"{TF}_-" for TF, score in tool_to_regulon_metadata['CellOracle'][['TF', 'coef_mean']].to_numpy()]
tool_to_regulon_metadata['FigR']['Consensus_name'] = [
	f"{TF}_+" if score >= 0 else f"{TF}_-" for TF, score in tool_to_regulon_metadata['FigR'][['TF', 'Corr']].to_numpy()]
tool_to_regulon_metadata['Pando']['Consensus_name'] = [
	 f"{TF}_+" if score >= 0 else f"{TF}_-" for TF, score in tool_to_regulon_metadata['Pando'][['TF', 'estimate']].to_numpy()]
tool_to_regulon_metadata['GRaNIE']['Consensus_name'] = [
	f"{TF}_+" if score == 'pos' else f"{TF}_-" for TF, score in tool_to_regulon_metadata['GRaNIE'][['TF', 'TF_peak.fdr_direction']].to_numpy()]
tool_to_regulon_metadata['SCENIC+']['Consensus_name'] = [x.split('_(')[0] for x in tool_to_regulon_metadata['SCENIC+']['Region_signature_name']]
tool_to_regulon_metadata['SCENIC+'].rename({'Gene': 'gene', 'Region': 'region'}, axis = 1, inplace = True)


def remove_second_sign(x):
        if 'extended' not in x:
                TF, first, second, n = x.split('_')
                return f'{TF}_{first}_{n}'
        else:
                TF, extended, first, second, n = x.split('_')
                return f'{TF}_{extended}_{first}_{n}'

tool_to_regulon_metadata['SCENIC+']['Gene_signature_name'] = [remove_second_sign(x) for x in tool_to_regulon_metadata['SCENIC+']['Gene_signature_name']]
tool_to_regulon_metadata['SCENIC+']['Region_signature_name'] = [remove_second_sign(x) for x in tool_to_regulon_metadata['SCENIC+']['Region_signature_name']]

selected_regulons = scplus_obj.uns['selected_eRegulons']['Region_based']

tool_to_regulon_metadata['SCENIC+'] = tool_to_regulon_metadata['SCENIC+'].query('Region_signature_name in @selected_regulons')

min_target_genes = 5
for tool in tool_to_regulon_metadata.keys():
	eRegulon_and_counts = tool_to_regulon_metadata[tool].groupby('Consensus_name')['gene'].apply(lambda x: len(set(x)))
	eRegulons_to_keep = eRegulon_and_counts.index[eRegulon_and_counts >= min_target_genes]
	tool_to_regulon_metadata[tool] = tool_to_regulon_metadata[tool].query('Consensus_name in @eRegulons_to_keep')

tool_to_regulon_metadata['SCENIC']['TF'] = [x.split('_')[0] for x in tool_to_regulon_metadata['SCENIC']['Consensus_name']]

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
	TFs_for_tools.loc[set([x.split('_')[0] for x in scplus_obj.uns['selected_eRegulons']['Gene_based']])].T,
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
g.savefig(os.path.join(out_dir, "scplus_TFs_found_by_methods.pdf"))
g.savefig(os.path.join(out_dir, "scplus_TFs_found_by_methods.png"), dpi = 300)

```

## Supplementary Figure S21 pannel i & j: enhancer activity

```python

out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/figures_reviewer_questions/mm_lines_enhancer_activity'

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

wdir = "/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/mm_lines_target_regions"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

f_scplus_obj = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/Figure_MM_lines/v2/scplus_obj.pkl'

import dill
scplus_obj = dill.load(open(f_scplus_obj, 'rb'))

import os
cheq_seq_activity = pd.read_table(
	os.path.join(wdir, 'CHEQ-seq_ATAC_activity.tsv'), index_col = 0).drop('phenotype', axis = 1)

starr_seq_activity = pd.read_table(
	os.path.join(wdir, 'STARR-seq_ATAC_activity.tsv'), index_col = 0).drop('phenotype', axis = 1)

cheq_seq_5_prime_ac = pd.read_table(
	os.path.join(wdir, 'CHEQ-seq_5_H3K27ac_activity.tsv'), index_col = 0).drop('phenotype', axis = 1)

cheq_seq_intron_ac = pd.read_table(
	os.path.join(wdir, 'CHEQ-seq_intron_H3K27ac_activity.tsv'), index_col = 0).drop('phenotype', axis = 1)

region_definitions = pd.read_table(
	os.path.join(wdir, 'mauduit_consensus_overlap.bed'), header = None).drop([0,1,2,4,8,9,10,11], axis = 1)
region_definitions.columns = ['name', 'chrom', 'start', 'end']
region_definitions['region_name'] = [f"{chrom}:{start}-{end}" for chrom, start, end in region_definitions[['chrom', 'start', 'end']].to_numpy()]

region_definitions['max_ac_5'] = [cheq_seq_5_prime_ac.loc[n].max() if n in cheq_seq_5_prime_ac.index else np.nan for n in region_definitions['name'].to_list()]
region_definitions['max_ac_intron'] = [cheq_seq_intron_ac.loc[n].max() if n in cheq_seq_intron_ac.index else np.nan for n in region_definitions['name'].to_list()]
region_definitions['max_cheq_atac'] = [cheq_seq_activity.loc[n].max() if n in cheq_seq_activity.index else np.nan for n in region_definitions['name'].to_list()]
region_definitions['max_starr_atac'] = [starr_seq_activity.loc[n].max() if n in starr_seq_activity.index else np.nan for n in region_definitions['name'].to_list()]

for c in region_definitions.columns:
	if c.startswith('max_'):
		region_definitions['log_'+c] = np.log10(region_definitions[c])

region_definitions = region_definitions.loc[[x in scplus_obj.region_names for x in region_definitions['region_name']]]

region_definitions['in_eRegulon'] = [region in list(scplus_obj.uns['eRegulon_metadata_filtered']['Region']) for region in region_definitions['region_name']]


db_fname = '/staging/leuven/stg_00002/lcb/saibar/Projects/PanCancer/1_runs_separateCancers/MMlines_1_sc/pycistopic_dbs_v10/mel_9sc.regions_vs_motifs.rankings.feather'

from scenicplus.triplet_score import calculate_TF_to_region_score, calculate_triplet_score

calculate_TF_to_region_score(
	scplus_obj, ctx_db_fname = db_fname, eRegulon_metadata_key = 'eRegulon_metadata_filtered')
calculate_triplet_score(
	scplus_obj, eRegulon_metadata_key = 'eRegulon_metadata_filtered')

region_definitions.loc[region_definitions['in_eRegulon'], 'Tf_to_region_score'] = [
	scplus_obj.uns['eRegulon_metadata_filtered'].set_index('Region').loc[region, 'TF_to_region_max_rank'].min()
	for region in region_definitions.loc[region_definitions['in_eRegulon'], 'region_name']]

region_definitions.loc[region_definitions['in_eRegulon'], 'triplet_score'] = [
	scplus_obj.uns['eRegulon_metadata_filtered'].set_index('Region').loc[region, 'triplet_score'].min()
	for region in region_definitions.loc[region_definitions['in_eRegulon'], 'region_name']]

plotting_parameters = {
	'data': region_definitions[['in_eRegulon', 'log_max_starr_atac']].dropna(),
	'x': 'in_eRegulon',
	'y': 'log_max_starr_atac'}

data = plotting_parameters['data']

from scipy.stats import mannwhitneyu
from itertools import combinations
pairs = list(combinations(set(region_definitions['in_eRegulon']), 2))

p_values = [
	mannwhitneyu(data.loc[data['in_eRegulon'], 'log_max_starr_atac'], data.loc[~data['in_eRegulon'], 'log_max_starr_atac']).pvalue
	for t1, t2 in pairs]

p_values_formated = [f'p={p:.1e}' for p in p_values]


from statannotations.Annotator import Annotator

import seaborn as sns

from matplotlib import ticker as mticker
fig, ax = plt.subplots(figsize = (4, 4))
sns.violinplot(
	inner = None,
	palette = {True: 'gray', False: 'darkgray'},
	ax = ax,
	**plotting_parameters)
sns.stripplot(
	color="black", size = 4,
	ax =ax, alpha = 0.4,
	**plotting_parameters)
ax.set_xlabel('eRegulon target')
ax.set_ylabel('STARR-seq LogFoldChange')
ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
ymin, ymax = ax.get_ylim()
tick_range = np.arange(np.floor(ymin), ymax)
ax.yaxis.set_ticks(tick_range)
ax.yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)
annotator = Annotator(ax, pairs, **plotting_parameters)
annotator.set_custom_annotations(p_values_formated)
annotator.annotate()
fig.tight_layout()
sns.despine(ax = ax)
fig.savefig(os.path.join(out_dir, 'activity_in_eRegulon_vs_out.pdf'))
fig.savefig(os.path.join(out_dir, 'activity_in_eRegulon_vs_out.png'))
fig1 = fig


from adjustText import adjust_text
from scipy.stats import spearmanr
fig, ax = plt.subplots(figsize = (4, 4))
data = region_definitions.loc[region_definitions['in_eRegulon'], ['name', 'max_starr_atac', 'region_name', 'triplet_score']].dropna()
ax.scatter(
	x = data['triplet_score'].to_numpy(),
	y = data['max_starr_atac'].to_numpy(),
	s = 10, color = '#2b2b2b')
texts = [plt.text(_x, _y, t) for _x, _y, t in data[['triplet_score', 'max_starr_atac', 'name']].to_numpy()]
adjust_text(texts)
ax.set_xlabel('Triplet Score Rank')
ax.set_ylabel('STARR-seq FoldChange')
ax.set_title(f"Spearmanr = {np.round(spearmanr(data['triplet_score'].to_numpy(), data['max_starr_atac'].to_numpy()).correlation, 3)}")
sns.despine(ax = ax)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'scat_triplet_max_starr.pdf'))
fig.savefig(os.path.join(out_dir, 'scat_triplet_max_starr.png'))
fig2 = fig

import quilter
fig = fig1 + fig2
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'vln_scat.pdf'))
fig.savefig(os.path.join(out_dir, 'vln_scat.png'))

```
