# Code for generating the pannels of figure 4 and accompanying supplementary figures

## Main Figure 4, pannel b: comparison of triplet ranking with chip-seq, Enformer and STARR-seq

```python

import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import quilter
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/figure_target_regions_DPCL'

#Load scplus_obj
import dill
f_scplus_obj = '/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/DPCL/scenicplus_final_autoreg/grnboost/scplus_obj.pkl'
scplus_obj = dill.load(open(f_scplus_obj, 'rb'))

#calculate triplet score and TF-to-region score
from mudata import AnnData
import mudata

adata_max_rank = mudata.read('/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/triplet_score/DPCL_adata_max_rank.h5ad')
df_max_rank = adata_max_rank.to_df()

eRegulon_metadata_key = 'eRegulon_metadata'
TF_region_iter = scplus_obj.uns[eRegulon_metadata_key][['TF', 'Region']].to_numpy()
TF_to_region_score = [df_max_rank.loc[region, TF] for TF, region in TF_region_iter]
scplus_obj.uns[eRegulon_metadata_key]['TF_to_region_max_rank'] = TF_to_region_score

from scenicplus.triplet_score import _rank_scores_and_assign_random_ranking_in_range_for_ties, _calculate_cross_species_rank_ratio_with_order_statistics
TF2G_score_key = 'TF2G_importance'
R2G_score_key = 'R2G_importance'
TF2R_score_key = 'TF_to_region_max_rank'
TF2G_score = scplus_obj.uns[eRegulon_metadata_key][TF2G_score_key].to_numpy()
R2G_score = scplus_obj.uns[eRegulon_metadata_key][R2G_score_key].to_numpy()
TF2R_score = scplus_obj.uns[eRegulon_metadata_key][TF2R_score_key].to_numpy()
#rank the scores
TF2G_rank = _rank_scores_and_assign_random_ranking_in_range_for_ties(TF2G_score)
R2G_rank = _rank_scores_and_assign_random_ranking_in_range_for_ties(R2G_score)
TF2R_rank = _rank_scores_and_assign_random_ranking_in_range_for_ties(-TF2R_score) #negate because lower score is better

import numpy as np
TF2G_rank_ratio = (TF2G_rank.astype(np.float64) + 1) / TF2G_rank.shape[0]
R2G_rank_ratio = (R2G_rank.astype(np.float64) + 1) / R2G_rank.shape[0]
TF2R_rank_ratio = (TF2R_rank.astype(np.float64) + 1) / TF2R_rank.shape[0]

rank_ratios = np.array([TF2G_rank_ratio, R2G_rank_ratio, TF2R_rank_ratio])
aggregated_rank = np.zeros((rank_ratios.shape[1],), dtype = np.float64)
for i in range(rank_ratios.shape[1]):
    aggregated_rank[i] = _calculate_cross_species_rank_ratio_with_order_statistics(rank_ratios[:, i])

scplus_obj.uns[eRegulon_metadata_key]['triplet_score'] = aggregated_rank.argsort().argsort()

#load STARR-seq data
starr_data_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/DPCL_starr_comparison' 
import pandas as pd
import os
line_to_starr_seq = {
	'K562': pd.read_table(os.path.join(starr_data_dir, 'consensus_peaks_ENCFF045TVA_K562.bed'), header = None),
	'HepG2': pd.read_table(os.path.join(starr_data_dir, 'consensus_peaks_ENCFF047LDJ_HepG2.bed'), header = None),
	'HCT116': pd.read_table(os.path.join(starr_data_dir, 'consensus_peaks_ENCFF428KHI_HCT116.bed'), header = None),
	'MCF7': pd.read_table(os.path.join(starr_data_dir, 'consensus_peaks_ENCFF826BPU_MCF7.bed'), header = None)}

for line in line_to_starr_seq.keys():
	line_to_starr_seq[line].columns = [
		'Chrom_1', 'Start_1', 'End_1', 'consensus_region_name',
		'Chrom_2', 'Start_2', 'End_2', 'peak',
		'score', 'strand', 'log2FoldChange', 'inputCount', 'outputCount', 'minusLog10PValue', 'minusLog10QValue', '?']

def fast_is_in(a, b: set):
	assert type(b) == set, "b should be of type set"
	is_in = np.zeros(len(a), dtype = bool)
	for i, val in enumerate(a):
		is_in[i] = val in b
	return is_in
	
eRegulon_regions = set(scplus_obj.uns['eRegulon_metadata']['Region'].to_list())

for line in line_to_starr_seq.keys():
	line_to_starr_seq[line]['in_eRegulon'] = fast_is_in(line_to_starr_seq[line]['consensus_region_name'].to_list(), eRegulon_regions)

for line in line_to_starr_seq.keys():
	line_to_starr_seq[line]['log_score'] = np.log10(line_to_starr_seq[line]['score'])

line_to_starr_all = pd.concat([line_to_starr_seq[line][['consensus_region_name', 'log2FoldChange', 'in_eRegulon']] for line in line_to_starr_seq.keys()])
line_to_starr_max = line_to_starr_all.groupby('consensus_region_name').agg({'log2FoldChange': max, 'in_eRegulon': sum})
line_to_starr_max['in_eRegulon'] = line_to_starr_max['in_eRegulon'].astype(bool)

import pickle
chip_seq_scores = pickle.load(open('/lustre1/project/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/DPCL_chip_seq_comparison/max_chip_DPCL.pkl', 'rb'))

chip_seq_scores = chip_seq_scores.max(1)
chip_seq_scores = pd.DataFrame({'max_chip':chip_seq_scores})


enformer_scores = pickle.load(open('/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/enformer_comparison/df_consensus_peaks_w_scores.pkl', 'rb'))
enformer_scores['Name'] = enformer_scores.index

TFs_of_interest = list(set(scplus_obj.uns['eRegulon_metadata']['TF']))
DPCLs = ['MCF-7', 'HepG2', 'PC-3', 'GM12878', 'Panc1', 'IMR-90', 'HCT116', 'K562']
targets_txt = 'https://raw.githubusercontent.com/calico/basenji/0.5/manuscripts/cross2020/targets_human.txt'
df_targets = pd.read_csv(targets_txt, sep='\t')
descriptions_of_interest = [
        (idx, x) for x, idx in zip(df_targets['description'], df_targets.index)
        if x.split(':')[0] == 'CHIP' and x.split(':')[-1] in DPCLs and x.split(':')[1] in TFs_of_interest]

idx_of_interest = np.array(sorted([d[0] for d in descriptions_of_interest]))

enformer_scores['max_track_of_interest'] = np.array(enformer_scores['Score'].tolist())[:, idx_of_interest].max(1)

enformer_scores_max = enformer_scores.groupby('Name')['max_track_of_interest'].max()

region_max_triplet_tf_2_region = scplus_obj.uns['eRegulon_metadata'][['Region', 'TF_to_region_max_rank', 'triplet_score']].groupby('Region').min()

enformer_score_and_scplus = pd.merge(enformer_scores_max, region_max_triplet_tf_2_region, left_index = True, right_index = True)
chip_score_and_scplus = pd.merge(chip_seq_scores, region_max_triplet_tf_2_region,  left_index = True, right_index = True)
starr_score_and_scplus = pd.merge(line_to_starr_max, region_max_triplet_tf_2_region, left_index = True, right_index = True)

from itertools import combinations
from scenicplus.utils import p_adjust_bh
from scipy.stats import mannwhitneyu
from statannotations.Annotator import Annotator
from matplotlib import ticker as mticker

def calculate_pvals(plot_params, test):
	pairs = list(combinations(set(plot_params['data'][plot_params['x']]), 2 ))
	p_values = [
		test(plot_params['data'].loc[plot_params['data'][plot_params['x']] == p1, plot_params['y']].to_numpy(), 
		     plot_params['data'].loc[plot_params['data'][plot_params['x']] == p2, plot_params['y']].to_numpy()).pvalue
		for p1, p2 in pairs]
	p_values = p_adjust_bh(p_values)
	formatted_pvalues = [f'p={pvalue:.1e}' if pvalue !=0 else f"****" for pvalue in p_values]
	return pairs, formatted_pvalues

def set_log_scale(ax):
	ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
	ymin, ymax = ax.get_ylim()
	tick_range = np.arange(np.floor(ymin), ymax)
	ax.yaxis.set_ticks(tick_range)
	ax.yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)

fig_corr_triplet, axs = plt.subplots(nrows = 3, ncols = 1, figsize = (4.75, 4.75))
sns.set_palette('BrBG')
to_show = 'max_chip'
top = chip_score_and_scplus.nsmallest(n = int(len(chip_score_and_scplus) * 0.1), columns = 'triplet_score')
top['Triplet Score'] = 'Top 10%'
bot = chip_score_and_scplus.nlargest( n = int(len(chip_score_and_scplus) * 0.1), columns = 'triplet_score')
bot['Triplet Score'] = 'Bottom 10%'
_not = chip_seq_scores.loc[list(set(chip_seq_scores.index) - set(chip_score_and_scplus.index))]
_not['Triplet Score'] = 'Not in eRegulon'
data = pd.concat([top, bot, _not])
data[to_show] = np.log10(data[to_show])
plot_params = {'data': data, 'x': 'Triplet Score', 'y': to_show}
sns.boxplot(ax = axs[0], **plot_params)
pairs, pvalues = calculate_pvals(plot_params, mannwhitneyu)
annotator = Annotator(axs[0], pairs[0:2][::-1], **plot_params)
annotator.configure(loc = 'outside')
annotator.set_custom_annotations(pvalues[0:2][::-1])
#annotator.annotate()
set_log_scale(axs[0])
axs[0].grid(color = 'lightgray')
axs[0].set_axisbelow(True)
axs[0].set_xlabel(None)
axs[0].set_ylabel('ChIP-seq Coverage')
axs[0].set_ylim((-0.5, 4.5))
#sns.despine(ax = axs[1, 0])
to_show = 'max_track_of_interest'
top = enformer_score_and_scplus.nsmallest(n = int(len(enformer_score_and_scplus) * 0.1), columns = 'triplet_score')
top['Triplet Score'] = 'Top 10%'
bot = enformer_score_and_scplus.nlargest( n = int(len(enformer_score_and_scplus) * 0.1), columns = 'triplet_score')
bot['Triplet Score'] = 'Bottom 10%'
_not = enformer_scores_max.loc[list(set(enformer_scores_max.index) - set(enformer_score_and_scplus.index))]
_not = pd.DataFrame({to_show:_not})
_not['Triplet Score'] = 'Not in eRegulon'
data = pd.concat([top, bot, _not])
data[to_show] = np.log10(data[to_show])
plot_params = {'data': data, 'x': 'Triplet Score', 'y': to_show}
sns.boxplot(ax = axs[1], **plot_params)
pairs, pvalues = calculate_pvals(plot_params, mannwhitneyu)
annotator = Annotator(axs[1], pairs[0:2][::-1], **plot_params)
annotator.configure(loc = 'outside')
annotator.set_custom_annotations(pvalues[0:2][::-1])
#annotator.annotate()
set_log_scale(axs[1])
axs[1].grid(color = 'lightgray')
axs[1].set_axisbelow(True)
axs[1].set_xlabel(None)
axs[1].set_ylabel('Enformer Score')
axs[1].set_ylim((-1, 2.5))
#sns.despine(ax = axs[1, 1])
to_show = 'log2FoldChange'
top = starr_score_and_scplus.nsmallest(n = int(len(starr_score_and_scplus) * 0.1), columns = 'triplet_score')
top['Triplet Score'] = 'Top 10%'
bot = starr_score_and_scplus.nlargest( n = int(len(starr_score_and_scplus) * 0.1), columns = 'triplet_score')
bot['Triplet Score'] = 'Bottom 10%'
_not = line_to_starr_max.loc[list(set(line_to_starr_max.index) - set(starr_score_and_scplus.index))]
_not['Triplet Score'] = 'Not in eRegulon'
data = pd.concat([top, bot, _not])
plot_params = {'data': data, 'x': 'Triplet Score', 'y': to_show}
sns.boxplot(ax = axs[2], **plot_params)
pairs, pvalues = calculate_pvals(plot_params, mannwhitneyu)
annotator = Annotator(axs[2], pairs[0:2][::-1], **plot_params)
annotator.configure(loc = 'outside')
annotator.set_custom_annotations(pvalues[0:2][::-1])
#annotator.annotate()
set_log_scale(axs[2])
axs[2].grid(color = 'lightgray')
axs[2].set_axisbelow(True)
axs[2].set_xlabel(None)
axs[2].set_ylabel('ChIP-seq Coverage')
axs[2].set_ylim((0,4))
#sns.despine(ax = axs[1, 2])
fig_corr_triplet.tight_layout()
fig_corr_triplet.savefig(os.path.join(out_dir, 'trplt_vs_chip_enf_starr_boxplot.pdf'))
fig_corr_triplet.savefig(os.path.join(out_dir, 'trplt_vs_chip_enf_starr_boxplot.png'))

```

## Main Figure 4, pannel c: ChIP-seq coverage and Enformer predictions on HNF4, FOXA2 and CEBPB target regions

```python

target_regions_directory = '/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/DPCL/methods_benchmark/regulons_bed_hq/'
import os
tool_to_target_regions_path = {
        'GRaNIE': os.path.join(target_regions_directory, 'GRaNIE_regions.tsv'),
        'Pando': os.path.join(target_regions_directory, 'Pando_regions.tsv'),
        'SCENIC+': os.path.join(target_regions_directory, 'SCENIC+_regions.tsv'),
        'CellOracle': os.path.join(target_regions_directory, 'CellOracle_regions.tsv')}

tool_to_target_regions_df = {
	tool: pd.read_table(
		tool_to_target_regions_path[tool], header = None, names = ['Chromosome', 'Start', 'End', 'TF'])
	for tool in tool_to_target_regions_path.keys()}

tool_to_TF_to_target_regions = {}
for tool in tool_to_target_regions_df.keys():
	tool_to_target_regions_df[tool].index = [f"{chrom}:{start}-{end}" for chrom, start, end in tool_to_target_regions_df[tool][['Chromosome', 'Start', 'End']].itertuples(index = False)]
	tool_to_target_regions_df[tool]['name'] = tool_to_target_regions_df[tool].index
	tool_to_target_regions_df[tool]['TF'] = [x.split('_')[0] for x in tool_to_target_regions_df[tool]['TF']]
	tmp = dict(tuple(tool_to_target_regions_df[tool].groupby('TF')['name']))
	tool_to_TF_to_target_regions[tool] = {
		TF: tmp[TF].to_numpy()
		for TF in tmp.keys()}

flatten_list = lambda l: [item for sublist in l for item in sublist]
TFs_of_interest = list(set(flatten_list([tool_to_TF_to_target_regions[tool].keys() for tool in tool_to_TF_to_target_regions.keys()])))

targets_txt = 'https://raw.githubusercontent.com/calico/basenji/0.5/manuscripts/cross2020/targets_human.txt'
DPCLs = ['MCF-7', 'HepG2', 'PC-3', 'GM12878', 'Panc1', 'IMR-90', 'HCT116', 'K562']
df_targets = pd.read_csv(targets_txt, sep='\t')
descriptions_of_interest = [
        (idx, x) for x, idx in zip(df_targets['description'], df_targets.index)
        if x.split(':')[0] == 'CHIP' and x.split(':')[-1] in DPCLs and x.split(':')[1] in TFs_of_interest]

TFs = np.array([x[1].split(':')[1] for x in descriptions_of_interest])
idx_a = np.array([x[0] for x in descriptions_of_interest])

TFs_to_label = ['HNF4A', 'FOXA2', 'CEBPB']

TF_to_bw = {
	'HNF4A': '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/enformer_comparison/ENCFF080FZD.bigWig',
	'FOXA1': '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/enformer_comparison/ENCFF529CBU.bigWig',
	'GATA2': '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/enformer_comparison/ENCFF781FFF.bigWig',
	'STAT5A': '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/enformer_comparison/ENCFF249CGE.bigWig',
	'CEBPB': '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/enformer_comparison/ENCFF003HJB.bigWig',
	'FOXA2': '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/enformer_comparison/ENCFF626IVY.bigWig'
}

all_scores = np.array(enformer_scores['Score'].tolist())

from scenicplus.utils import Groupby
grouper = Groupby(enformer_scores['Name'])
max_score_across_enformer_bins = np.zeros(grouper.n_keys)

def reshape(x, l):
	if x.shape[0] >= l:
		return x[0:l]
	else:
		if x.shape[0] % 2 == 0 and l % 2 == 0:
			return np.pad(x, (l - x.shape[0] - 1), 'edge')
		else:
			return np.pad(x, (l - x.shape[0], 0), 'edge')

def region_name_to_range(r):
	chrom, start, end = r.replace(':', '-').split('-')
	return chrom, int(start), int(end)

import pyBigWig
from tqdm import tqdm

#heatmap w ChIP-seq
nbins = 5
kernel_size = 20
kernel = np.ones(kernel_size) / kernel_size
figs = {}
for TF in TFs_to_label:
	kernel_size = int(len(tool_to_TF_to_target_regions['SCENIC+'][TF]) / 50) 
	kernel = np.ones(kernel_size) / kernel_size
	print(TF)
	idx = idx_a[TFs == TF]
	i = idx[0]
	triplet_score_regions = scplus_obj.uns[eRegulon_metadata_key].query('TF == @TF')[['Region', 'TF_to_region_max_rank']].groupby('Region').max()
	scores_over_consensus_regions_TF = all_scores[:, i][np.array([reshape(x, 5) for x in grouper.indices], dtype = int)]
	percentiles = np.percentile(scores_over_consensus_regions_TF.max(1), [1, 99], interpolation='nearest')
	upper_percentile_region = grouper.keys[np.where(scores_over_consensus_regions_TF.max(1) == percentiles[1])][0]
	lower_percentile_region = grouper.keys[np.where(scores_over_consensus_regions_TF.max(1) == percentiles[0])][0]
	percentiles = np.percentile(scores_over_consensus_regions_TF.max(1), [1, 99])
	y_trues = []
	for tool in tool_to_TF_to_target_regions.keys():
		if TF in tool_to_TF_to_target_regions[tool].keys():
			regions_w_scores = set(tool_to_TF_to_target_regions[tool][TF]) & set(grouper.keys)
			y_true = fast_is_in(grouper.keys, regions_w_scores)
			y_trues.append(y_true * 1)
		else:
			y_trues.append(np.zeros(len(grouper.keys)))
	y_trues = np.array(y_trues).T
	idx_to_keep = y_trues.sum(1) >= 1
	y_trues = y_trues[idx_to_keep]
	scores_over_consensus_regions_TF = scores_over_consensus_regions_TF[idx_to_keep]
	sorted_idx = np.argsort(-scores_over_consensus_regions_TF.max(1))
	y_trues = y_trues[sorted_idx]
	scores_over_consensus_regions_TF = scores_over_consensus_regions_TF[sorted_idx]
	bw = pyBigWig.open(TF_to_bw[TF])
	scores_over_consensus_regions_TF_ChIP = np.array(
		[
			bw.stats(*region_name_to_range(grouper.keys[idx_to_keep][sorted_idx][i]), nBins = nbins) 
			for i in tqdm(range(sum(idx_to_keep)), total = sum(idx_to_keep))],
		dtype = float)
	upper_percentile_chip = max(bw.stats(*region_name_to_range(upper_percentile_region), nBins = nbins))
	lower_percentile_chip = max(bw.stats(*region_name_to_range(lower_percentile_region), nBins = nbins))
	fig, axs = plt.subplots(figsize = ((122.5 / 10) / len(TFs_to_label), 4.75), ncols = 4, gridspec_kw = {'width_ratios': [4, 4, 2, 2]}, sharey = True, constrained_layout = True)
	sns.heatmap(
		scores_over_consensus_regions_TF_ChIP, ax = axs[0],
		cbar = True, yticklabels = False, cmap = 'Spectral_r', xticklabels = [0, 100, 100 * 2, 100 * 3, 100 * 4],
		vmin = lower_percentile_chip, vmax = upper_percentile_chip,
		cbar_kws = dict(use_gridspec=True, location="left", shrink = 0.2, label = 'ChIP-seq Coverage', fraction = 0.05, pad = 0.01))
	sns.heatmap(
		scores_over_consensus_regions_TF, ax = axs[1],
		cbar = True, yticklabels = False, cmap = 'Spectral_r', xticklabels = [0, 128, 128 * 2, 128 * 3, 128 * 4],
		vmin = percentiles[0], vmax = percentiles[1],
		cbar_kws = dict(use_gridspec=True, location="left", shrink = 0.2, label = 'Enformer Score', fraction = 0.05, pad = 0.01))
	sns.heatmap(
		y_trues, ax = axs[2], cmap = 'binary', vmin = 0, vmax = 1, yticklabels = False, cbar = False, xticklabels = tool_to_TF_to_target_regions.keys())
	yx = np.array(
		[
			(i, triplet_score_regions.loc[region]['TF_to_region_max_rank'])
			for i, region in enumerate(grouper.keys[idx_to_keep][sorted_idx]) if region in triplet_score_regions.index
		]).T
	convolved_TF2R = np.convolve(yx[1, :] + 1, kernel, mode = 'full')[0:yx.shape[1]]
	axs[3].scatter(y = yx[0, :], x = convolved_TF2R, s = 2, c = 'black')
	axs[3].set_xlabel('TF2R ranking', fontsize = 'x-small')
	axs[3].set_xticks([int(min(convolved_TF2R)), int(max(convolved_TF2R))])
	sns.despine(ax = axs[3])
	fig.suptitle(TF)
	figs[TF] = fig
	fig.savefig(os.path.join(out_dir, f'{TF}_heatmap_enf_chip.pdf'))

fig_heatmap = figs['HNF4A'] + figs['FOXA1'] + figs['GATA2'] + figs['STAT5A']
fig_heatmap.set_size_inches((172.5 / 10, 47.5 / 10))
fig_heatmap.tight_layout()
fig_heatmap.savefig(os.path.join(out_dir, 'heatmap_enf_chip.png'))
fig_heatmap.savefig(os.path.join(out_dir, 'heatmap_enf_chip.pdf'))

scplus_obj.uns['eRegulon_metadata'].to_csv(os.path.join(out_dir, 'eRegulon_metadata.tsv'), sep = '\t', header = True, index = False)

```

## Main Figure 4, pannel f and Supplementary Figure S20: Enformer importance scores on best HNF4 target region

```python

out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/figure_target_regions_DPCL'

import os
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = '0' # CHOOSE GPU HERE
import tensorflow as tf
from tensorflow.python.keras.backend import set_session

config = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=8, 
                        inter_op_parallelism_threads=8, 
                        allow_soft_placement=True,
                        device_count = {'CPU': 8})
session = tf.compat.v1.Session(config=config)
set_session(session)
tf.config.threading.set_intra_op_parallelism_threads(8)
tf.config.threading.set_inter_op_parallelism_threads(8)

import sys
sys.path.insert(0, '/lustre1/project/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/enformer_deep_explainers')
import enformer_funcs as ef

import kipoiseq
from kipoiseq import Interval
import pyfaidx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from pybedtools import BedTool
import numpy as np
from tqdm import tqdm

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'

SEQUENCE_LENGTH = 393_216

# Define reference genome and load Enformer model
species = 'human'
genomefile_path_hg = "/staging/leuven/stg_00002/lcb/resources/human/hg38/hg38.fa"
fasta_extractor = ef.FastaStringExtractor(genomefile_path_hg)
model = ef.Enformer("https://tfhub.dev/deepmind/enformer/1")

targets_txt = 'https://raw.githubusercontent.com/calico/basenji/0.5/manuscripts/cross2020/targets_human.txt'
DPCLs = ['MCF-7', 'HepG2', 'PC-3', 'GM12878', 'Panc1', 'IMR-90', 'HCT116', 'K562']
df_targets = pd.read_csv(targets_txt, sep='\t')
descriptions_of_interest = [
        (idx, x) for x, idx in zip(df_targets['description'], df_targets.index)
        if x.split(':')[0] == 'CHIP' and x.split(':')[-1] in DPCLs]


def region_name_to_chrom_start_end(r):
	chrom, start, end = r.replace(':', '-').split('-')
	return chrom, int(start), int(end)

import pandas as pd

eRegulon_metadata = pd.read_csv(
	os.path.join(out_dir, 'eRegulon_metadata.tsv'),
	sep = '\t')

selected_region_HEPG2 = eRegulon_metadata.query('TF == "HNF4A"').sort_values('triplet_score')['Region'].tolist()[0]
TFs_targettting_selected_region = set(eRegulon_metadata.query('Region == @selected_region_HEPG2')['TF'])
target_genes_of_selected_region = set(eRegulon_metadata.query('Region == @selected_region_HEPG2')['Gene'])

track_nr = [x for x in descriptions_of_interest if x[1].split(':')[1] == 'HNF4A' and x[1].split(':')[2] == 'HepG2'][0][0]

def plot_deepexplainer_givenax_givenregion(region, track_nr, ax):
	chrom, start, end = region_name_to_chrom_start_end(region)
	target_interval = kipoiseq.Interval(chrom, start - 57094, end + 57094)
	sequence_one_hot = ef.one_hot_encode(fasta_extractor.extract(target_interval.resize(SEQUENCE_LENGTH)))
	predictions = model.predict_on_batch(sequence_one_hot[np.newaxis])[species][0]
	target_mask = np.zeros_like(predictions)
	for idx in [447, 448, 449]: #only for the middle part of the sequence
	    target_mask[idx, track_nr] = 1
	    #target_mask[idx, 183] = 1 # you can add more than one track
	contribution_scores = model.contribution_input_grad(sequence_one_hot.astype(np.float32), target_mask, output_head = species).numpy()
	pooled_contribution_scores = ef.tf.nn.avg_pool1d(np.abs(contribution_scores)[np.newaxis, :, np.newaxis], 128, 128, 'VALID')[0, :, 0].numpy()[1088:-1088]
	# Select the center 500bp region and reshape it to (500,4) for practical purposes
	contrib_scores = np.repeat(contribution_scores[196608-250:196608+250],4) 
	contrib_scores = np.reshape(contrib_scores, (500,4))
	ef.plot_deepexplainer_givenax(contribution_scores = contrib_scores*sequence_one_hot[196608-250:196608+250], ax = ax)

def calculate_saturation_mutagenesis_for_region(region, track_nr):
	chrom, start, end = region_name_to_chrom_start_end(region)
	target_interval = kipoiseq.Interval(chrom, start - 57094, end + 57094)
	sequence_one_hot = ef.one_hot_encode(fasta_extractor.extract(target_interval.resize(SEQUENCE_LENGTH)))
	original_pred = np.mean(model.predict_on_batch(sequence_one_hot[np.newaxis])[species][0][447:450, track_nr])
	x = ef.create_saturation_mutagenesis_x(sequence_one_hot[196608-250:196608+250])
	x['Prediction'] = []
	for i in tqdm(range(500*3), total = 500 * 3):
		seq = np.vstack((np.vstack((sequence_one_hot[:196608-250],x['X'][i])),sequence_one_hot[196608+250:]))
		x['Prediction'].append((model.predict_on_batch(seq[np.newaxis])[species][0]))
	return original_pred, x

def plot_mutagenesis_givenax_fast(mut_scores, original_prediction, track_nr, ax, title='Enfomer ISM', dotsize = 1):
    seq_shape = (500,4)
    mutagenesis_X = mut_scores
    arr_a = np.zeros(seq_shape[0])
    arr_c = np.zeros(seq_shape[0])
    arr_g = np.zeros(seq_shape[0])
    arr_t = np.zeros(seq_shape[0])
    mut_preds = np.array(mutagenesis_X['Prediction'])[:,447:450,track_nr]
    delta_pred = original_prediction - np.mean(mut_preds, axis=1)
    for i,mut in enumerate(mutagenesis_X["ids"]):
        if mut.endswith("A"):
            arr_a[int(mut.split("_")[0])]=delta_pred[i]
        if mut.endswith("C"):
            arr_c[int(mut.split("_")[0])]=delta_pred[i]
        if mut.endswith("G"):
            arr_g[int(mut.split("_")[0])]=delta_pred[i]
        if mut.endswith("T"):
            arr_t[int(mut.split("_")[0])]=delta_pred[i]
    arr_a[arr_a == 0] = None
    arr_c[arr_c == 0] = None
    arr_g[arr_g == 0] = None
    arr_t[arr_t == 0] = None
    #ax = fig.add_subplot(ntrack, 1, track_no)
    ax.set_ylabel('In silico\nMutagenesis\n')
    ax.set_title(title)
    ax.scatter(range(seq_shape[0]), -1 * arr_a, label='A', color='green', s = dotsize)
    ax.scatter(range(seq_shape[0]), -1 * arr_c, label='C', color='blue', s = dotsize)
    ax.scatter(range(seq_shape[0]), -1 * arr_g, label='G', color='orange', s = dotsize)
    ax.scatter(range(seq_shape[0]), -1 * arr_t, label='T', color='red', s = dotsize)
    ax.legend()
    ax.axhline(y=0, linestyle='--', color='gray')
    ax.set_xlim((0, seq_shape[0]))
    _ = ax.set_xticks(np.arange(0, seq_shape[0]+1, 10))

    return ax


fig, axs = plt.subplots(nrows = 10, figsize=(50, 50))
for ax, region in zip(axs, eRegulon_metadata.query('Region_signature_name == "HNF4A_extended_+_+_(4305r)"').nsmallest(10, 'triplet_score')['Region']):
	print(region)
	plot_deepexplainer_givenax_givenregion(region, track_nr, ax)
	ax.set_title(region)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'top_10_triplet_HNF4A.pdf'))
	
track_nr = [x for x in descriptions_of_interest if x[1].split(':')[1] == 'FOXA2' and x[1].split(':')[2] == 'HepG2'][0][0]
fig, axs = plt.subplots(nrows = 10, figsize=(50, 50))
for ax, region in zip(axs, eRegulon_metadata.query('Region_signature_name == "FOXA2_extended_+_+_(2070r)"').nsmallest(10, 'triplet_score')['Region']):
	print(region)
	plot_deepexplainer_givenax_givenregion(region, track_nr, ax)
	ax.set_title(region)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'top_10_triplet_FOXA2.pdf'))


selected_region_HEPG2 = eRegulon_metadata.query('TF == "HNF4A"').sort_values('triplet_score')['Region'].tolist()[0]
TFs_targettting_selected_region = set(eRegulon_metadata.query('Region == @selected_region_HEPG2')['TF'])
target_genes_of_selected_region = set(eRegulon_metadata.query('Region == @selected_region_HEPG2')['Gene'])
track_nr = [x for x in descriptions_of_interest if x[1].split(':')[1] == 'HNF4A' and x[1].split(':')[2] == 'HepG2'][0][0]

original_pred, x = calculate_saturation_mutagenesis_for_region(selected_region_HEPG2, track_nr)

fig, axs = plt.subplots(nrows = 2, figsize = (17.25, 4.75), sharex = True)
plot_deepexplainer_givenax_givenregion(selected_region_HEPG2, track_nr, axs[0])
plot_mutagenesis_givenax_fast(x, original_pred, track_nr, axs[1], dotsize = 4)
axs[0].set_xticks(list(range(0, 600, 100)))
sns.despine(ax = axs[0])
sns.despine(ax = axs[1])
fig.tight_layout()
fig.savefig(os.path.join(out_dir, f'{selected_region_HEPG2}_enformer_de.pdf'))
fig.savefig(os.path.join(out_dir, f'{selected_region_HEPG2}_enformer_de.png'))

track_nr = [x for x in descriptions_of_interest if x[1].split(':')[1] == 'FOXA2' and x[1].split(':')[2] == 'HepG2'][0][0]

original_pred, x = calculate_saturation_mutagenesis_for_region(selected_region_HEPG2, track_nr)

fig, axs = plt.subplots(nrows = 2, figsize = (17.25, 4.75), sharex = True)
plot_deepexplainer_givenax_givenregion(selected_region_HEPG2, track_nr, axs[0])
plot_mutagenesis_givenax_fast(x, original_pred, track_nr, axs[1], dotsize = 4)
axs[0].set_xticks(list(range(0, 600, 100)))
sns.despine(ax = axs[0])
sns.despine(ax = axs[1])
fig.tight_layout()
fig.savefig(os.path.join(out_dir, f'{selected_region_HEPG2}_enformer_de_FOXA2.pdf'))
fig.savefig(os.path.join(out_dir, f'{selected_region_HEPG2}_enformer_de_FOXA2.png'))

```

## Main Figure 4, pannel d: network

```python


import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import quilter
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

out_dir = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/figure_target_regions_DPCL'

#Load scplus_obj
import dill
f_scplus_obj = '/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/DPCL/scenicplus_final_autoreg/grnboost/scplus_obj.pkl'
scplus_obj = dill.load(open(f_scplus_obj, 'rb'))

#calculate triplet score and TF-to-region score
from mudata import AnnData
import mudata

adata_max_rank = mudata.read('/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/triplet_score/DPCL_adata_max_rank.h5ad')
df_max_rank = adata_max_rank.to_df()

eRegulon_metadata_key = 'eRegulon_metadata'
TF_region_iter = scplus_obj.uns[eRegulon_metadata_key][['TF', 'Region']].to_numpy()
TF_to_region_score = [df_max_rank.loc[region, TF] for TF, region in TF_region_iter]
scplus_obj.uns[eRegulon_metadata_key]['TF_to_region_max_rank'] = TF_to_region_score

from scenicplus.triplet_score import _rank_scores_and_assign_random_ranking_in_range_for_ties, _calculate_cross_species_rank_ratio_with_order_statistics
TF2G_score_key = 'TF2G_importance'
R2G_score_key = 'R2G_importance'
TF2R_score_key = 'TF_to_region_max_rank'
TF2G_score = scplus_obj.uns[eRegulon_metadata_key][TF2G_score_key].to_numpy()
R2G_score = scplus_obj.uns[eRegulon_metadata_key][R2G_score_key].to_numpy()
TF2R_score = scplus_obj.uns[eRegulon_metadata_key][TF2R_score_key].to_numpy()
#rank the scores
TF2G_rank = _rank_scores_and_assign_random_ranking_in_range_for_ties(TF2G_score)
R2G_rank = _rank_scores_and_assign_random_ranking_in_range_for_ties(R2G_score)
TF2R_rank = _rank_scores_and_assign_random_ranking_in_range_for_ties(-TF2R_score) #negate because lower score is better

import numpy as np
TF2G_rank_ratio = (TF2G_rank.astype(np.float64) + 1) / TF2G_rank.shape[0]
R2G_rank_ratio = (R2G_rank.astype(np.float64) + 1) / R2G_rank.shape[0]
TF2R_rank_ratio = (TF2R_rank.astype(np.float64) + 1) / TF2R_rank.shape[0]

rank_ratios = np.array([TF2G_rank_ratio, R2G_rank_ratio, TF2R_rank_ratio])
aggregated_rank = np.zeros((rank_ratios.shape[1],), dtype = np.float64)
for i in range(rank_ratios.shape[1]):
    aggregated_rank[i] = _calculate_cross_species_rank_ratio_with_order_statistics(rank_ratios[:, i])

scplus_obj.uns[eRegulon_metadata_key]['triplet_score'] = aggregated_rank.argsort().argsort()


TFs_of_interest = ['HNF4A', 'FOXA2', 'CEBPB']

columns_of_interest = ['TF', 'Region', 'Gene', 'triplet_score', 'TF_to_region_max_rank', 'R2G_importance']

tmp = scplus_obj.uns['eRegulon_metadata'].query("TF in @TFs_of_interest & TF2G_rho > 0")[columns_of_interest].nsmallest(10, 'triplet_score')
extra_data = []
for region in set(tmp['Region']):
	extra_data.append(scplus_obj.uns['eRegulon_metadata'].query("TF in @TFs_of_interest & TF2G_rho > 0").query('Region == @region'))
for gene in set(tmp['Gene']):
	extra_data.append(scplus_obj.uns['eRegulon_metadata'].query("TF in @TFs_of_interest & TF2G_rho > 0").query('Gene == @gene'))
data_of_interest = pd.concat([*extra_data, tmp])[columns_of_interest].drop_duplicates()



from scenicplus.diff_features import get_differential_features
get_differential_features(scplus_obj, 'ACC_Cell_type', use_hvg = False, contrast_type = ['DEGs', 'DARs'], log2fc_thr = np.NINF, adjpval_thr = np.inf)

data_of_interest['region_FC'] = [scplus_obj.uns['DARs']['ACC_Cell_type']['HepG2'].loc[
	region, 'Log2FC'] if region in scplus_obj.uns['DARs']['ACC_Cell_type']['HepG2'].index else np.nan for region in data_of_interest['Region']]
data_of_interest['gene_FC'] = [scplus_obj.uns['DEGs']['ACC_Cell_type']['HepG2'].loc[
	gene, 'Log2FC'] if gene in scplus_obj.uns['DEGs']['ACC_Cell_type']['HepG2'].index else np.nan for gene in data_of_interest['Gene']]

data_of_interest['TF_to_region_max_rank'] = (data_of_interest['TF_to_region_max_rank'].max() - data_of_interest['TF_to_region_max_rank']) / data_of_interest['TF_to_region_max_rank'].max()
data_of_interest['R2G_importance'] = (data_of_interest['R2G_importance'] - data_of_interest['R2G_importance'].min()) / (data_of_interest['R2G_importance'].max() - data_of_interest['R2G_importance'].min())

data_of_interest['triplet_score_inverted'] = abs( 
	(data_of_interest['triplet_score'] - data_of_interest['triplet_score'].max()) / (data_of_interest['triplet_score'].max() - data_of_interest['triplet_score'].min()) )
data_of_interest['TF_to_region_max_rank_invered'] = abs(
	(data_of_interest['TF_to_region_max_rank'] - data_of_interest['TF_to_region_max_rank'].max()) / (data_of_interest['TF_to_region_max_rank'].max() - data_of_interest['TF_to_region_max_rank'].min()))

scplus_obj.uns['metadata_for_network'] = data_of_interest

from scenicplus.networks import *

nx_tables = create_nx_tables(scplus_obj, eRegulon_metadata_key = 'metadata_for_network')

nx_tables['Node']['Gene'] = nx_tables['Node']['Gene'].merge(data_of_interest[['Gene', 'gene_FC']], on = 'Gene')
nx_tables['Node']['Gene'].index = nx_tables['Node']['Gene']['Gene'].to_numpy()
nx_tables['Node']['Gene'] = nx_tables['Node']['Gene'].drop_duplicates()

nx_tables['Node']['Region'] = nx_tables['Node']['Region'].merge(data_of_interest[['Region', 'region_FC']], on = 'Region')
nx_tables['Node']['Region'].index = nx_tables['Node']['Region']['Region'].to_numpy()
nx_tables['Node']['Region'] = nx_tables['Node']['Region'].drop_duplicates()


TF_to_color = {
	'HNF4A': '#97CC04', 
	'FOXA2': '#820263', 
	'CEBPB': '#41D3BD'}

G_kk, pos_kk, edge_tables_kk, node_tables_kk = create_nx_graph(
	nx_tables,
	color_edge_by = {
		'TF2R': {'variable': 'TF_to_region_max_rank_invered', 'continuous_color': 'binary', 'v_min': 0, 'v_max': 0.5},
		'R2G': {'variable': 'R2G_importance', 'continuous_color': 'binary', 'v_min': 0, 'v_max': 0.2}},
	color_node_by = {'TF': {'variable': 'TF', 'category_color' : TF_to_color},
			 'Gene': {'variable': 'gene_FC', 'continuous_color': 'Reds', 'v_max': 4, 'v_min': 0},
			 'Region': {'variable': 'region_FC', 'continuous_color': 'Blues', 'v_max': 4}},
	width_edge_by = {
		'TF2R': {'variable': 'TF_to_region_max_rank_invered', 'max_size': 2, 'min_size': 0.5},
		'R2G':  {'variable': 'R2G_importance', 'max_size': 2, 'min_size': 0.5}},
	size_node_by = {'TF': {'variable': 'fixed_size', 'fixed_size': 120},
                        'Gene': {'variable': 'fixed_size', 'fixed_size': 60},
                        'Region': {'variable': 'fixed_size', 'fixed_size': 40}},
	shape_node_by = {'TF': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
			 'Gene': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
			 'Region': {'variable': 'fixed_shape', 'fixed_shape': 'diamond'}},
        label_size_by = {'TF': {'variable': 'fixed_label_size', 'fixed_label_size': 20.0},
			 'Gene': {'variable': 'fixed_label_size', 'fixed_label_size': 10.0},
			 'Region': {'variable': 'fixed_label_size', 'fixed_label_size': 0.0}},
	layout='concentrical_layout', scale_position_by = 300)
from adjustText import adjust_text
from matplotlib.path import Path
from matplotlib.patches import FancyArrowPatch
from sklearn.cluster import AgglomerativeClustering
from scenicplus.utils import Groupby
from scenicplus.plotting.grn_plot import _line_two_points, _distance
d = 0.05
pos_kk = concentrical_layout(G_kk, dist_genes = 1, dist_TF = 0.1)
df_pos = pd.DataFrame(pos_kk, index = ['x', 'y']).T
for i in range(4):
	min_dist = 0.1
	pos_too_close = np.array(np.where(pairwise_distances(df_pos.to_numpy()) < min_dist)).T
	pos_too_close = pos_too_close[np.where(np.diff(pos_too_close) != 0 )[0], :]
	pos_too_close.sort()
	pos_too_close = np.unique(pos_too_close, axis = 0)
	for node_idx in pos_too_close[:, 0]:
		df_pos.iloc[node_idx]['x'] = df_pos.iloc[node_idx]['x'] + min_dist * 2
		df_pos.iloc[node_idx]['y'] = df_pos.iloc[node_idx]['y'] + min_dist * 2
df_pos.loc['SERPINF1']['y'] = df_pos.loc['SERPINF1']['y'] + 0.5
df_pos.loc['WDR81']['y'] = df_pos.loc['WDR81']['y'] + 0.5
df_pos.loc['SLC7A1']['x'] = df_pos.loc['SLC7A1']['x'] - 0.5
df_pos.loc['SLC7A1']['y'] = df_pos.loc['SLC7A1']['y'] + 0.5
from adjustText import adjust_text
fig = plt.figure(figsize = (4.75, 4.75), frameon = False)
ax = fig.add_axes([0, 0, 1, 1])
ax.axis('off')
X_regions = df_pos.loc[[x for x in df_pos.index if x.startswith('chr')], 'x'].to_numpy()
Y_regions = df_pos.loc[[x for x in df_pos.index if x.startswith('chr')], 'y'].to_numpy()
color_regions = node_tables_kk.loc[[x for x in df_pos.index if x.startswith('chr')], 'color'].to_numpy()
size_regions = node_tables_kk.loc[[x for x in df_pos.index if x.startswith('chr')], 'size'].to_numpy().astype(float)
ax.scatter(
	x = X_regions, y = Y_regions, color = color_regions, s = size_regions, facecolors = 'none', lw = 2, zorder = 3)
X_genes = df_pos.loc[[x for x in df_pos.index if not x.startswith('chr')], 'x'].to_numpy()
Y_genes = df_pos.loc[[x for x in df_pos.index if not x.startswith('chr')], 'y'].to_numpy()
color_genes = node_tables_kk.loc[[x for x in df_pos.index if not x.startswith('chr')], 'color'].to_numpy()
size_genes = node_tables_kk.loc[[x for x in df_pos.index if not x.startswith('chr')], 'size'].to_numpy().astype(float)
ax.scatter(
	x = X_genes, y = Y_genes, color = color_genes, s = size_genes, zorder = 3)
for source, target, color, width in edge_tables_kk.loc[[not x.startswith('chr') for x in edge_tables_kk['target']]][['source', 'target', 'color', 'width']].to_numpy():
	path = Path([df_pos.loc[source].to_list(), df_pos.loc[target].to_list()])
	con = FancyArrowPatch(path=path, arrowstyle='->', lw = width, color = color)
	ax.add_artist(con)
for TF in node_tables_kk.query('group == "TF"')['label']:
	regions = list(set(edge_tables_kk.loc[edge_tables_kk['source'] == TF, 'target']))
	region_positions = df_pos.loc[regions].to_numpy()
	posTF = df_pos.loc[TF].to_numpy()
	clustering = AgglomerativeClustering(distance_threshold = 2, n_clusters = None).fit(region_positions)
	grouper = Groupby(clustering.labels_)
	for ind in grouper.indices:
		region_pos = region_positions[ind]
		p_mean = region_pos.mean(0)
		m, b = _line_two_points(posTF, p_mean, return_func=False)
		# get a new point which is on this virtual line and d units towards posTF
		p_new = [p_mean[0] - d * np.sqrt(1 / (1 + m**2)),
		     p_mean[1] - m * d * np.sqrt(1 / (1 + m**2))]
		if _distance(p_new, posTF) > _distance(p_mean, posTF):
			p_new = [p_mean[0] + d * np.sqrt(1 / (1 + m**2)),
				 p_mean[1] + m * d * np.sqrt(1 / (1 + m**2))]
		# draw a straight line from posTF to p_new and curved lines from p_new to the region positions
		for target_point, width in zip(region_pos, edge_tables_kk.loc[edge_tables_kk['source'] == TF, 'width'].to_numpy()):
			verts = [
			    tuple(posTF),  # start of straight line
			    tuple(p_new),  # end of straight line
			    tuple(p_mean),  # control point for bezier curve
			    tuple(target_point)  # end point
			]
			codes = [
			    Path.MOVETO,
			    Path.LINETO,
			    Path.CURVE3,
			    Path.CURVE3
			]
			path = Path(verts, codes)
			con = FancyArrowPatch(path=path, arrowstyle='->', color = node_tables_kk.loc[TF, 'color'], lw = width)
			ax.add_artist(con)	
to_label = [node for node in df_pos.index if not node.startswith('chr')]
texts = [ax.text(df_pos.loc[node, 'x'], df_pos.loc[node, 'y'], node, fontweight = 'bold') for node in to_label]
adjust_text(texts)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'network.pdf'))
fig.savefig(os.path.join(out_dir, 'network.png'))

```

## Main Figure 4, pannel e: Genome browser plot

```python

import pyBigWig
def plot_coverage_bw(ax, bw_f, color, region, ymax, lw = 2, binsize = 10):
	chrom, start, end = region.replace(':', '-').split('-')
	start = int(start)
	end = int(end)
	x = np.array(range(start, end, binsize))
	bw = pyBigWig.open(bw_f)
	y = bw.stats(chrom, start, end, nBins = int((end - start) / binsize))
	y = [v if v is not None else 0 for v in y]
	y = np.nan_to_num(y)
	ax.fill_between(x[0:len(y)], y1=y, y2=0, step="mid",  linewidth=lw, color=color, edgecolor = 'black')
	ax.patch.set_alpha(0)
	ax.hlines(y = 0, xmin = start, xmax = end, colors = color, lw = lw)
	sns.despine(ax = ax, top = True, left = True, bottom = True, right = True)
	ax.set_ylim((0, ymax))
	

path_to_bw = '/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/DPCL/pycisTopic/consensus_peak_calling/pseudobulk_bw_files'
line_to_bw = {_file.split('.')[0]: os.path.join(path_to_bw, _file) for _file in os.listdir(path_to_bw)}

line_to_bw['HNF4A'] = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/enformer_comparison/ENCFF080FZD.bigWig'
line_to_bw['FOXA2'] = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/enformer_comparison/ENCFF626IVY.bigWig'
line_to_bw['CEBPB'] = '/staging/leuven/stg_00002/lcb/sdewin/PhD/papers/SCENICPLUS_BravoDeWinter_etal_2022/revisions_nat_meth/enformer_comparison/ENCFF003HJB.bigWig'

line_to_color = {
	'HepG2': '#51ae4f', 
	'PC3': '#f687c1', 
	'IMR90': '#a454ad', 
	'HCT116': '#377fbb', 
	'GM12878': '#e01e1e', 
	'Panc1': '#aa5c2b', 
	'MCF7': '#c8c923', 
	'K562': '#fe7f00',
	'HNF4A': '#97CC04', 
	'FOXA2': '#820263', 
	'CEBPB': '#41D3BD'}

line_to_ymax = {
	'HepG2': 5,
	'PC3': 5,
	'IMR90': 5,
	'HCT116': 5,
	'GM12878': 5,
	'Panc1': 5,
	'MCF7': 5,
	'K562': 5,
	'HNF4A': 40,
	'FOXA2': 40,
	'CEBPB': 40}

region_of_interest = 'chr4:87973864-88110972'
import matplotlib.patches as mpatches
import pyranges as pr
r2g = pd.read_csv('/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/DPCL/scenicplus_final_autoreg/grnboost/r2g.rho.bed', skiprows = 1, sep = '\t', header = None)
r2g.columns = ['Chromosome', 'Start', 'End', 'name', 'score', 'value', 'exp', 'color', 'sourceChrom', 'sourceStart', 'sourceEnd', 'sourceName', 'sourceStrand', 'targetChrom', 'targetStart', 'targetEnd', 'targetName', 'targetStrand']
r2g = pr.PyRanges(r2g)

chrom, start, end = region_of_interest.replace(':', '-').split('-')
pr_region_of_interest = pr.PyRanges(pd.DataFrame({'Chromosome': chrom, 'Start': int(start), 'End': int(end)}, index = [1]))

r2g_region = r2g.intersect(pr_region_of_interest).df
r2g_region.loc[r2g_region['name'] == 'SPP1_1']

fig, axs = plt.subplots(figsize = (7.25, 4.75), nrows = len(line_to_color.keys()) + 2, frameon = False, sharex = True, sharey = False,
	gridspec_kw = dict(hspace = 0))
for i, line in enumerate(line_to_color.keys()):
	print(line)
	plot_coverage_bw(axs[i], line_to_bw[line], line_to_color[line], region_of_interest, lw = 0.5, binsize = 100, ymax = line_to_ymax[line])
	axs[i].get_xaxis().set_visible(False)
x = np.array(range(int(start), int(end), 100))
for _, r2g in r2g_region.loc[r2g_region['name'] == 'SPP1_1'].sort_values('value').iterrows():
	posA = (int(r2g['sourceStart']), 0)
	posB = (int(r2g['targetStart']), 0)
	# this to ensure arcs are always down
	sign = '-' if posA[0] > posB[0] else ''
	color = r2g['color']
	color = [int(x) / 255 for x in color.split(',')]
	if (posA[0] > x.min() and posA[0] < x.max()) and (posB[0] > x.min() and posB[0] < x.max()):
		print('line')
		arrow = mpatches.FancyArrowPatch(posA=posA,
						 posB=posB,
						 connectionstyle=f"arc3,rad={sign}0.08",
						 color=color,
						 lw=1)
		axs[i+1].add_patch(arrow)
		axs[i+1].set_ylim((-1, 0))
		sns.despine(ax = axs[i+1], top = True, bottom = True, left = True, right = True)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'coverage.pdf'))
fig.savefig(os.path.join(out_dir, 'coverage.png'))

```

