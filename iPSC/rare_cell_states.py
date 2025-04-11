# Commented out IPython magic to ensure Python compatibility.
import palantir
import anndata as ad
import scanpy as sc
import pandas as pd
import os
import scipy.io

# Plotting
import matplotlib
import matplotlib.pyplot as plt

# warnings
import warnings
from numba.core.errors import NumbaDeprecationWarning

warnings.filterwarnings(action="ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings(
    action="ignore", module="scanpy", message="No data for colormapping"
)

# Inline plotting
# %matplotlib inline

matrix_file = "/content/d8/D8_matrix.mtx.gz"
genes_file = "/content/d8/D8_genes.tsv.gz"
barcodes_file = "/content/d8/D8_barcodes.tsv.gz"

expression_matrix = scipy.io.mmread(matrix_file).tocsc()
barcodes = pd.read_csv(barcodes_file, header=None, sep="\t")[0].values
genes = pd.read_csv(genes_file, header=None, sep="\t")[1].values

ad8 = sc.AnnData(X=expression_matrix.T)
ad8.var_names = genes
ad8.obs_names = barcodes

print(ad8)

matrix_file = "/content/d14/D14_matrix.mtx.gz"
genes_file = "/content/d14/D14_genes.tsv.gz"
barcodes_file = "/content/d14/D14_barcodes.tsv.gz"

expression_matrix = scipy.io.mmread(matrix_file).tocsc()
barcodes = pd.read_csv(barcodes_file, header=None, sep="\t")[0].values
genes = pd.read_csv(genes_file, header=None, sep="\t")[1].values

ad14 = sc.AnnData(X=expression_matrix.T)
ad14.var_names = genes
ad14.obs_names = barcodes

print(ad14)

matrix_file = "/content/d21/D21_matrix.mtx.gz"
genes_file = "/content/d21/D21_genes.tsv.gz"
barcodes_file = "/content/d21/D21_barcodes.tsv.gz"

expression_matrix = scipy.io.mmread(matrix_file).tocsc()
barcodes = pd.read_csv(barcodes_file, header=None, sep="\t")[0].values
genes = pd.read_csv(genes_file, header=None, sep="\t")[1].values

ad21 = sc.AnnData(X=expression_matrix.T)
ad21.var_names = genes
ad21.obs_names = barcodes

print(ad21)

matrix_file = "/content/d28/D28_matrix.mtx.gz"
genes_file = "/content/d28/D28_genes.tsv.gz"
barcodes_file = "/content/d28/D28_barcodes.tsv.gz"

expression_matrix = scipy.io.mmread(matrix_file).tocsc()
barcodes = pd.read_csv(barcodes_file, header=None, sep="\t")[0].values
genes = pd.read_csv(genes_file, header=None, sep="\t")[1].values

ad28 = sc.AnnData(X=expression_matrix.T)
ad28.var_names = genes
ad28.obs_names = barcodes

print(ad28)

matrix_file = "/content/d35/D35_matrix.mtx.gz"
genes_file = "/content/d35/D35_genes.tsv.gz"
barcodes_file = "/content/d35/D35_barcodes.tsv.gz"

expression_matrix = scipy.io.mmread(matrix_file).tocsc()
barcodes = pd.read_csv(barcodes_file, header=None, sep="\t")[0].values
genes = pd.read_csv(genes_file, header=None, sep="\t")[1].values

ad35 = sc.AnnData(X=expression_matrix.T)
ad35.var_names = genes
ad35.obs_names = barcodes

print(ad35)

ad8.obs_names = ["8_" + bc.replace("-1", "") for bc in ad8.obs_names]
ad14.obs_names = ["14_" + bc.replace("-1", "") for bc in ad14.obs_names]
ad21.obs_names = ["21_" + bc.replace("-1", "") for bc in ad21.obs_names]
ad28.obs_names = ["28_" + bc.replace("-1", "") for bc in ad28.obs_names]
ad35.obs_names = ["35_" + bc.replace("-1", "") for bc in ad35.obs_names]

adata_list = [ad8, ad14, ad21, ad28, ad35]
for ada in adata_list:
    ada.var_names_make_unique()

adata = ad.concat(adata_list, label="batch", keys=["D8", "D14", "D21", "D28", "D35"])
print(adata)

dnlin = pd.read_csv("DN_linage_barcodes.csv", encoding="utf-8")

dnlin.shape

dnlin.head()

dnlin["barcode_processed"] = dnlin["X"].str.replace("^day", "", regex=True)
common_barcodes = set(dnlin["barcode_processed"]).intersection(set(adata.obs_names))
dnlin_sub = dnlin[dnlin["barcode_processed"].isin(common_barcodes)]
adata = adata[adata.obs_names.isin(common_barcodes)].copy()

adata.obs["G2M_Score"] = dnlin_sub.set_index("barcode_processed").loc[adata.obs_names, "G2M_Score"].values
adata.obs["S_Score"] = dnlin_sub.set_index("barcode_processed").loc[adata.obs_names, "S_Score"].values
adata.obs["Phase"] = dnlin_sub.set_index("barcode_processed").loc[adata.obs_names, "Phase"].values
adata.obs["Day"] = dnlin_sub.set_index("barcode_processed").loc[adata.obs_names, "day"].values
adata.obs["cell_type"] = dnlin_sub.set_index("barcode_processed").loc[adata.obs_names, "cell_type"].values

adata

adata.obs

adata.X = adata.X.astype("float32")
sc.pp.normalize_per_cell(adata)

palantir.preprocess.log_transform(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=1500, flavor="seurat")

sc.pp.pca(adata)
adata

dm_res = palantir.utils.run_diffusion_maps(adata, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(adata)

sc.pp.neighbors(adata)
sc.tl.umap(adata)
adata

sc.pl.embedding(
    adata,
    basis="umap",
    frameon=False,
)

imputed_X = palantir.utils.run_magic_imputation(adata)

import numpy as np

# تبدیل Index به یک لیست و جایگزینی مقدار مورد نظر
adata.var_names = np.where(adata.var_names == "CASC10", "MIR1915HG", adata.var_names)

sc.pl.embedding(
    adata,
    basis="umap",
    layer="MAGIC_imputed_data",
    color=["MIAT", "LINC00632", "KCNQ1OT1", "ZFAS1", "AL139246.5", "TP53TG1",
              "AL391807.1", "JPX", "MIR1915HG", "SHH", "FOXA2", "LMX1B",
              "LMX1A", "NHLH1", "PITX3", "TH", "NR4A2", "NHLH2", "UCHL1", "NES",
              "SYP", "NEUROD1", "SOX2", "NANOG"],
    ncols=4,
    frameon=False,
)
plt.show()

sc.pl.umap(
    adata,
    color="cell_type",
    # Setting a smaller point size to get prevent overlap
    size=10,
)

sc.pl.umap(
    adata,
    color="Phase",
    # Setting a smaller point size to get prevent overlap
    size=10,
)

palantir.plot.plot_diffusion_components(adata)

print(adata.obs_names)

max_umap1_idx = np.argmax(adata.obsm["X_umap"][:, 0])
max_umap2_idx = np.argmax(adata.obsm["X_umap"][:, 1])

max_umap1_cell = adata.obs_names[max_umap1_idx]
max_umap2_cell = adata.obs_names[max_umap2_idx]

print(f"max UMAP1 or end point: {max_umap1_cell}")
print(f"max UMAP2 or start point: {max_umap2_cell}")

start_cell = "14_GTGGGTCCACATGGGA"
terminal_states = pd.Series(
    ["DN"],
    index=["35_GCATGCGTCATGGTCA"]
)

palantir.plot.highlight_cells_on_umap(adata, [max_umap1_cell, max_umap2_cell])

pr_res = palantir.core.run_palantir(
    adata, start_cell, num_waypoints=500, terminal_states=terminal_states
)

palantir.plot.plot_palantir_results(adata, s=3)
plt.show()

masks = palantir.presults.select_branch_cells(adata, q=.01, eps=.01)

palantir.plot.plot_trajectory(adata, "DN")
plt.show()

gene_trends = palantir.presults.compute_gene_trends(
    adata,
    expression_key="MAGIC_imputed_data",
)

genes = ["MIAT", "LINC00632", "KCNQ1OT1", "ZFAS1", "AL139246.5", "TP53TG1",
              "AL391807.1", "JPX", "MIR1915HG", "SHH", "FOXA2", "LMX1B",
              "LMX1A", "NHLH1", "PITX3", "TH", "NR4A2", "NHLH2", "UCHL1", "NES",
              "SYP", "NEUROD1", "SOX2", "NANOG"]
palantir.plot.plot_gene_trends(adata, genes)
plt.show()

palantir.plot.plot_gene_trend_heatmaps(adata, genes)
plt.show()

palantir.plot.plot_trend(adata, "DN", "TP53TG1", color="palantir_entropy", position_layer="MAGIC_imputed_data")
palantir.plot.plot_trend(adata, "DN", "TP53TG1", color="palantir_pseudotime", position_layer="MAGIC_imputed_data")
palantir.plot.plot_trend(adata, "DN", "TP53TG1", color="n_counts", position_layer="MAGIC_imputed_data")
plt.show()

!pip install --q mellon

import mellon
model = mellon.DensityEstimator()
log_density = model.fit_predict(adata.obsm["DM_EigenVectors"])

predictor = model.predict

adata.obs["mellon_log_density"] = log_density
adata.obs["mellon_log_density_clipped"] = np.clip(
    log_density, *np.quantile(log_density, [0.05, 1])
)

adata

adata.obs

sc.pl.scatter(
    adata, color=["mellon_log_density", "mellon_log_density_clipped"], basis="umap"
)

fig, (ax1, ax2) = plt.subplots(1, 2, width_ratios=[3, 2], figsize=[12, 4])
sc.pl.violin(adata, "mellon_log_density", "cell_type", rotation=45, ax=ax1, show=False)
sc.pl.scatter(adata, color="cell_type", basis="umap", ax=ax2, show=False)
plt.show()

palantir.plot.plot_branch(
    adata,
    branch_name="DN",
    position="mellon_log_density",
    color="cell_type",
    masks_key="branch_masks",
    s=100,
)
plt.show()
palantir.plot.plot_branch(
    adata,
    branch_name="DN",
    position="mellon_log_density_clipped",
    color="cell_type",
    masks_key="branch_masks",
    s=100,
)
plt.show()

ct_colors = pd.Series(
    adata.uns["cell_type_colors"], index=adata.obs["cell_type"].values.categories
)

plt.figure(figsize=[8, 3])
plt.scatter(
    adata.obs["mellon_log_density_clipped"],
    adata.obs["mellon_log_density"],
    s=10,
    color=ct_colors[adata.obs["cell_type"]],
)
plt.show()

palantir.plot.plot_branch(
    adata,
    branch_name="DN",
    position="mellon_log_density",
    color= "TP53TG1",
    color_layer="MAGIC_imputed_data",
    masks_key="branch_masks",
    s=100,
)
plt.show()

palantir.plot.plot_trend(
    adata,
    "DN",
    "TP53TG1",
    masks_key="branch_masks",
    color="TP53TG1",
    position="mellon_log_density",
    color_layer="MAGIC_imputed_data",
)
plt.show()

palantir.utils.run_local_variability(adata)

palantir.plot.plot_trend(
    adata,
    "DN",
    "TP53TG1",
    masks_key="branch_masks",
    color="TP53TG1",
    position="mellon_log_density",
    color_layer="local_variability",
)
plt.show()

for gene in genes:
    palantir.plot.plot_trend(
        adata,
        "DN",
        gene,
        masks_key="branch_masks",
        color=gene,
        position="mellon_log_density",
        color_layer="MAGIC_imputed_data",
    )
    plt.show()

for gene in genes:
    palantir.plot.plot_trend(
        adata,
        "DN",
        gene,
        masks_key="branch_masks",
        color=gene,
        position="mellon_log_density",
        color_layer="local_variability",
        s=100,
    cmap="Purples",
    edgecolor="black",
    linewidth=0.1,
    )
    plt.show()

# 1. Re-computing density with low intrinsic dimensionality to make sure weights are representative
model = mellon.DensityEstimator(d_method="fractal")
log_density = model.fit_predict(adata.obsm["DM_EigenVectors"])

adata.obs["mellon_log_density_lowd"] = log_density

# 2. Computing scores for each gene
score_key = "change_scores"
palantir.utils.run_low_density_variability(
    adata,
    cell_mask="branch_masks",
    density_key="mellon_log_density_lowd",
    score_key=score_key,
)
plt.show()

scores = adata.var["change_scores_DN"]
scores.sort_values(ascending=False)

# The scores can be visualized as a histogram. Some important genes are highlgited below
highlight_gene = ["MIAT", "AL139246.5", "TP53TG1", "SOX2", "NES"]

palantir.plot.gene_score_histogram(adata, "change_scores_DN", highlight_gene)
plt.show()

# The scores can be visualized as a histogram. Some important genes are highlgited below
palantir.plot.gene_score_histogram(adata, "change_scores_DN", genes)
plt.show()

selected_trends = adata.varm["gene_trends_DN"].loc[genes, :]

selected_trends

import seaborn as sns
sns.clustermap(
    selected_trends,
    z_score=0,
    cmap="RdBu_r",
    col_cluster=False,
    vmin=-2,
    vmax=2,
    figsize=[7, 7],
    xticklabels=0,
)
plt.show()

# Alternatively sort by maximum exprression bin for clean visualization
plot_order = selected_trends.idxmax(axis=1).sort_values().index
sns.clustermap(
    selected_trends.loc[plot_order, :],
    z_score=0,
    cmap="RdBu_r",
    col_cluster=False,
    row_cluster=False,
    vmin=-2,
    vmax=2,
    figsize=[7, 7],
    xticklabels=0,
)
plt.show()

percentile = 0.95
percentile_cutoff = np.quantile(scores, percentile)
scores = scores.sort_values(ascending=False)
dopa_change_genes = scores.index[scores >= percentile_cutoff]
selected_trends = adata.varm["gene_trends_DN"].loc[dopa_change_genes, :]
selected_trends.shape

# Visualize the expression and variability of top 5 genes
sc.pl.embedding(
    adata,
    basis="umap",
    color=dopa_change_genes[:5],
    layer="MAGIC_imputed_data",
    ncols=5,
)
sc.pl.embedding(
    adata,
    basis="umap",
    color=dopa_change_genes[:5],
    layer="local_variability",
    ncols=5,
)
plt.show()

sns.clustermap(
    selected_trends,
    z_score=0,
    cmap="RdBu_r",
    col_cluster=False,
    vmin=-2,
    vmax=2,
    figsize=[7, 7],
    xticklabels=0,
)
plt.show()

# Alternatively sort by maximum exprression bin for clean visualization
plot_order = selected_trends.idxmax(axis=1).sort_values().index
sns.clustermap(
    selected_trends.loc[plot_order, :],
    z_score=0,
    cmap="RdBu_r",
    col_cluster=False,
    row_cluster=False,
    vmin=-2,
    vmax=2,
    figsize=[7, 12],
    xticklabels=0,
)
plt.show()

adata.obs

adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

adata

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
miat_expression = adata[:, "MIAT"].X.toarray().flatten() if hasattr(adata[:, "MIAT"].X, "toarray") else adata[:, "MIAT"].X.flatten()
# استخراج مقادیر G2M_Score
g2m_score = adata.obs["G2M_Score"].values
correlation = np.corrcoef(miat_expression, g2m_score)[0, 1]
correlation, p_value = pearsonr(miat_expression, g2m_score)
df = pd.DataFrame({"MIAT Expression": miat_expression, "G2M Score": g2m_score})
plt.figure(figsize=(8, 6))
sns.scatterplot(data=df, x="MIAT Expression", y="G2M Score", alpha=0.5)
sns.regplot(data=df, x="MIAT Expression", y="G2M Score", scatter=False, color="red")
plt.title(f"Correlation between MIAT Expression and G2M Score\n(Pearson r: {correlation:.2f}, p-value: {p_value:.2e})")
plt.xlabel("MIAT Expression")
plt.ylabel("G2M Score")
plt.show()

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
miat_expression = adata[:, "MIAT"].X.toarray().flatten() if hasattr(adata[:, "MIAT"].X, "toarray") else adata[:, "MIAT"].X.flatten()
# استخراج مقادیر S_Score
g2m_score = adata.obs["S_Score"].values
correlation = np.corrcoef(miat_expression, g2m_score)[0, 1]
correlation, p_value = pearsonr(miat_expression, g2m_score)
df = pd.DataFrame({"MIAT Expression": miat_expression, "G2M Score": g2m_score})
plt.figure(figsize=(8, 6))
sns.scatterplot(data=df, x="MIAT Expression", y="G2M Score", alpha=0.5)
sns.regplot(data=df, x="MIAT Expression", y="G2M Score", scatter=False, color="red")
plt.title(f"Correlation between MIAT Expression and G2M Score\n(Pearson r: {correlation:.2f}, p-value: {p_value:.2e})")
plt.xlabel("MIAT Expression")
plt.ylabel("S Score")
plt.show()

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
miat_expression = adata[:, "TP53TG1"].X.toarray().flatten() if hasattr(adata[:, "MIAT"].X, "toarray") else adata[:, "MIAT"].X.flatten()
# استخراج مقادیر G2M_Score
g2m_score = adata.obs["G2M_Score"].values
correlation = np.corrcoef(miat_expression, g2m_score)[0, 1]
correlation, p_value = pearsonr(miat_expression, g2m_score)
df = pd.DataFrame({"TP53TG1 Expression": miat_expression, "G2M Score": g2m_score})
plt.figure(figsize=(8, 6))
sns.scatterplot(data=df, x="TP53TG1 Expression", y="G2M Score", alpha=0.5)
sns.regplot(data=df, x="TP53TG1 Expression", y="G2M Score", scatter=False, color="red")
plt.title(f"Pearson r: {correlation:.2f}, p-value: {p_value:.2e})")
plt.xlabel("TP53TG1 Expression")
plt.ylabel("G2M Score")
plt.show()

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
miat_expression = adata[:, "TP53TG1"].X.toarray().flatten() if hasattr(adata[:, "MIAT"].X, "toarray") else adata[:, "MIAT"].X.flatten()
# استخراج مقادیر G2M_Score
g2m_score = adata.obs["S_Score"].values
correlation = np.corrcoef(miat_expression, g2m_score)[0, 1]
correlation, p_value = pearsonr(miat_expression, g2m_score)
df = pd.DataFrame({"TP53TG1 Expression": miat_expression, "S Score": g2m_score})
plt.figure(figsize=(8, 6))
sns.scatterplot(data=df, x="TP53TG1 Expression", y="S Score", alpha=0.5)
sns.regplot(data=df, x="TP53TG1 Expression", y="S Score", scatter=False, color="red")
plt.title(f"Pearson r: {correlation:.2f}, p-value: {p_value:.2e})")
plt.xlabel("TP53TG1 Expression")
plt.ylabel("S Score")
plt.show()
