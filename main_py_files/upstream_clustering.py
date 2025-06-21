import sopa
import spatialdata
import scanpy as sc
import squidpy as sq
import seaborn as sns
import spatialdata as sd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


""" ####################################################### """
""" _____________  Introducing The Parameters ______________"""
""" ####################################################### """

Tissue = "Region4"
TableFolder = "REGION4_TABLES"

out_path = "/Volumes/ProstateCancerEvo_SpatialAnalysis/dbs/Ongoing/{Tissue}/{TableFolder}"
os.makedirs(out_path, exist_ok=True)

phen_adata_file = "/Volumes/ProstateCancerEvo_SpatialAnalysis/dbs/Ongoing/Region2/REGION2_TABLES/Region_2_Xenium_Phen_HE_Integrated.Xenium_Process_Table.h5ad"
xenium_p_file = "/Volumes/ProstateCancerEvo_SpatialAnalysis/dbs/Ongoing/Region2/REGION2_TABLES/Region_2_Xenium_Phen_HE_Integrated.Protein_PhenCycTable.V2.h5ad"

clust_dir = os.path.join(out_path, "CLUST_UMAP_QC_ANN")

os.makedirs(clust_dir, exist_ok=True)

adata_xenium_p = sc.read_h5ad(phen_adata_file)
adata_pheno = sc.read_h5ad(xenium_p_file)

top_genes = 2000
flavor_chosen = "seurat"

""" ####################################################### """
""" __________________________ QC __________________________"""
""" ####################################################### """


sc.pp.calculate_qc_metrics(adata_xenium_p, percent_top=(10, 20, 50, 150), inplace=True)

cprobes = (
    adata_xenium_p.obs["control_probe_counts"].sum() / adata_xenium_p.obs["total_counts"].sum() * 100
)
cwords = (
    adata_xenium_p.obs["control_codeword_counts"].sum() / adata_xenium_p.obs["total_counts"].sum() * 100
)


print(f"Negative DNA probe count % : {cprobes}")
print(f"Negative decoding count % : {cwords}")


fig, axs = plt.subplots(1, 3, figsize=(15, 3))

axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata_xenium_p.obs["total_counts"],
    kde=False,
    ax=axs[0],
)


axs[1].set_title("Area of segmented cells")
sns.histplot(
    adata_xenium_p.obs["cell_area"],
    kde=False,
    ax=axs[1],
)

axs[2].set_title("Nucleus ratio")
sns.histplot(
    adata_xenium_p.obs["nucleus_area"] / adata_xenium_p.obs["cell_area"],
    kde=False,
    ax=axs[2],
)

fig.tight_layout()
fig.savefig(os.path.join(clust_dir, "QC.png"), dpi=300, bbox_inches="tight")

""" ###################################################################################### """
""" ___________ PhenoCycler Z - Score Annotation in terms of Biomarkers on UMAP ___________"""
""" ###################################################################################### """

marker_cell_dict_cells = {
'1': 'CD8', 
'2': 'CD31', 
'3': 'CD11c', 
'4': 'CD34', 
'5': 'CD20', 
'6': 'CD4', 
'7': 'CD21', 
'8': 'Ki67', 
'9': 'CD68', 
'10': 'E-cadherin', 
'11': 'pten', 
'12': 'NF-H',  
'13': 'FOXp3', 
'14': 'CD45', 
'15': 'TFAM', 
'16': 'PDL-1', 
'17': 'SOX2', 
'18': 'Vimentin', 
'19': 'HIF1A', 
'20': 'CD44', 
'21': 'Keratin 8/18', 
'22': 'aSMA',  
'24': 'K14'
}

#Find the enriched cells by specific marker genes from PhenoCycler data
sopa.utils.higher_z_score(adata_pheno, marker_cell_dict_cells)

# Record as cell_type
adata_xenium_p.obs["cell_type"] = (
    adata_pheno
         .obs["cell_type"]
         .values
)

# now that obs has two columns: cell_id & cell_type
adata_xenium_p.obs["cell_type"]


""" ###################################################################################### """
""" _______________________ Pre - Filtering Of Less Variable Genes  _______________________"""
""" ###################################################################################### """



# 1. Normalize to counts per 10 000 spots per cell (CP10K)
sc.pp.normalize_total(adata_xenium_p, target_sum=1e4)

# 2. Log‐transform
sc.pp.log1p(adata_xenium_p)

# 3. Find highly variable genes (HVGs)
sc.pp.highly_variable_genes(
    adata_xenium_p,
    n_top_genes=top_genes,
    flavor=flavor_chosen
)
adata_xenium_p = adata_xenium_p[:, adata_xenium_p.var.highly_variable] ### NOT: DANININ KUYRUGU BURADA KOPUYOR
sc.pp.scale(adata_xenium_p, max_value=10) # Doesn't allow to dominate outlier genes

sc.tl.pca(adata_xenium_p, svd_solver = "argpack")
sc.pl.pca_variance_ratio(adata_xenium_p, log=True)

plt.tight_layout()
plt.savefig(os.path.join(clust_dir, "PCA_Variance_Ratio.png"), dpi=300, bbox_inches="tight")



# 2. LEIDEN + UMAP
sc.pp.neighbors(adata_xenium_p, n_neighbors=15, n_pcs=40) # includes elbow plot calculation
sc.tl.leiden(adata_xenium_p, resolution=1.0, key_added="leiden")
sc.tl.umap(adata_xenium_p)

""" ###################################################################################### """
""" _______________________ Pre - Filtering Of Less Variable Genes  _______________________"""
""" ###################################################################################### """

sc.pl.umap(
    adata_xenium_p,
    color="leiden",
    palette="tab20",
    title="UMAP embedding — Leiden clusters",
    legend_loc="on data",
    frameon=False
)

plt.tight_layout()
plt.savefig(os.path.join(clust_dir, "UMAP.png"), dpi=300, bbox_inches="tight")


""" ###################################################################################### """
""" _____________ PhenoCycler Biomarkers enriched UMAP done by Transcriptome  _____________"""
""" ###################################################################################### """


cat_order = adata_xenium_p.obs['cell_type'].cat.categories.tolist()

palette = sns.color_palette("tab20", n_colors=len(cat_order)).as_hex()
adata_xenium_p.obs['cell_type'] = adata_xenium_p.obs['cell_type'].astype('category')
cat_order = adata_xenium_p.obs['cell_type'].cat.categories.tolist()

colors = palette  # list of colors in the same order as cat_order

# 3) Plot UMAP into a new Figure/Axes
fig, ax = plt.subplots(figsize=(4,4), dpi=200)
sc.pl.umap(
    adata_xenium_p,
    color='cell_type',
    palette=colors,
    ax=ax,
    show=False,
    frameon=False
)

if ax.get_legend() is not None:
    ax.get_legend().remove()

handles = [
    mpatches.Patch(color=colors[i], label=cat_order[i])
    for i in range(len(cat_order))
]
leg = ax.legend(
    handles=handles,
    title='Cell Type',
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    fontsize=6,
    title_fontsize=8
)

ax.tick_params(labelsize=3)
ax.title.set_fontsize(5)

plt.tight_layout()
plt.savefig(os.path.join(clust_dir, "UMAP_enrichedPhenoCyclerBiomarkers.png"), dpi=300, bbox_inches="tight")


#### SPLITS THE CELLS IN AND OUT CLUSTER COM
# 5. Rank “Xenium spot counts” per Leiden group
sc.tl.rank_genes_groups(
    adata_xenium_p,
    groupby="leiden",
    method="wilcoxon",
    key_added="rank_leiden"
)

sc.pl.rank_genes_groups(
    adata_xenium_p,
    key="rank_leiden",
    n_genes=10,
    sharey=False
)

plt.savefig(os.path.join(clust_dir, "VarianceExplained_Genes_Sorted_Top10.png"), dpi=300, bbox_inches="tight")

""" ###################################################################################### """
""" _____________________________ Leiden Cluster Specific Groups  _________________________"""
""" ###################################################################################### """

# 1) Flatten all DE results into a DataFrame
df = sc.get.rank_genes_groups_df(
    adata_xenium_p,
    group=None,           # None = include every Leiden cluster
    key="rank_leiden"     # the key you used when running rank_genes_groups
)                     

# 2) For each cluster, take the first 7 genes
top7 = (
    df
      .groupby("group")["names"]
      .apply(lambda genes: genes.iloc[:7].tolist())
)

# 3) Convert to a  table
df_top7 = pd.DataFrame({
    f"cluster_{grp}": genes 
    for grp, genes in top7.items()
})

df_top7 = df_top7[sorted(df_top7.columns, key=lambda x: int(x.split("_")[1]))]

df_top7

# 6. Visualize top markers
sc.pl.rank_genes_groups_heatmap(
    adata_xenium_p,
    key="rank_leiden",
    groupby="leiden",
    n_genes=5,
    swap_axes=False,
    show_gene_labels=True,
    dendrogram=True
)


plt.tight_layout()
plt.savefig(os.path.join(clust_dir, "LeidenGroups_Heatmap_Top7.png"), dpi=300, bbox_inches="tight")

sc.pl.rank_genes_groups_dotplot(
    adata_xenium_p,
    key="rank_leiden",
    groupby="leiden",
    n_genes=5,
    standard_scale="var"
)

plt.tight_layout()
plt.savefig(os.path.join(clust_dir, "LeidenGroups_DotPlot_Top7.png"), dpi=300, bbox_inches="tight")

# Record the UMAP introduced anndata
adata_xenium_p.write_h5ad(f"/Volumes/ProstateCancerEvo_SpatialAnalysis/dbs/Ongoing/Region2/REGION2_TABLES/{Tissue}_Xenium_Phen_HE_Integrated.Xenium_Process_Table.V2.h5ad")


""" ###################################################################################### """
""" _______________________________ UMAP Build By PhenoCycler  ____________________________"""
""" ###################################################################################### """

sc.pp.normalize_total(adata_pheno)
sc.pp.log1p(adata_pheno)
sc.pp.pca(adata_pheno)
sc.pp.neighbors(adata_pheno)
sc.tl.umap(adata_pheno)
sc.tl.leiden(adata_pheno, resolution=0.1)


sc.pl.umap(
    adata_pheno,
    color='cell_type',
    palette="tab20",
    title='UMAP by z-score enriched protein markers',
    legend_loc='right margin',
    frameon=False
)

plt.tight_layout()
plt.savefig(os.path.join(clust_dir, "UMAP_doneByPhenoCycler.png"), dpi=300, bbox_inches="tight")