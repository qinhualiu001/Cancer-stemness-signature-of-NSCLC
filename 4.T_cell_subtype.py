import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import scanpy as sc
import os

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300, facecolor="white")

adata = sc.read("GSE148071_cell_type.h5")
adata

sub_adata = adata[adata.obs['cell_type'].isin(["T"])]

#sub_adata = immune_sub_adata.raw.to_adata()
#sub_adata.raw = sub_adata
#sc.pp.highly_variable_genes(sub_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sub_adata = sub_adata[:, sub_adata.var.highly_variable]
#sc.pp.scale(sub_adata, max_value=10)

sc.tl.pca(sub_adata, svd_solver="arpack")
sc.pp.neighbors(sub_adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(sub_adata)
sc.tl.leiden(sub_adata, resolution=1)
sc.pl.umap(sub_adata, color=["leiden"],save="_T_leiden_resolution1.png",legend_loc="on data")

sc.tl.rank_genes_groups(sub_adata, "leiden", method="wilcoxon")
#save the DEG analysis results and mannual annotate the cell type with CellMarker2 database
result = sub_adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names
DEG_res = pd.DataFrame(
    {
        group + "_" + key[:1]: result[key][group]
        for group in groups
        for key in ["names", "pvals","scores","pvals_adj","logfoldchanges"]
    }
)
# 定义每个工作表的列数
columns_per_sheet = 4
# 创建 ExcelWriter 对象
writer = pd.ExcelWriter('T_leiden_DEG_res.xlsx')
# 将数据写入多个工作表
for i in range(0, len(DEG_res.columns), columns_per_sheet):
    print(i)
    sheet_name = f'Sheet{i // columns_per_sheet + 1}'  # 工作表名称
    start_col = i  # 起始列索引
    end_col = i + columns_per_sheet  # 结束列索引
    # 提取要保存的列
    subset = DEG_res.iloc[:, start_col:end_col]
    # 写入工作表
    subset.to_excel(writer, sheet_name=sheet_name, index=False)

# 保存 Excel 文件
writer.close()


###### cell type annotation
cell_type_annotation = {'0':'Naive T cell','1':'CD8+ effector T cell',
'2':'CD4+ T cell','3':'CD8+ effector T cell','4':'CD8+ effector T cell',
	'5':'CD4+ T cell','6':'CD8+ effector T cell','7':'Proliferative T cell','8':'Epithelial cell',
	'9':'CD8+ effector T cell','10':'Epithelial cell','11':'Naive CD4+ T cell'}
sub_adata.obs['cell_type'] = sub_adata.obs['leiden'].map(cell_type_annotation)

sub_adata = sub_adata[sub_adata.obs["cell_type"]!="Epithelial cell"]
sc.pl.umap(sub_adata, color=["cell_type"],save="_cell_type_T.png",title="Cell type")

sub_adata.write("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071_T_cell.h5")

adata = sc.read("GSE148071_T_cell.h5")
adata
sc.tl.rank_genes_groups(adata, "cell_type", method="wilcoxon")

marker_genes_dict = {
    "CD4+ T cell": ['TXN','CTSL','SCPEP1','ANXA1','S100A11'],
    "CD8+ effector T cell": ['CD8A','GZMA','GZMB'],
    "Naive CD4+ T cell":["CD55"],
    "Naive T cell": ['IL7R',"LEF1"],
    "Proliferative T cell": ['MKI67','TOP2A']
}
sc.pl.heatmap(adata, marker_genes_dict, groupby="cell_type",show_gene_labels=True, cmap="viridis", dendrogram=False,save="_T_marker_heatmap.png")
sc.pl.stacked_violin(
    adata, marker_genes_dict, groupby="cell_type", swap_axes=False, dendrogram=True,save="cell_type_T.png"
)
sc.pl.dotplot(adata, marker_genes_dict, "cell_type", dendrogram=False,save="_dotplot_celltype_T.png")
