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

#immune_sub_adata = adata[adata.obs['cell_type'].isin(["T","B","Follicular B","Mast","Myeloid"])]
immune_sub_adata = adata[adata.obs['cell_type'].isin(["Myeloid"])]

sub_adata = immune_sub_adata.raw.to_adata()
sub_adata.raw = sub_adata
sc.pp.highly_variable_genes(sub_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sub_adata = sub_adata[:, sub_adata.var.highly_variable]
sc.pp.scale(sub_adata, max_value=10)
sc.tl.pca(sub_adata, svd_solver="arpack")

sc.pp.neighbors(sub_adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(sub_adata)
sc.tl.leiden(sub_adata, resolution=0.5)
sc.pl.umap(sub_adata, color=["leiden"],save="_Myeloid_leiden_resolution1.png",legend_loc="on data")

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
writer = pd.ExcelWriter('Myeloid_leiden_DEG_res.xlsx')
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


cell_type_annotation = {'0':'Macrophage cell','1':'M2 macrophage','2':'M1 macrophage','3':'Natural killer cell','4':'Neutrophil cell',
	'5':'Monocyte cell','6':'Monocyte cell','7':'Neutrophil cell','8':'Epithelial cell','9':'Epithelial cell',
	'10':'Red blood cell','11':'Natural killer cell','12':'Epithelial cell','13':'Epithelial cell','14':'Dendritic cell',
	'15':'Macrophage cell','16':'Epithelial cell','17':'Epithelial cell','18':'Epithelial cell','19':'Epithelial cell'}
sub_adata.obs['cell_type'] = sub_adata.obs['leiden'].map(cell_type_annotation)

sub_adata = sub_adata[sub_adata.obs["cell_type"]!="Epithelial cell"]
sc.pl.umap(sub_adata, color=["cell_type"],save="_cell_type_Myeloid.png",title="Cell type")
marker_genes_dict = {
    "Dendritic cell": ['TXN','CTSL','SCPEP1','ANXA1'],
    "M1 macrophage": ['CTSB','CCL4','CD86'],
    "M2 macrophage":["MMP12","MRC1"],
    "Macrophage cell": ['APOE','CTSB'],
    "Monocyte cell": ['FCN1','MNDA'],
    "Natural killer cell": ['FCGR3A'],
    "Neutrophil cell": ['IFITM2','CXCL8','CSF3R'],
    "Red blood cell": ['HBB','HBA1','HBA2'],
}
sc.pl.heatmap(sub_adata, marker_genes_dict, groupby="cell_type",show_gene_labels=True, cmap="viridis", dendrogram=False,save="_myeloid_marker_heatmap.png")
sc.pl.dotplot(sub_adata, marker_genes_dict, "cell_type", dendrogram=False,save="_dotplot_celltype_Myeloid.png")

sub_adata.write("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071_Myeloid_cell.h5")
