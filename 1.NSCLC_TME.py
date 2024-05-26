import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import scanpy as sc
import os

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300, facecolor="white")

# a = sc.read('/p300s/jiapl_group/kanghongen/projects/TISCH/anndata_h5ad/NSCLC_GSE176021_aPD1.h5ad')
# b = sc.read("/p300s/jiapl_group/kanghongen/projects/TISCH/anndata_h5ad/NSCLC_GSE151537.h5ad")
# c = sc.read("/p300s/jiapl_group/kanghongen/projects/TISCH/anndata_h5ad/NSCLC_GSE146100.h5ad")

# 定义目录路径
directory = '/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071'
# 创建一个空的数据框用于存储合并后的数据
merged_df = pd.DataFrame()
file_list = sorted(os.listdir(directory))
patient_list = []
# 遍历目录下的所有文件
patient_count = 1
for filename in file_list:
    print(filename)
    if filename.endswith('.txt'):
        # 构建文件的完整路径
        filepath = os.path.join(directory, filename)
        # 读取文件数据并使用第一列作为行名
        df = pd.read_csv(filepath, sep='\t', index_col=0)
        n = df.shape[1]
        strings = ['patient_' + str(patient_count)] * n
        patient_list = patient_list + strings
        # 合并数据框
        merged_df = pd.concat([merged_df, df], axis=1)
        patient_count = patient_count + 1

merged_df.shape
#(29527, 89887)
len(patient_list)
#89887

#构建AnnData
count_t = pd.DataFrame(merged_df.values.T,index=merged_df.columns,columns=merged_df.index)
counts=csr_matrix(count_t.to_numpy(),dtype=np.float32)
data = {'patient':patient_list}
meta = pd.DataFrame(index=merged_df.columns,data=data)
var = pd.DataFrame(index=merged_df.index,data={'gene':merged_df.index})
adata = ad.AnnData(X=counts,obs=meta,var=var)
# AnnData object with n_obs × n_vars = 89887 × 29527
#     obs: 'patient'
#     var: 'gene'
adata.X.todense()[0:3,0:3]

#Preprocessing
sc.pl.highest_expr_genes(adata, n_top=20,save="_HighExpr.png")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=5000)
sc.pp.filter_cells(adata, max_counts=30000)
sc.pp.filter_genes(adata, min_cells=3)

##QC
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save="_QC.png"
)
#Remove cells that have too many mitochondrial genes expressed or too many total counts
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt",save="_1.png")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts",save="_2.png")
adata = adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.pct_counts_mt < 30, :].copy()

#Total-count normalize 
sc.pp.normalize_total(adata, target_sum=1e4)
#Logarithmize 
sc.pp.log1p(adata)
#Identify highly-variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata,save=".png")


adata.raw = adata
adata.uns['rwa_var'] = adata.var

adata = adata[:, adata.var.highly_variable]
#Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
#Scale the data to unit variance.
sc.pp.scale(adata, max_value=10)

#Principal component analysis
sc.tl.pca(adata, svd_solver="arpack")
sc.pl.pca(adata, color="patient",save="_patient.png")

sc.pl.pca_variance_ratio(adata, log=True,save="_pca_variance_ratio.png")
results_file="/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071_processed.h5"
adata.write(results_file)


adata=sc.read("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071_processed.h5")
#Computing the neighborhood graph and UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["patient"],save="_patient.png")

#remove the effect of patient
import harmonypy
harmony_out = harmonypy.run_harmony(adata.obsm["X_pca"].astype(np.float64), adata.obs, 'patient')
adata.obsm['X_pca_harmony'] = harmony_out.Z_corr.T

#after harmony, check the patient effect
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["patient"],save="_patient_after_harmony.png")

# Clustering the neighborhood graph
sc.tl.leiden(adata, resolution=1)
sc.pl.umap(adata, color=["leiden"],save="_leiden_resolution1.png",legend_loc="on data")
sc.pl.umap(adata, color=["leiden"],save="_leiden_resolution1_2.png")


#Finding marker genes
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,save='_DEG.png')

results_file="/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071_batch_removed.h5"
adata.write(results_file)


#save the DEG analysis results and mannual annotate the cell type with CellMarker2 database
result = adata.uns["rank_genes_groups"]
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
writer = pd.ExcelWriter('/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/leiden_DEG_res.xlsx')
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


############################################
adata = sc.read("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071_batch_removed.h5")

cell_type_annotation = {'0':'Tumor','1':'Natural Killer','2':'Epithelial','3':'Tumor','4':'Fibroblast',
    '5':'T','6':'Tumor','7':'Fibroblast','8':'B','9':'Plasma',
    '10':'Neutrophils','11':'Fibroblast','12':'Tumor','13':'Ciliated','14':'Alveolar',
    '15':'Artery','16':'Mast','17':'Tumor','18':'Tumor','19':'Tumor',
    '20':'Myeloid','21':'Monocytes',}
adata.obs['cell_type'] = adata.obs['leiden'].map(cell_type_annotation)
adata.obs['cell_type'].to_csv('/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/cell_type_annotation.tsv',sep='\t',index=True,header=False)

##infercnv cell type 
cell_type_annotation = {'0':'malignant','1':'Natural Killer','2':'malignant','3':'malignant','4':'malignant',
    '5':'T','6':'malignant','7':'Fibroblast','8':'B','9':'Plasma',
    '10':'Neutrophils','11':'Neutrophils','12':'malignant','13':'Ciliated','14':'Fibroblast',
    '15':'Artery','16':'Mast','17':'malignant','18':'malignant','19':'malignant',
    '20':'Myeloid','21':'malignant'}
adata.obs['cell_type_infercnv'] = adata.obs['leiden'].map(cell_type_annotation)
mask = adata.obs['cell_type_infercnv'] == 'malignant'
def custom_replace(text):
    text = text.replace("_", "")
    return(text)
adata.obs.loc[mask, 'cell_type_infercnv'] = adata.obs.loc[mask, 'cell_type_infercnv'].astype(str) + '_' + adata.obs.loc[mask, 'patient'].astype(str).apply(custom_replace)
adata.obs
adata.obs['cell_type_infercnv'].to_csv('/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/cell_type_infercnv_annotation.tsv',sep='\t',index=True,header=False)

sc.pl.umap(adata, color=["cell_type"],save="_cell_type.png",legend_loc="on data")

sc.pl.umap(adata, color=["ICK","FOXE1","DAPL1","EXOSC7","KCTD1","COQ3","CALML3","GJB6","C3orf67","SLC4A11","TTC32","RGMA"],save="_final_signature.png",legend_loc="on data")

###save counts matrix and cell type annotation to run inferCNV and CytoTRACE
count_t.iloc[0:3,0:3]
sub_counts = count_t.loc[adata.obs.index.tolist()]
sub_counts.T.to_csv("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071_counts.tsv",sep='\t',header=True,index=True)
