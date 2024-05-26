import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import scanpy as sc
import os

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300, facecolor="white")

adata = sc.read("GSE148071_batch_removed.h5")
adata

sc.pl.umap(adata, color=["leiden"],save="_leiden_resolution_1.png",legend_loc="on data")
sc.pl.umap(adata, color=["leiden"],save="_leiden_resolution_1_2.png")

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
writer = pd.ExcelWriter('leiden_DEG_res.xlsx')
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

cell_type_annotation = {'0':'Tumor','1':'Myeloid','2':'Tumor','3':'Tumor','4':'Basal',
	'5':'Myeloid','6':'Epithelial','7':'T','8':'Tumor','9':'Alveolar',
	'10':'Fibroblast','11':'B','12':'Follicular B','13':'Tumor','14':'Tumor',
	'15':'Tumor','16':'Ciliated','17':'Fibroblast','18':'Endothelial','19':'Mast',
	'20':'Tumor','21':'Tumor','22':'Tumor','23':'Tumor','24':'Tumor','25':'Epithelial','26':'Tumor','27':'Tumor',
    '28':'Basal','29':'Tumor','30':'Follicular B','31':'Tumor','32':'Tumor','33':'Club'}
adata.obs['cell_type'] = adata.obs['leiden'].map(cell_type_annotation)

ax = sc.pl.umap(adata, color=["cell_type"],title="Cell type",show=False)
ax.legend(loc='right',bbox_to_anchor=(0.97,.25,0.5,0.5),ncol=1,fontsize=9)
plt.savefig('./figures/_cell_type_new.png',bbox_inches='tight')


import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
# Inital setting for plot size
from matplotlib import rcParams
rcParams["figure.figsize"] = (3,3)
#sc.pl.umap(adata, color=["patient"],title="Patient",save="_patient.png")
ax = sc.pl.umap(adata, color="patient", title="Patient", show=False)
# Modify legend
ax.legend(loc='right',bbox_to_anchor=(1.03,.22,0.5,0.5),ncol=2,fontsize=6)
plt.savefig('./figures/umap_patient.png',bbox_inches='tight')


#plot marker gene
sc.pl.umap(adata, color=["EGFR","EPCAM","CLDN5","CAPS","LUM","COL1A1","ACTA2","MYH11","LYZ","HLA-DRA","CD2","CD79A","CLDN18","GATA2","CSF3R","FDCSP"],save="_genes.png",legend_loc="on data")

#plot the cell type fraction
import matplotlib.pyplot as plt
#adata.obs.groupby(['cell_type']).size()
a = adata.obs.groupby(['patient','cell_type']).size()
cell_fraction = a.reset_index(name='Count').pivot_table(index='cell_type', columns='patient', values='Count', fill_value=0)
cell_fraction.sum(axis=0)
data = cell_fraction.T
plt.figure(figsize=(15, 7))
labels = data.index
width=0.8
bottom_y = [0]*len(labels)
for i in data.columns:
    y = data[i]/data.sum(axis=1)
    plt.bar(labels,y,width,bottom=bottom_y,label=i)
    bottom_y = y+bottom_y
plt.legend(loc='right',bbox_to_anchor=(.64,.22,0.5,0.5),ncol=1)
plt.tick_params(axis='x', rotation=90)
plt.xlabel('Patients')
# 定义一个函数，将0-1之间的值转换为百分比形式
def to_percent(y, position):
    # 将数值乘以100并四舍五入到2位小数，然后添加百分号
    return f"{y * 100:.1f}%"

# 使用FuncFormatter类将纵坐标转换为百分比形式
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(to_percent))
plt.ylabel('Cell type fraction')
plt.grid(False)
plt.savefig('./figures/patient_cellfraction.png',bbox_inches='tight')


##CytoTRACE
cytotrace = pd.read_csv('CytoTRACE_output/CytoTRACE_results.txt',sep='\t')
adata.obs = adata.obs.merge(cytotrace.set_index('cell_barcode'), left_index=True,right_on='cell_barcode')
sc.pl.umap(adata, color=["CytoTRACE"],save="_CytoTRACE.png")

#Finding marker genes of cell type
sc.tl.rank_genes_groups(adata, "cell_type", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,save='_DEG_celltype.png')

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
writer = pd.ExcelWriter('celltype_DEG_res.xlsx')
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



marker_genes_dict = {
    "Tumor": ['AKR1C1','EPCAM','MYC','KRT17'],
    "Endothelial": ['FLT1','VWF','PECAM1','KDR','CDH5'],
    "Basal":["SPRR1B","KRT6A","CSTA","KRT16"],
    "Epithelial": ['CLDN4','TACSTD2','CEACAM6','ELF3'],
    "Fibroblast": ['COL1A1','COL1A2','DCN','LUM'],
    "Myeloid": ['TYROBP','FCER1G','CTSS','LAPTM5','AIF1'],
    "T": ['CD2','CD3D','TRBC2','CD3G','CD3D','GZMA'],
    "B": ['MZB1','CD79A','IGHG1','IGKC','JCHAIN'],
    "Alveolar": ['NAPSA','MUC1','CLDN18','AQP4'],
    "Mast" : ['TPSB2','TPSAB1','CPA3','MS4A2'],
    "Club": ['SCGB1A1','SCGB3A1','SCGB3A2'],
    "Follicular B":['HLA-DRB1','HLA-DPA1','HLA-DRA','CD74'],
    "Ciliated":['TSPAN1','CAPS','TPPP3','PIFO'],
}
sc.pl.heatmap(adata, marker_genes_dict, groupby="cell_type",show_gene_labels=True, cmap="viridis", dendrogram=True,save="_marker_heatmap.png")

sc.pl.dotplot(adata, marker_genes_dict, "cell_type", dendrogram=True,save="_dotplot_celltype.png")

sc.pl.rank_genes_groups_dotplot(
    adata,
    n_genes=4,
    values_to_plot="logfoldchanges",
    min_logfoldchange=3,
    vmax=7,
    vmin=-7,
    cmap="bwr",
    save='_foldchange_marker.png'
)

subtype_dict={'patient_16':'LUAD','patient_28':'LUAD','patient_9':'LUAD','patient_8':'LUAD',
'patient_20':'LUAD','patient_24':'LUAD','patient_2':'LUAD','patient_13':'LUAD',
'patient_33':'LUAD','patient_34':'LUAD','patient_38':'LUAD','patient_29':'LUAD',
'patient_18':'LUSC','patient_12':'LUAD','patient_42':'NSCLC','patient_5':'LUAD',
'patient_32':'LUAD','patient_21':'LUAD','patient_35':'LUAD','patient_39':'LUAD',
'patient_1':'LUSC','patient_3':'LUSC','patient_7':'LUSC','patient_10':'LUSC',
'patient_11':'NSCLC','patient_17':'LUSC','patient_19':'LUSC','patient_23':'LUSC',
'patient_25':'LUSC','patient_26':'LUSC','patient_27':'LUSC','patient_31':'LUSC',
'patient_36':'LUSC','patient_4':'LUSC','patient_22':'LUSC','patient_30':'LUSC',
'patient_37':'LUSC','patient_14':'LUSC','patient_15':'LUSC','patient_41':'LUSC',
'patient_40':'LUSC','patient_6':'LUSC'}

adata.obs["subtype"] = adata.obs["patient"].map(subtype_dict)
sc.pl.umap(adata, color=["subtype"],save="_subtype.png")

adata.obs['cell_type'].value_counts()
# cell_type
# Tumor           42014
# Myeloid         14539
# Basal            5312
# Fibroblast       5270
# Epithelial       5073
# T                4357
# Alveolar         4287
# B                3093
# Follicular B     2524
# Ciliated         1299
# Endothelial      1179
# Mast              858
# Club               14
# Name: count, dtype: int64


## 随机选择一些细胞作为cibersortx的输入
# 定义一个函数，根据分组的大小来决定取多少行
def sample_func(group):
    if len(group) > 20000:
        return group.sample(n=min(4000, len(group)), random_state=1)
    elif len(group) > 500:
        return group.sample(n=500, random_state=1)
    else:
        return group.sample(len(group), random_state=1)

# 根据'cell_type'列进行分组，并在每个分组中应用sample_func函数
selected_rows = adata.obs.groupby('cell_type').apply(sample_func)

# 将结果转换为DataFrame，以便更好地查看结果
selected_rows = selected_rows.droplevel(0)

selected_rows.to_csv('cibersortx_input_selected_cells.tsv',sep='\t',header=True,index=True)

#selected_rows = pd.read_csv('cibersortx_input_selected_cells.tsv',sep='\t',index_col=0)
#all_counts = pd.read_csv('GSE148071_counts.tsv',sep='\t',index_col=0)


adata_raw = adata.raw.to_adata()
len(adata_raw.X[0].todense().tolist()[0])
##find cancer stemness signature
from scipy.stats import spearmanr
#只用cancer cell 算与CytoTRACE的相关性，然后筛选基因
sub_adata = adata_raw[adata_raw.obs['cell_type']=="Tumor"]
sub_adata
# 创建一个空列表来保存每一列的Spearman相关系数
spearman_correlations = []
spearman_correlations_p = []
# 遍历数组的每一列
for column in sub_adata.X.T:
    # 计算当前列与向量的Spearman相关系数
    corr, p = spearmanr(column.todense().tolist()[0], sub_adata.obs['CytoTRACE'])
    # 将相关系数添加到列表中
    spearman_correlations.append(corr)
    spearman_correlations_p.append(p)
    
gene_cor_with_cytotrace_res = pd.DataFrame({'gene':sub_adata.var.index.tolist(),'spearman_cor':spearman_correlations,'cor_p':spearman_correlations_p})
gene_cor_with_cytotrace_res.to_csv('gene_cor_with_cytotrace_res.tsv',sep='\t',header=True)

results_file="/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071_cell_type.h5"
adata.write(results_file)

