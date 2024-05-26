##find cancer stemness signature
from scipy.stats import spearmanr
#只用cancer cell 算与CytoTRACE的相关性，然后筛选基因
sub_adata = adata_raw[adata_raw.obs['cell_type']=="Tumor"]
sub_adata

spearmanr(sub_adata.X.T[0].todense().tolist()[0], sub_adata.obs['CytoTRACE'])

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
