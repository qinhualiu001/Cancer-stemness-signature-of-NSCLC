# export PATH=/software/biosoft/software/python/python2020/bin:$PATH
# source activate infercnv

library('infercnv')

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="GSE148071_counts.tsv",
                                    annotations_file="cell_type_infercnv_annotation.tsv",
                                    delim="\t",
                                    gene_order_file="hg38_gencode_v27.txt",
                                    ref_group_names=c("B","T","Natural Killer","Mast","Fibroblast","Artery","Neutrophils","Ciliated"))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_output_dir",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T
)

#infercnv results analysis: https://zhuanlan.zhihu.com/p/376438151
#

