library(data.table)

counts = fread("GSE148071_counts.tsv")
counts = as.data.frame(counts)

selected_cells = fread("cibersortx_input_selected_cells.tsv")
selected_cells = as.data.frame(selected_cells)

rownames(counts) = counts$V1
rownames(selected_cells) = selected_cells$V1

sub_counts = counts[selected_cells$V1]
cell_type=selected_cells$cell_type
sub_counts = rbind(cell_type,sub_counts)
rownames(sub_counts)[1] = "GeneSymbol"

fwrite(sub_counts,"cibersortx_input_counts.tsv",sep='\t',col.names = FALSE, row.names = TRUE)

lusc_deconvolution = fread("cibersortx_results/CIBERSORTx_ lusc_tcga_Results.txt")
luad_deconvolution = fread("cibersortx_results/CIBERSORTx_ luad_tcga_Results.txt")

lusc_tcga = fread('cibersortx_input/TCGA_LUSC_data.txt')
luad_tcga = fread('cibersortx_input/TCGA_LUAD_data.txt')


#combine lusc and luad
all_tcga = cbind(lusc_tcga,luad_tcga[-1])
all_tcga = as.data.frame(all_tcga)
rownames(all_tcga) = all_tcga$Gene
all_deconvolution = rbind(lusc_deconvolution,luad_deconvolution)


#tumor cell markers:
tumor_markers <- read_excel("celltype_DEG_res.xlsx", 
                            sheet = "Sheet13")
tumor_markers$fdr =  p.adjust(tumor_markers$Tumor_p,method="fdr")

filter_gene_tumor = tumor_markers[which((tumor_markers$fdr<1e-5) & (tumor_markers$Tumor_l>=0.5)),]


#correlation with T cell infiltration
#for (i in tumor_markers) {
#  res = cor(as.numeric(lusc_tcga[which(lusc_tcga$Gene==i),-1]),lusc_deconvolution$T,method = 'spearman')
#  print(c(i,res))
#}
#for (i in tumor_markers) {
#  res = cor(as.numeric(luad_tcga[which(luad_tcga$Gene==i),-1]),luad_deconvolution$T,method = 'spearman')
#  print(c(i,res))
#}

tumor_markers = filter_gene_tumor$Tumor_n
filter_tumor_genes_cor_with_T = c()
for (i in tumor_markers) {
  res = cor((all_deconvolution$T+all_deconvolution$B+all_deconvolution$`Follicular B`+all_deconvolution$Mast+all_deconvolution$Myeloid),as.numeric(all_tcga[i,all_deconvolution$Mixture]),method = 'spearman')
  #res = cor(all_deconvolution$T,as.numeric(all_tcga[i,all_deconvolution$Mixture]),method = 'spearman')
  filter_tumor_genes_cor_with_T = c(filter_tumor_genes_cor_with_T,res)
}
filter_tumor_genes_cor_with_T = data.frame(gene=tumor_markers,cor=filter_tumor_genes_cor_with_T)
filter_tumor_genes_cor_with_T = filter_tumor_genes_cor_with_T[which(!is.na(filter_tumor_genes_cor_with_T$cor)),]
filter_tumor_genes_cor_with_T = filter_tumor_genes_cor_with_T[which(filter_tumor_genes_cor_with_T$cor<(-0.2)),]
filter_tumor_genes_cor_with_T = filter_tumor_genes_cor_with_T[order(filter_tumor_genes_cor_with_T$cor),]

signature = filter_tumor_genes_cor_with_T$gene
signature = intersect(signature,stem_signature)
signature
intersect(signature,filter_survival_related_gene$gene)
#RIMKLA




gene_cor_with_T = c()
##gene cor with CD8+ T cell fraction
for (i in all_tcga$Gene) {
  res = cor((all_deconvolution$T+all_deconvolution$Myeloid+all_deconvolution$B+all_deconvolution$`Follicular B`+all_deconvolution$Mast),as.numeric(all_tcga[i,all_deconvolution$Mixture]),method = 'spearman')
  gene_cor_with_T = c(gene_cor_with_T,res)
}

gene_cor_with_T = data.frame(gene= all_tcga$Gene,cor=gene_cor_with_T)
gene_cor_with_T = gene_cor_with_T[which(!is.na(gene_cor_with_T$cor)),]
gene_cor_with_T = gene_cor_with_T[order(gene_cor_with_T$cor,decreasing = T),]

signature = tail(gene_cor_with_T$gene,20)


gene_list = gene_cor_with_T$cor
names(gene_list) = gene_cor_with_T$gene
res <- GSEA(
  gene_list,
  TERM2GENE = gmt
)
dotplot(res)
gseaplot2(res, title = res$Description[1],geneSetID = 1)
res$Description[which(res$NES==max(res$NES))]
#REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL
res$Description[which(res$NES==min(res$NES))]
#REACTOME_RRNA_PROCESSING

signature = gmt$gene[gmt$term=="REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL"]

