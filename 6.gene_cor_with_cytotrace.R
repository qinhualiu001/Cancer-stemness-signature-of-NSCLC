library(data.table)
library(gridExtra)
library(grid)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(pROC)
library("GSVA")
library(readxl)

gene_cor_with_cytotrace_res = fread("gene_cor_with_cytotrace_res.tsv")
gene_cor_with_cytotrace_res = gene_cor_with_cytotrace_res[which(!is.na(gene_cor_with_cytotrace_res$spearman_cor))]
sum(gene_cor_with_cytotrace_res$spearman_cor>0.2)
gene_cor_with_cytotrace_res = gene_cor_with_cytotrace_res[order(gene_cor_with_cytotrace_res$spearman_cor,decreasing = T),]
gene_cor_with_cytotrace_res$fdr = p.adjust(gene_cor_with_cytotrace_res$cor_p,method="fdr")

filter_gene_cytotrace = gene_cor_with_cytotrace_res[which((gene_cor_with_cytotrace_res$fdr<1e-5) & (gene_cor_with_cytotrace_res$spearman_cor>0.1)),]
filter_gene_cytotrace = filter_gene_cytotrace[order(filter_gene_cytotrace$spearman_cor,decreasing = T),]

stem_signature = intersect(filter_gene_cytotrace$gene,filter_gene_tumor$Tumor_n)
stem_signature

intersect(stem_signature,DEG_res$gene[order(DEG_res$log10p)[1:200]])
# "PTPRZ1" "MATN2"  "SEMA6A" "CHST1"  "FZD10"  "SFRP1" 


##GSEA of gene_cor_with_cytotrace
gene_cor_with_cytotrace_res= gene_cor_with_cytotrace_res[order(gene_cor_with_cytotrace_res$spearman_cor,decreasing = T),]
gene_list = gene_cor_with_cytotrace_res$spearman_cor
names(gene_list) = gene_cor_with_cytotrace_res$gene
gmt <- read.gmt("pathway/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
res <- GSEA(
  gene_list,
  TERM2GENE = gmt
)
dotplot(res)
gseaplot2(res, title = res$Description[1],geneSetID = 1)
res$Description[which(res$NES==max(res$NES))]
#GSEA: REACTOME_MITOCHONDRIAL_TRANSLATION

res$Description[which(res$NES==min(res$NES))]
#REACTOME_DISEASES_ASSOCIATED_WITH_SURFACTANT_METABOLISM

signature = gmt$gene[gmt$term=="REACTOME_DISEASES_ASSOCIATED_WITH_SURFACTANT_METABOLISM"]


