Sys.setenv(RETICULATE_PYTHON="/xtdisk/jiapl_group/kanghongen/software/miniconda/envs/R4.2/bin/python") 
library("CytoTRACE")
library(data.table)

lung_cancer_counts = fread("GSE148071_counts.tsv",header=T)
lung_cancer_counts = as.data.frame(lung_cancer_counts)
rownames(lung_cancer_counts) = lung_cancer_counts$V1
lung_cancer_counts = lung_cancer_counts[,-1]

annotation = fread("cell_type_annotation.tsv",header=F)
cell_type_annotation = annotation$V2
names(cell_type_annotation) = annotation$V1
cell_type_annotation[1:4]

results <- CytoTRACE(mat = lung_cancer_counts, ncores=10)
saveRDS(results,"./CytoTRACE_output/CytoTRACE_results.rds")
res = data.frame(cell_barcode=names(results$CytoTRACE),CytoTRACE = as.numeric(results$CytoTRACE),CytoTRACErank=results$CytoTRACErank)
fwrite(res,file="./CytoTRACE_output/CytoTRACE_results.txt",sep="\t")

plotCytoGenes(results, numOfGenes = 10,outputDir = "./CytoTRACE_output/")
# plotCytoTRACE(results, phenotype = cell_type_annotation, gene = "TOP2A",outputDir = "./CytoTRACE_output/")

