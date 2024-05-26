library("Seurat")
library("harmony")
library("data.table")
library("ggplot2")

# read the single cell dataset "GSE131907"
GSE131907_log2tpm = readRDS("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds")
GSE131907_annotation = fread("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE131907_Lung_Cancer_cell_annotation.txt")
GSE131907_log2tpm = as.matrix(GSE131907_log2tpm)
GSE131907_annotation = as.data.frame(GSE131907_annotation)
rownames(GSE131907_annotation)=GSE131907_annotation$Index
GSE131907 = CreateSeuratObject(counts=GSE131907_log2tpm,project = "GSE131907")
GSE131907@assays$RNA$counts[1:3,1:3]
GSE131907@meta.data[1:3,1:3]
GSE131907_annotation_sorted = GSE131907_annotation[match(rownames(GSE131907@meta.data),GSE131907_annotation$Index),]
combined_meta = cbind(GSE131907@meta.data,GSE131907_annotation_sorted)
GSE131907@meta.data = combined_meta
save(GSE131907,file="/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE131907.Rdata")


# read the single cell dataset "GSE131907"
all_files = list.files(path="/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071",pattern="txt")
GSE148071_all = data.frame()
for(i in 1:42){
  patient = fread(paste0("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071/",all_files[i]))
  genes = patient$V1
  patient = patient[,-1]
  patient = as.data.frame(patient)
  rownames(patient) = genes
  patient = t(patient)
  GSE148071_all = rbind.data.frame(GSE148071_all,patient)
}
GSE148071 = CreateSeuratObject(counts=t(as.matrix(GSE148071_all)), project = "GSE148071")
GSE148071@meta.data[1:3,1:3]
# add annotation of patient source
GSE148071_patient_cells_num = c()
for(i in 1:42){
  patient = fread(paste0("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071/",all_files[i]))
  num = ncol(patient)-1
  GSE148071_patient_cells_num = c(GSE148071_patient_cells_num,rep(paste0("patient_",i),num))
}
GSE148071@meta.data$patient = GSE148071_patient_cells_num
save(GSE148071,file="/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071.Rdata")



load("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE131907.Rdata")
#Quality control
GSE131907[["percent.mt"]] <- PercentageFeatureSet(GSE131907, pattern = "^MT-")
head(GSE131907@meta.data, 5)
##filter poor quality cells
GSE131907 <- subset(GSE131907, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20 & nCount_RNA < 150000 & nCount_RNA > 1000)
## only maintain tLung, mLN, tL/B, nLung, nLN
GSE131907 <- subset(GSE131907, subset = `Sample_Origin` %in% c("nLN","mLN","tL/B","tLung","nLung"))
##GSE131907 has been normalized and log2 transformated
GSE131907@assays$RNA$data = GSE131907@assays$RNA$counts

load("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/GSE148071.Rdata")
#Quality control
GSE148071[["percent.mt"]] <- PercentageFeatureSet(GSE148071, pattern = "^MT-")
head(GSE148071@meta.data, 5)
##filter poor quality cells
GSE148071 <- subset(GSE148071, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20 & nCount_RNA < 150000 & nCount_RNA > 1000)
GSE148071 <- NormalizeData(GSE148071, normalization.method = "LogNormalize", scale.factor = 10000)


overlape_genes = intersect(rownames(GSE148071),rownames(GSE131907))
# union_genes = union(rownames(GSE131907),rownames(GSE148071))
counts = cbind(GSE148071@assays$RNA$data[overlape_genes,],GSE131907@assays$RNA$data[overlape_genes,])

Lung_merged = CreateSeuratObject(counts=counts,project = "Lung_merged")
meta_data = data.frame(dataset=c(rep("GSE148071",ncol(GSE148071)),rep("GSE131907",ncol(GSE131907))),
  patient=c(GSE148071@meta.data$patient,GSE131907@meta.data$Sample_Origin),
  percent.mt=c(GSE148071@meta.data$percent.mt,GSE131907@meta.data$percent.mt),
  celltype_original=c(rep("NA",ncol(GSE148071)),GSE131907@meta.data$Cell_type),
  cell_subtype_original=c(rep("NA",ncol(GSE148071)),GSE131907@meta.data$Cell_subtype))
merged_meta_data = cbind(Lung_merged@meta.data,meta_data)
Lung_merged@meta.data = merged_meta_data
Lung_merged@assays$RNA$data = Lung_merged@assays$RNA$counts

Lung_merged <- FindVariableFeatures(Lung_merged, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Lung_merged)
Lung_merged <- ScaleData(Lung_merged, features = all.genes)
Lung_merged <- RunPCA(Lung_merged, npcs=20, features = VariableFeatures(object = Lung_merged))

plot <- DimPlot(object=Lung_merged,reduction="pca",pt.size=.1,group.by="dataset")
ggsave("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/PCA_before_integrated.jpg", plot)

pc.num=1:20
test_Lung_merged <- RunTSNE(Lung_merged, dims = pc.num,check_duplicates = FALSE)
plot = DimPlot(test_Lung_merged, reduction = "tsne", group.by='dataset')
ggsave("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/tSNE_before_integrated.jpg", plot)

test_Lung_merged<- RunUMAP(Lung_merged, dims = pc.num)
plot = DimPlot(test_Lung_merged, reduction = "umap", group.by='dataset')
ggsave("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/UMAP_before_integrated.jpg", plot)

Lung_merged <- RunHarmony(Lung_merged, group.by.vars = "dataset")
plot = DimPlot(Lung_merged, reduction = "harmony", pt.size=.1, group.by='dataset')
ggsave("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/after_harmony_integrated.jpg", plot)

Lung_merged <- RunUMAP(Lung_merged, reduction = "harmony",dims=1:20)
plot = DimPlot(Lung_merged, reduction = "umap", pt.size=.1, group.by='dataset')
ggsave("/p300s/jiapl_group/kanghongen/projects/TISCH/Lung_cancer/UMAP_after_harmony_integrated.jpg", plot)



