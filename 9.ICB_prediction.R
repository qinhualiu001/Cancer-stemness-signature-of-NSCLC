##ICB response prediction
load("./ICB_datasets/Jung_2019.Rdata")
load('./ICB_datasets/Braun_2020.Rdata')
load('./ICB_datasets/Gide_2019.Rdata')
load('./ICB_datasets/Kim_2018.Rdata')
load('./ICB_datasets/Mariathasan_2018.Rdata')
load('./ICB_datasets/Riaz_2017.Rdata')

library(ggplot2)
library(ggpubr)
library(ggsci)
library(pROC)
library("GSVA")

my_comparisons = list(c("NonResponder","Responder"))

##ssGSEA based signature_score
ssgsea_score = function(data,name){
  signature = c('SLC4A11','TTC32', 'RGMA')
  signature_score = gsva(data$TPM, list(signature), method="ssgsea", verbose=T)
  expMarker = merge(data$Samples,data.frame(Sample=colnames(signature_score),Signature_score=as.numeric(signature_score)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  expMarker$Class = expMarker$Resp_NoResp
  expMarker$Class[which(expMarker$Class=="No_Response")] = "NonResponder" 
  expMarker$Class[which(expMarker$Class=="Response")] = "Responder" 
  g = ggboxplot(expMarker,x="Class",y="Signature_score",color = "Class",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("figures/",name,"/ssgsea_score_wilcox_test.png.png"),width = 6,height = 6, bg="white")
  #res = wilcox.test(signature_score ~ Resp_NoResp, data = expMarker,alternative = "less")
  res = wilcox.test(Signature_score ~ Class, data = expMarker,alternative = "greater")
  auc = roc(Class ~ Signature_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("figures/",name,"/ssgsea_score_AUC.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}

ssgsea_score(Jung_2019,"Jung_2019")
#p_value        AUC 
#0.01083882 0.78289474

ssgsea_score(Gide_data,"Gide_data")
#p_value       AUC 
#0.1640731 0.5619748 

ssgsea_score(Riaz_data,"Riaz_data")
#p_value       AUC 
#0.4020025 0.5148237 

ssgsea_score(Braun_data,"Braun_data")
#p_value       AUC 
#0.2958173 0.5281487 

ssgsea_score(Kim_data,"Kim_data")
#p_value       AUC 
#0.3017013 0.5530303

ssgsea_score(Mariathasan_data,"Mariathasan_data")
#p_value       AUC 
#0.4588534 0.5041560 




##################
##innate anti-PD-1 resistance (IPRES) signature
# We row-normalized
# the GSVA scores of each gene set in the IPRES signature across the samples
# from the four cohorts
# The IPRES (enrichment) score was defined as the average
# Z score across all gene sets in the IPRES signature

IPRES_score = function(data,name){
  library("qusage")
  IPRES_signatures = read.gmt("marker/IPRES_signatures.gmt")
  gsva.es <- gsva(data$TPM, IPRES_signatures, method="ssgsea", verbose=T)
  gsva.es = scale(t(gsva.es))
  IPRES_score = apply(gsva.es,1,mean)
  expMarker = merge(data$Samples,data.frame(Sample=names(IPRES_score),IPRES_score=IPRES_score),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  expMarker$Class = expMarker$Resp_NoResp
  expMarker$Class[which(expMarker$Class=="No_Response")] = "NonResponder" 
  expMarker$Class[which(expMarker$Class=="Response")] = "Responder" 
  library(ggpubr)
  g = ggboxplot(expMarker,x="Class",y="IPRES_score",color = "Class",
                palette = "npg",add = "jitter") +  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("figures/",name,"/IPRES_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(IPRES_score ~ Class, data = expMarker,alternative = "greater")
  auc = roc(Class ~ IPRES_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("figures/",name,"/IPRES_score_AUC.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}
IPRES_score(Jung_2019,"Jung_2019")
#p_value       AUC 
#0.2407302 0.5921053 

C_ECM_score = function(data,name){
  library(readr)
  C_ECM_genes <- read_delim("marker/C_ECM.txt", delim = ";", 
                            escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)
  geneSet = list()
  geneSet[["C-ECM"]] = C_ECM_genes$X1
  C_ECM_score <- gsva(data$TPM, geneSet, method="ssgsea", verbose=T)
  expMarker = merge(data$Samples,data.frame(Sample=colnames(C_ECM_score),C_ECM_score=as.numeric(C_ECM_score[1,])),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  expMarker$Class = expMarker$Resp_NoResp
  expMarker$Class[which(expMarker$Class=="No_Response")] = "NonResponder" 
  expMarker$Class[which(expMarker$Class=="Response")] = "Responder" 
  library(ggpubr)
  g = ggboxplot(expMarker,x="Class",y="C_ECM_score",color = "Class",
                palette = "npg",add = "jitter") +  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("figures/",name,"/C-ECM_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(C_ECM_score ~ Class, data = expMarker,alternative = "greater")
  auc = roc(Class ~ C_ECM_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("figures/",name,"/C-ECM_score_AUC.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}
C_ECM_score(Jung_2019,"Jung_2019")
#p_value       AUC 
#0.1188352 0.6513158 

##########pan-fibroblast TGFβ response signature (F-TBRS), based on PCA
F_TBRS_score = function(data,name){
  F_TBRS_genes = c("ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CCN2", "CTPS1",
                   "RFLNB", "FSTL3", "HSPB1", "IGFBP3", "PXDC1", "SEMA7A", "SH3PXD2A", "TAGLN", 
                   "TGFBI", "TNS1", "TPM1")
  m = data$TPM
  m <- t(scale( t( m ),
                center=TRUE, 
                scale=TRUE)
  )
  m2 = m[intersect(F_TBRS_genes,rownames(data$TPM)),]
  ##' Calculate score across genes and samples
  gsScore <- function(gm, summarizationFunction="PC") {
    if (summarizationFunction == "PC") {
      pc <- prcomp(t(gm),
                   retx=TRUE)
      gss <- pc$x[,1] * sign(cor(pc$x[,1], colMeans(gm)))
    } else {
      gss <- colMeans(gm)
    }
    return(gss)
  }
  F_TBRS_score = gsScore(m2) 
  expression = data.frame(Sample=names(F_TBRS_score),F_TBRS=as.numeric(F_TBRS_score))
  expMarker = merge(data$Samples,expression,by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  expMarker$Class = expMarker$Resp_NoResp
  expMarker$Class[which(expMarker$Class=="No_Response")] = "NonResponder" 
  expMarker$Class[which(expMarker$Class=="Response")] = "Responder" 
  g = ggboxplot(expMarker,x="Class",y="F_TBRS",color = "Class",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("figures/",name,"/F-TBRS_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(F_TBRS ~ Class, data = expMarker,alternative = "greater")
  auc = roc(Class ~ F_TBRS, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("figures/",name,"/F-TBRS_score_AUC.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}
F_TBRS_score(Jung_2019,"Jung_2019")
#p_value        AUC 
#0.07343851 0.68421053 

TIDE = function(data,name){
  library(readr)
  file_path = paste0("marker/",name,"_TIDE_output.csv")
  res <- read_csv(file_path)
  expMarker = merge(data$Samples,data.frame(Sample=res$Patient,TIDE=res$TIDE),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  expMarker$Class = expMarker$Resp_NoResp
  expMarker$Class[which(expMarker$Class=="No_Response")] = "NonResponder" 
  expMarker$Class[which(expMarker$Class=="Response")] = "Responder" 
  g = ggboxplot(expMarker,x="Class",y='TIDE',color = "Class",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("figures/",name,"/TIDE_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(TIDE ~ Class, data = expMarker,alternative = "greater")
  auc = roc(Class ~ TIDE, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("figures/",name,"/TIDE_score_AUC.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}
TIDE(Jung_2019,"Jung_2019")
#p_value        AUC 
#0.01376349 0.77631579 

############ EMT/Stroma_core signature: 8 genes were included in the EdgeSeq expression panel 
##Patients with high CD8 infiltration and low EMT/Stromal core gene expression had the highest response rates and longest PFS and OS
EMT_Stroma_core_marker = function(data,name){
  EMT_Stroma_core_signature = c("FLNA","EMP3","CALD1","FN1","FOXC2","LOX","FBN1","TNC")
  expression = data$TPM[intersect(EMT_Stroma_core_signature,rownames(data$TPM)),]
  expression = as.data.frame(t(expression))
  expression["EMT_Stroma_core_signature"] = apply(expression, 1, function(x){return(mean(as.numeric(x)))})
  expression["Sample"] = rownames(expression)
  expMarker = merge(data$Samples,expression[c("Sample","EMT_Stroma_core_signature")],by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  expMarker$Class = expMarker$Resp_NoResp
  expMarker$Class[which(expMarker$Class=="No_Response")] = "NonResponder" 
  expMarker$Class[which(expMarker$Class=="Response")] = "Responder" 
  g = ggboxplot(expMarker,x="Class",y="EMT_Stroma_core_signature",color = "Class",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("figures/",name,"/ESCS_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(EMT_Stroma_core_signature ~ Class, data = expMarker,alternative = "greater")
  auc = roc(Class ~ EMT_Stroma_core_signature, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("figures/",name,"/ESCS_score_AUC.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}
EMT_Stroma_core_marker(Jung_2019,"Jung_2019")
#p_value       AUC 
#0.2092119 0.6052632 

df = data.frame(Biomarker=c("Cancer stemness score","IPRES score","C-ECM score","F-TBRS score","TIDE score","ESCS score"),
                neg_log10_p_value = -log10(c(0.01083882,0.2407302,0.1188352,0.07343851,0.01376349,0.2092119)),
                AUC = c(0.78289474,0.5921053,0.6513158,0.68421053,0.77631579,0.6052632))

df_long <- tidyr::gather(df, key = "Type", value = "value", -Biomarker)
df_long$Type[which(df_long$Type=="neg_log10_p_value")] = "-log10(p)"

ggplot(df_long, aes(x = Biomarker, y = value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Biomarker", y = "Value", fill = "Type") +
  scale_fill_manual(values=c("#F27BBD","#8B93FF")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
  
ggsave(filename = "figures/Jung_2019/benchmark.png",width = 4,height = 4)




