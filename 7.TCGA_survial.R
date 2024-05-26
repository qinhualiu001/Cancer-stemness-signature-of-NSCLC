#TCGA survial
library(survminer)
library(survival)
library(png)

load("TCGA_data/TCGA_LUAD_data.Rdata")
load("TCGA_data/TCGA_LUSC_data.Rdata")



HR_95CI <- function(x){ 
  x <- summary(x)
  HR <-signif(x$coef[2], digits=2)
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
  res<- c(HR,HR.confint.lower,HR.confint.upper)
  names(res)<-c("Hazard Ratio","95% Upper CI","95% lower CI")
  return(res)
}
dir.create("figures/Results_survival_analysis")
dir.create("figures/Results_survival_analysis/OS")
dir.create("figures/Results_survival_analysis/PFS")
Biomarker_OS <- function(data,name){
  print(name)
  dir.create(paste0("figures/Results_survival_analysis/OS/",name))
  dir.create(paste0("figures/Results_survival_analysis/OS/",name,""))
  OS_results <- c()
  #signature = c("SLC4A11","TTC32","RGMA" )
  signature = c("ICK","FOXE1","DAPL1","EXOSC7","KCTD1","COQ3","CALML3","GJB6","C3orf67","SLC4A11","TTC32","RGMA" )
  signature_score = gsva(data$TPM, list(signature), method="ssgsea", verbose=T)
  data$Landscape = data.frame(Sample=colnames(signature_score),marker=as.numeric(signature_score))
  expMarker <- merge.data.frame(data$Clinical, data$Landscape,by="Sample",all.y = T)
  expMarker <- expMarker[which(!is.na(expMarker$OS)),]
  expMarker$OS <- as.numeric(expMarker$OS)
  expMarker$OS_CNSR <- as.numeric(expMarker$OS_CNSR)
  expMarker$Class <- ifelse(expMarker$marker > mean(expMarker$marker),"Signature High","Signature Low")
  fit = survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("figures/Results_survival_analysis/OS/",name,"/marker.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,pval.method = TRUE,surv.median.line="hv",conf.int=T,
                    ggtheme = theme_classic2(base_family = "Times New Roman"),
                    font.main = c(16, "bold"),
                    font.title= c(16, "bold"),
                    font.subtitle= c(16, "bold"),
                    font.caption= c(16, "bold"),
                    font.x = c(16, "bold"),
                    font.y = c(16, "bold"),
                    font.tickslab = c(14, "bold"),
                    font.legend = c(16, "bold"))

  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
    ggsave(paste0("figures/Results_survival_analysis/OS/",name,"/forest_marker.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  return(result)
}

res = Biomarker_OS(TCGA_LUSC_data,"TCGA_LUSC_data")

res

all_clinical = rbind(TCGA_LUSC_data$Clinical,TCGA_LUAD_data$Clinical)

survival_related_gene = data.frame()
for(i in all_tcga$Gene){
  a = all_tcga[which(all_tcga$Gene==i),]
  b = data.frame(Sample=colnames(a)[-1],gene=as.numeric(a[-1]))
  b = b[which(!is.na(b$gene)),]
  c = merge(b,all_clinical,by="Sample",all.y = T)
  c = c[which(!is.na(c$gene)),]
  c = c[which(!is.na(c$OS)),]
  c$Class <- ifelse(c$gene > median(c$gene),"marker High","marker Low")
  if(length(unique(c$Class))>1){
    c$OS <- as.numeric(c$OS)
    c$OS_CNSR <- as.numeric(c$OS_CNSR)  
    fit = survfit(Surv(OS,OS_CNSR) ~ Class,data=c)
    p = surv_pvalue(fit,data=c)
    res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=c)
    hr = HR_95CI(res.cox)
    result = c(hr,p$pval,i)
    names(result) = c(names(hr),"p_value","gene")
    survival_related_gene = rbind(survival_related_gene,result)
    colnames(survival_related_gene) = names(result)
  }
}

filter_survival_related_gene = survival_related_gene[which((survival_related_gene$p_value<0.05) & (survival_related_gene$`Hazard Ratio`>1)),]
#filter_survival_related_gene$gene = gsub(".","_",filter_survival_related_gene$gene,fixed = T)


signature = intersect(stem_signature,filter_tumor_genes_cor_with_T$gene)

signature = intersect(signature,filter_survival_related_gene$gene)
signature
#[1] "ICK"     "FOXE1"   "DAPL1"   "EXOSC7"  "KCTD1"   "TTC32"   "COQ3"   
#[8] "CALML3"  "GJB6"    "SLC4A11" "C3orf67" "RGMA" 


# venn plot of the signature genes
library("VennDiagram")

venn.diagram(x=list(stem_signature,filter_tumor_genes_cor_with_T$gene,filter_survival_related_gene$gene),
             scaled = F, # 根据比例显示大小
             alpha= 0.5, #透明度
             lwd=1,lty=1,col=c('#153448','#1A4D2E',"#874CCC"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex = 2, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=c('#153448','#1A4D2E',"#874CCC"), # 填充色
             category.names = c("Cancer stemness signature", "Immune related signature","Survival related siagnture") , #标签名
             cat.dist = 0.02, # 标签距离圆圈的远近
             cat.pos = c(-10, 10, -180), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex = 2, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             output=TRUE,
             filename='./figures/venn.png',# 文件保存
             imagetype="png",  # 类型（tiff png svg）
             resolution = 300,  # 分辨率
             compression = "lzw"# 压缩算法
)





