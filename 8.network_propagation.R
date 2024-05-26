library(igraph)
library(org.Hs.eg.db) 
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

network = readRDS('Combined_STRINGv11_OTAR281119_FILTER.rds')
net=graph_from_data_frame(d=network,directed=F)
E(net)$weight=as.numeric(as.character(network[,"combined_score"]))
net.clean=igraph::simplify(net,
                           remove.loops = T,
                           remove.multiple = T ,
                           edge.attr.comb = c(weight="max","ignore"))
net.clean

#utils
# Gene symbol search function
keytypes(org.Hs.eg.db)

# Gene symbol search function
ENSEMBLIDtoSYMBOL <- function(x){
  vec=mapIds(org.Hs.eg.db,keys=x,column="SYMBOL",keytype="ENSEMBL",multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){vec[y]=x[y]}
  return(vec)
}

ENSEMBLtoENTREZ <- function(x){
  vec=mapIds(org.Hs.eg.db,keys=x,column="ENTREZID",keytype="ENSEMBL",multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){vec[y]=x[y]}
  return(vec)
}

#### seed to run page_rank
all_nodes = names(V(net.clean))
all_nodes=data.frame(ensembl=all_nodes,symbol= ENSEMBLIDtoSYMBOL(all_nodes),seed=0)
all_nodes$seed[which(is.element(all_nodes$symbol,signature))] = 1

#all_nodes$seed[which(is.element(all_nodes$symbol,c("PDCD1","CD274")))] = 1
page.rank=page_rank(net.clean, personalized=as.numeric(all_nodes[,"seed"]), weights=E(net.clean)$weight)
signature_network = ENSEMBLIDtoSYMBOL(names(sort(page.rank$vector,decreasing = T)[1:100]))
signature = signature_network[1:3]
signature
#SLC4A11,TTC32, RGMA,

#plot
signature_network[1:10]

sub_net = subgraph(net.clean,vids = names(signature_network[1:20]))

plot(sub_net,vertex.label=signature_network[1:20])




intersect(tail(DEG_res$gene,500),signature)


all_nodes = cbind(all_nodes,page.rank$vector)
colnames(all_nodes)[ncol(all_nodes)]="page.rank"

deg=igraph::degree(net.clean)

##Network filter
node.filter = all_nodes[as.numeric(all_nodes[,"page.rank"])>quantile(as.numeric(all_nodes[,"page.rank"]))[4],]

colnames(node.filter)[1] = "ENSP"
edge.filter=network[  as.character(network[,1])%in%node.filter[,"ENSP"] & 
                           as.character(network[,2])%in%node.filter[,"ENSP"] ,]

node.filter=node.filter[node.filter[,"ENSP"]%in%c(as.character(edge.filter[,1]),as.character(edge.filter[,2])),]
edge.string=edge.filter[,1:2]
net=graph_from_data_frame(d=edge.string,vertices=node.filter,directed=F)
E(net)$weight=as.numeric(as.character(edge.filter[,"combined_score"]))

net.clean=igraph::simplify(net,
                           remove.loops = T,
                           remove.multiple = T ,
                           edge.attr.comb = c(weight="max","ignore"))

cwt=cluster_walktrap(	net.clean, 
                      weights = E(net.clean)$weight, 
                      steps = 6,
                      merges = TRUE, 
                      modularity = TRUE, 
                      membership = TRUE)

degree=igraph::degree(net.clean)

node.filter=cbind(node.filter,degree,cwt$membership,cwt$modularity)

colnames(node.filter)[(ncol(node.filter)-2):ncol(node.filter)]=c("node.degree","cluster.walktrap","modularity.walktrap")
clust=as.matrix(as.data.frame(table(cwt$membership)))
View(node.filter)





# GSEA of page_rank results
gene_list = sort(page.rank$vector,decreasing = T)
names(gene_list) = ENSEMBLIDtoSYMBOL(names(sort(page.rank$vector,decreasing = T)))

#gmt <- read.gmt("pathway/h.all.v2023.2.Hs.symbols.gmt")
gmt <- read.gmt("pathway/c2.cp.reactome.v2023.2.Hs.symbols.gmt")

res <- GSEA(
  gene_list,
  TERM2GENE = gmt
)
dotplot(res)
#p = gseaplot2(res, title = res$Description[which(res$NES==max(res$NES))],geneSetID = 1)
gseaplot2(res, title = res$Description[1],geneSetID = 1)
#ggsave(p,filename = "figures/gsea.png",width = 10,height = 6)
res$Description[1]
#REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE
res$Description[which(res$NES==max(res$NES))]
#REACTOME_COMPLEX_I_BIOGENESIS
res$Description[which(res$NES==min(res$NES))]
#REACTOME_ANTIGEN_PROCESSING_UBIQUITINATION_PROTEASOME_DEGRADATION

signature = gmt$gene[gmt$term=="REACTOME_ANTIGEN_PROCESSING_UBIQUITINATION_PROTEASOME_DEGRADATION"]
intersect(signature,filter_gene_icb$gene)




