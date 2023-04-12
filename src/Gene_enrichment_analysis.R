library(clusterProfiler)
library(org.Hs.eg.db)
library(DO.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)

gg <- read.table("gene_list.txt")
gg_table=bitr(gg$V1,fromType='SYMBOL',toType='ENTREZID',OrgDb='org.Hs.eg.db')

##GO analysis
ego <- enrichGO(
  gene          = gg_table$ENTREZID,
  keyType = "ENTREZID",
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.25,
  qvalueCutoff  = 1,
  readable      = TRUE)
write.table(ego,"GO")

barplot(ego, showCategory = 10)
dotplot(ego, showCategory = 10)
emapplot(pairwise_termsim(ego), showCategory = 30)
cnetplot(pairwise_termsim(ego), showCategory = 10)


## KEGG analysis
kk <- enrichKEGG(gene = gg_table$ENTREZID, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)


barplot(kk, showCategory = 10)
dotplot(kk, showCategory = 10)
emapplot(pairwise_termsim(kk), showCategory = 30)
cnetplot(pairwise_termsim(kk), showCategory = 5)
