## merged spearman and lasso associations between cell type and tumor microbiome are needed for 32 cancer type
## candidate_assciations.sh can be used to get the above results, col 2 is cancer name

library(dplyr)

## get asso list
cell_taxa_asso <- read.table("cell_taxa_asso")
cell_taxa_asso_new <- cell_taxa_asso %>%mutate_all(funs(gsub("k__.+g__", "", .)))
asso_list <- data.frame(paste(cell_taxa_asso_new[,1], cell_taxa_asso_new[,3], sep="--"))
colnames(asso_list)[1] <- "asso"

## count and define pancancer asso
asso_count <- count(asso_list_new,asso)
pancancer_asso <- filter(asso_count, n >= 8)
pancancer_asso_new <- pancancer_asso[!grepl("Other",pancancer_asso$asso), ]

associations <- cbind(asso_list, cell_taxa_asso_new[,2], cell_taxa_asso_new[,4])
colnames(associations)[1] <- "V1"
colnames(associations)[2] <- "V2"
colnames(associations)[3] <- "V3"

cancername <- c( "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC",
           "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG",
           "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD",
           "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
           "UCEC","UCS","UVM")


## read pancancer asso
pan_cancer_cell_microbes_asso_list <- data.frame(pancancer_asso_new$asso)

pan_cancer_cell_microbes_asso_list <- as.factor(pan_cancer_cell_microbes_asso_list[,1])
  
cell_associations <- data.frame(matrix(nrow=length(pan_cancer_cell_microbes_asso_list),ncol=length(cancername)))
colnames(cell_associations) = cancername
rownames(cell_associations) = pan_cancer_cell_microbes_asso_list

#associations <- read.table("cell-microbe-asso", sep = "\t")

## get spearman rho
pick_rho <- function (asso, cancer){
    pick_asso <- subset(associations, associations$V1==asso & associations$V2==cancer)
    return (pick_asso$V3)
}


for (asso in pan_cancer_cell_microbes_asso_list){
  for (cancer in cancername){
    if (length(pick_rho (asso, cancer)) == 0 ){
      cell_associations[asso, cancer] <- "0"
    }
    else{cell_associations[asso, cancer] <- pick_rho (asso, cancer)}
  }
}
write.table(cell_associations,"cell_associations",sep="\t")

