
library(lavaan)
library(bruceR)
library(data.table)
library(tidyverse)
library(stargazer)
rm(list = ls()) ## check to make sure no precomputed dataframe before clearing
## make a empty dataframe to accept data
mediation_all_tab<- as.data.frame(matrix(ncol=0,nrow=0))
setwd("C:/Lab/project/pan-cancer-host-microbiome-associations/src")
project <- "TCGA-BLCA"
result.loc <- paste0("../result/SparseCCA_xena/",project,"/")


input_dirname <- paste0("../data/Preprocessed/xena_tcga/",project,"/case/")
filenames <- list.files(input_dirname)
cat("Read data from ",input_dirname," .....")


feat.rppa <- as.matrix(fread(paste0(input_dirname,grep("RPPA",filenames,value = TRUE))),rownames=1)
cat('Feature microbiome otu dimensions:', dim(feat.rppa), '\n')

# rppa_interest_list <-c("BCL2","IRS1","ERALPHA_pS118","ERALPHA") 
rppa_interest_list <-c("HER3","HER2")

feat.rppa.interest <- as.data.frame(feat.rppa[,rppa_interest_list])
# feat.rppa.interest <- sapply(feat.rppa.interest, as.numeric)
rm(feat.rppa)

feat.otu.relt <- read.table(paste0(input_dirname,grep("Kraken",filenames,value = TRUE)))
colnames(feat.otu.relt) <- str_split_fixed(colnames(feat.otu.relt),'g__',2)[,2]
taxa_interest <- c("Lachnoclostridium","Sutterella")
feat.otu.relt.interest <- as.data.frame(feat.otu.relt[,taxa_interest])



cat('Feature microbiome otu dimensions:', dim(feat.otu.relt.interest), '\n')
rm(feat.otu.relt)

feat.gene.expr <- read.table(paste0(input_dirname,grep("Geneexpr",filenames,value = TRUE)))
gene_interest <- c("ESR1","ERBB3")
feat.gene.expr.interest <- as.data.frame(feat.gene.expr[,gene_interest])
rm(feat.gene.expr)

feat.meth <- as.matrix(fread(paste0(input_dirname,grep("CpG",filenames,value = TRUE))),rownames=1)
meth_interest <- c("cg04794420","cg07151117")
feat.meth.interest <- as.data.frame(feat.meth[,meth_interest])

idx <- Reduce(intersect, list(rownames(feat.otu.relt.interest),
                              rownames(feat.gene.expr.interest),
                              rownames(feat.rppa.interest),
                              rownames(feat.meth.interest)))

feat.otu.relt.interest <- feat.otu.relt.interest[idx,]
feat.gene.expr.interest <- feat.gene.expr.interest[idx,]      
feat.rppa.interest <- feat.rppa.interest[idx,]
feat.meth.interest <- feat.meth.interest[idx,]

colnames(feat.rppa.interest) <- paste0("HostP_",colnames(feat.rppa.interest))
colnames(feat.gene.expr.interest) <- paste0("HostT_",colnames(feat.gene.expr.interest))
colnames(feat.meth.interest) <- paste0("HostM_",colnames(feat.meth.interest))

BLCA_BG_deconvolution <- read.table(file="C:/Lab/project/pan-cancer-host-microbiome-associations/data/cell_type/BLCA_BG_deconvolution.t.tsv",sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
## select cell type column

cell_interest <- BLCA_BG_deconvolution[idx,grep("_std$", colnames(BLCA_BG_deconvolution), invert = TRUE, value = TRUE)]

# 
# cell_type_interest <- c("Lymphocytes","Tregs","CD8_T_cells_PD1_high","CD8_T_cells_PD1_low","T_cells","Tregs","Macrophages_M1","Macrophages_M2","Macrophages","CD8_T_cells","NK_cells","Plasma_B_cells","CD4_T_cells")

# cell_interest <- LUSC_BG_deconvolution[idx,cell_type_interest]
# cell_interest <- BLCA_BG_deconvolution[idx,cell_type_interest]

Meddata <- data.frame(feat.otu.relt.interest,feat.meth.interest,feat.gene.expr.interest,feat.rppa.interest,cell_interest)

Meddata_scale <- apply(Meddata,2,scale)

for (i in 1:ncol(cell_interest)){
  cat("Analysis with ",colnames(cell_interest)[i],".\n")
  X_name <- "Lachnoclostridium"
  Y_name <- colnames(cell_interest)[i]
  M1_name <- "HostM_cg04794420"
  M2_name <- "HostT_ERBB3"
  M3_name <- "HostP_HER3"
  Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                                    meds=c(M1_name,M2_name,M3_name),
                                    med.type="serial",
                                    ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
  summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
  # print(summary_tab[8,4] < 0.05)
  summary_tab$Path <- rownames(summary_tab)
  # replace
  summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
  summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
  summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
  summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
  summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
  summary_tab <- summary_tab %>%
    mutate(Significance = case_when(
      (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
      (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
      (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
      (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
      TRUE ~ "ns"
    ))
  summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
  ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
  write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
              quote = F, row.names = T, sep = "\t")
}



for (i in 1:ncol(cell_interest)){
  cat("Analysis with ",colnames(cell_interest)[i],".\n")
  X_name <- "Lachnoclostridium"
  Y_name <- colnames(cell_interest)[i]
  M1_name <- "HostM_cg07151117"
  M2_name <- "HostT_ERBB3"
  M3_name <- "HostP_HER3"
  Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                                             meds=c(M1_name,M2_name,M3_name),
                                             med.type="serial",
                                             ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
  summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
  # print(summary_tab[8,4] < 0.05)
  summary_tab$Path_simple <- rownames(summary_tab)
  summary_tab$Path <- rownames(summary_tab)
  # replace
  summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
  summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
  summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
  summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
  summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
  summary_tab <- summary_tab %>%
    mutate(Significance = case_when(
      (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
      (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
      (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
      (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
      TRUE ~ "ns"
    ))
  summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
  ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
  write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
              quote = F, row.names = T, sep = "\t")
}



X_name <- "Lachnoclostridium"
M1_name <- "HostM_cg04794420"
M2_name <- "HostT_ERBB3"
Y_name <- "HostP_HER3"

Lachnoclostridium_simple_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                                  meds=c(M1_name,M2_name),
                                  med.type="serial",
                                  ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"


summary_tab <- Lachnoclostridium_simple_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")

summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)

X_name <- "Lachnoclostridium"
M1_name <- "HostM_cg07151117"
M2_name <- "HostT_ERBB3"
Y_name <- "HostP_HER3"

Lachnoclostridium_simple_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                                                  meds=c(M1_name,M2_name),
                                                  med.type="serial",
                                                  ci="boot", nsim=999, seed=2023)# or omit "mod.type", default is "2-way"

summary_tab <- Lachnoclostridium_simple_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")

summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)


#########################################

file_list<- ls()
rm(list=file_list[which(file_list != 'mediation_all_tab')])
project <- "TCGA-LUSC"
result.loc <- paste0("../result/SparseCCA_xena/",project,"/")


input_dirname <- paste0("../data/Preprocessed/xena_tcga/",project,"/case/")
filenames <- list.files(input_dirname)
cat("Read data from ",input_dirname," .....")


feat.rppa <- as.matrix(fread(paste0(input_dirname,grep("RPPA",filenames,value = TRUE))),rownames=1)
cat('Feature microbiome otu dimensions:', dim(feat.rppa), '\n')

# rppa_interest_list <-c("BCL2","IRS1","ERALPHA_pS118","ERALPHA") 
rppa_interest_list <-c("HER3","ERALPHA")

feat.rppa.interest <- as.data.frame(feat.rppa[,rppa_interest_list])
# feat.rppa.interest <- sapply(feat.rppa.interest, as.numeric)
rm(feat.rppa)

feat.otu.relt <- read.table(paste0(input_dirname,grep("Kraken",filenames,value = TRUE)))
colnames(feat.otu.relt) <- str_split_fixed(colnames(feat.otu.relt),'g__',2)[,2]
taxa_interest <- c("Lachnoclostridium","Sutterella")
feat.otu.relt.interest <- as.data.frame(feat.otu.relt[,taxa_interest])



cat('Feature microbiome otu dimensions:', dim(feat.otu.relt.interest), '\n')
rm(feat.otu.relt)

feat.gene.expr <- read.table(paste0(input_dirname,grep("Geneexpr",filenames,value = TRUE)))
gene_interest <- c("ESR1","ERBB3")
feat.gene.expr.interest <- as.data.frame(feat.gene.expr[,gene_interest])
rm(feat.gene.expr)

feat.meth <- as.matrix(fread(paste0(input_dirname,grep("CpG",filenames,value = TRUE))),rownames=1)
meth_interest <- c("cg04794420","cg07746998")
feat.meth.interest <- as.data.frame(feat.meth[,meth_interest])

idx <- Reduce(intersect, list(rownames(feat.otu.relt.interest),
                              rownames(feat.gene.expr.interest),
                              rownames(feat.rppa.interest),
                              rownames(feat.meth.interest)))

feat.otu.relt.interest <- feat.otu.relt.interest[idx,]
feat.gene.expr.interest <- feat.gene.expr.interest[idx,]      
feat.rppa.interest <- feat.rppa.interest[idx,]
feat.meth.interest <- feat.meth.interest[idx,]

colnames(feat.rppa.interest) <- paste0("HostP_",colnames(feat.rppa.interest))
colnames(feat.gene.expr.interest) <- paste0("HostT_",colnames(feat.gene.expr.interest))
colnames(feat.meth.interest) <- paste0("HostM_",colnames(feat.meth.interest))
# Omics_data <- data.frame(feat.otu.relt.interest,feat.gene.expr.interest,feat.rppa.interest)
# 
# X <-  feat.otu.relt.interest[,"Lachnoclostridium"]
# 
# Y <- feat.rppa.interest[,"HostP_ERALPHA"]
# M2 <- feat.gene.expr.interest[,"HostT_ESR1"]
# M1 <- feat.meth.interest[,"cg07746998"]

# 
LUSC_BG_deconvolution <- read.table(file="C:/Lab/project/pan-cancer-host-microbiome-associations/data/cell_type/LUSC_BG_deconvolution.t.tsv",sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)


# cell_type_interest <- c("Lymphocytes","Tregs","CD8_T_cells_PD1_high","CD8_T_cells_PD1_low","T_cells","Tregs","Macrophages_M1","Macrophages","CD8_T_cells","CD4_T_cells")
cell_interest <- LUSC_BG_deconvolution[idx,grep("_std$", colnames(LUSC_BG_deconvolution), invert = TRUE, value = TRUE)]
# 
# cell_interest <- LUSC_BG_deconvolution[idx,cell_type_interest]


# Meddata <- data.frame(X, M1,M2,Y_cell,Y)
Meddata <- data.frame(feat.otu.relt.interest,feat.meth.interest,feat.gene.expr.interest,feat.rppa.interest,cell_interest)

Meddata_scale <- apply(Meddata,2,scale)


for (i in 1:ncol(cell_interest)){
  cat("Analysis with ",colnames(cell_interest)[i],".\n")
  X_name <- "Lachnoclostridium"
  Y_name <- colnames(cell_interest)[i]
  M1_name <- "HostM_cg07746998"
  M2_name <- "HostT_ESR1"
  M3_name <- "HostP_ERALPHA"
  Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                                             meds=c(M1_name,M2_name,M3_name),
                                             med.type="serial",
                                             ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
  summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
  summary_tab$project <- project
  # print(summary_tab[8,4] < 0.05)
  summary_tab$Treat <- X_name
  summary_tab$Y <- Y_name
  summary_tab$Mediator <- paste(c(M1_name,M2_name,M3_name),collapse = ",")
  summary_tab$Path_simple <- rownames(summary_tab)
  summary_tab$Path <- rownames(summary_tab)
  # replace
  summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
  summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
  summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
  summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
  summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
  summary_tab <- summary_tab %>%
    mutate(Significance = case_when(
      (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
      (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
      (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
      (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
      TRUE ~ "ns"
    ))
  summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
  ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
  mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
  write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
              quote = F, row.names = T, sep = "\t")
}

X_name <- "Lachnoclostridium"
M1_name <- "HostM_cg07746998"
M2_name <- "HostT_ESR1"
Y_name <- "HostP_ERALPHA"


Lachnoclostridium_simple_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                                                  meds=c(M1_name,M2_name),
                                                  med.type="serial",
                                                  ci="boot", nsim=999, seed=24)# or omit "mod.type", default is "2-way"

summary_tab <- Lachnoclostridium_simple_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")

summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)


###########
# 
# Ruegeria_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
#                                   meds=c(M1_name,M2_name),
#                                   med.type="serial",
#                                   ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
# 
# 
# model <- bruceR::PROCESS(Meddata_scale, y="Macrophages_std", x="Lachnoclostridium",
#                          meds=c("HostM_cg04794420","HostT_ERBB3"),
#                          med.type="serial",
#                          ci="boot", nsim=999, seed=42)  # or omit "mod.type", default is "2-way"
# # 
# # model <- PROCESS_update(Meddata_scale, y="Macrophages_std", x="Lachnoclostridium",
# #                          meds=c("HostM_cg04794420","HostT_ERBB3"),
# #                          med.type="serial",
# #                          ci="boot", nsim=999, seed=42)  # or omit "mod.type", default is "2-way"
# # 
# # bruceR::PROCESS(Meddata_scale, y="HostT_ERBB3", x="Lachnoclostridium",
# #                 meds=c("HostM_cg04794420"),
# #                 med.type="serial",
# #                 ci="boot", nsim=999, seed=42) 
# 
# 
# model <- bruceR::PROCESS(Meddata_scale, y="HostP_HER3", x="Lachnoclostridium",
#                          meds=c("HostM_cg04794420","HostT_ERBB3"),
#                          med.type="serial",
#                          ci="boot", nsim=999, seed=42) # or omit "mod.type", default is "2-way"
# 
# X_name <- "Lachnoclostridium"
# Y_name <- "HostP_HER3"
# mediator_name <- paste("HostM_cg04794420","&","HostT_ERBB3")
# 
# summary_tab <- model$results[[1]]$lavaan.mediation
# # print(summary_tab[8,4] < 0.05)
# summary_tab$Path <- rownames(summary_tab)
# # summary_tab  <- summary_tab%>% mutate(signif = case_when(
# #   BootLLCI!=0 | BootULCI != 0 &pval <= 0.05 ~"*",
# #   BootLLCI!=0 | BootULCI != 0 &pval < 0.01 ~"**",
# #   BootLLCI!=0 | BootULCI != 0 &pval < 0.001 ~"***"
# # ))
# # 
# # ifelse(summary_tab$BootLLCI !=0, ,
# #        ifelse(sparseCCA_nodes_rppa_comp$feat_group == "taxa", "triangle",
# #               ifelse(sparseCCA_nodes_rppa_comp$feat_group == "meth", "circle",
# #                      ifelse(sparseCCA_nodes_rppa_comp$feat_group == "rppa", "diamond", NA))))
# 
# write.table(summary_tab, file = paste("../result/mediation/",X_name,"_affects_",Y_name, "_through_",mediator_name,".txt",sep = ""),
#             quote = F, row.names = T, sep = "\t")


####################
# rm(list = ls()) ## check to make sure no precomputed dataframe before clearing
file_list<- ls()
rm(list=file_list[which(file_list != 'mediation_all_tab')])

setwd("C:/Lab/project/pan-cancer-host-microbiome-associations/src")
project <- "TCGA-BRCA"
result.loc <- paste0("../result/SparseCCA_xena/",project,"/")



input_dirname <- paste0("../data/Preprocessed/xena_tcga/",project,"/case/")
filenames <- list.files(input_dirname)
cat("Read data from ",input_dirname," .....")


feat.rppa <- as.matrix(fread(paste0(input_dirname,grep("RPPA",filenames,value = TRUE))),rownames=1)
cat('Feature microbiome otu dimensions:', dim(feat.rppa), '\n')

# rppa_interest_list <-c("BCL2","IRS1","ERALPHA_pS118","ERALPHA") 
rppa_interest_list <-c("HER3","BCL2","IRS1","CYCLIND1","ERALPHA_pS118")

feat.rppa.interest <- as.data.frame(feat.rppa[,rppa_interest_list])
# feat.rppa.interest <- sapply(feat.rppa.interest, as.numeric)
rm(feat.rppa)

feat.otu.relt <- read.table(paste0(input_dirname,grep("Kraken",filenames,value = TRUE)))
colnames(feat.otu.relt) <- str_split_fixed(colnames(feat.otu.relt),'g__',2)[,2]
taxa_interest <- c("Lachnoclostridium","Ruegeria","Simplexvirus")
feat.otu.relt.interest <- as.data.frame(feat.otu.relt[,taxa_interest])



cat('Feature microbiome otu dimensions:', dim(feat.otu.relt.interest), '\n')
rm(feat.otu.relt)

feat.gene.expr <- read.table(paste0(input_dirname,grep("Geneexpr",filenames,value = TRUE)))
gene_interest <- c("ESR1","BCL2","IRS1","ESR1","CCND1")
feat.gene.expr.interest <- as.data.frame(feat.gene.expr[,gene_interest])
rm(feat.gene.expr)

feat.meth <- as.matrix(fread(paste0(input_dirname,grep("CpG",filenames,value = TRUE))),rownames=1)
meth_interest <- c("cg24740531","cg04751089","cg19411146","cg24900983","cg00953256","cg07781399","cg14513281","cg13608094","cg19964454")
feat.meth.interest <- as.data.frame(feat.meth[,meth_interest])
rm(feat.meth)
idx <- Reduce(intersect, list(rownames(feat.otu.relt.interest),
                              rownames(feat.gene.expr.interest),
                              rownames(feat.rppa.interest),
                              rownames(feat.meth.interest)))

feat.otu.relt.interest <- feat.otu.relt.interest[idx,]
feat.gene.expr.interest <- feat.gene.expr.interest[idx,]      
feat.rppa.interest <- feat.rppa.interest[idx,]
feat.meth.interest <- feat.meth.interest[idx,]

colnames(feat.rppa.interest) <- paste0("HostP_",colnames(feat.rppa.interest))
colnames(feat.gene.expr.interest) <- paste0("HostT_",colnames(feat.gene.expr.interest))
colnames(feat.meth.interest) <- paste0("HostM_",colnames(feat.meth.interest))

# 
BRCA_BG_deconvolution <- read.table(file="C:/Lab/project/pan-cancer-host-microbiome-associations/data/cell_type/BRCA_BG_deconvolution.t.tsv",sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)


# cell_type_interest <- c("Lymphocytes","Tregs","CD8_T_cells_PD1_high","CD8_T_cells_PD1_low","T_cells","Tregs","Macrophages_M1","Macrophages","CD8_T_cells","CD4_T_cells")
cell_interest <- BRCA_BG_deconvolution[idx,grep("_std$", colnames(BRCA_BG_deconvolution), invert = TRUE, value = TRUE)]
# 


Meddata <- data.frame(feat.otu.relt.interest,feat.meth.interest,feat.gene.expr.interest,feat.rppa.interest,cell_interest)

Meddata_scale <- apply(Meddata,2,scale)

# 
# for (i in 1:ncol(cell_interest)){
#   cat("Analysis with ",colnames(cell_interest)[i],".\n")
#   X_name <- "Ruegeria"
#   Y_name <- colnames(cell_interest)[i]
#   M1_name <- "HostM_cg24740531"
#   M2_name <- "HostT_BCL2"
#   M3_name <- "HostP_BCL2"
#   Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
#                                              meds=c(M1_name,M2_name,M3_name),
#                                              med.type="serial",
#                                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
#   summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
#   summary_tab$project <- project
#   # print(summary_tab[8,4] < 0.05)
#   summary_tab$Treat <- X_name
#   summary_tab$Y <- Y_name
#   summary_tab$Mediator <- paste(c(M1_name,M2_name,M3_name),collapse = ",")
#   summary_tab$Path_simple <- rownames(summary_tab)
#   summary_tab$Path <- rownames(summary_tab)
#   # replace
#   summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
#   summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
#   summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
#   summary_tab <- summary_tab %>%
#     mutate(Significance = case_when(
#       (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
#       TRUE ~ "ns"
#     ))
#   summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
#   ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
#   mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
#   write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
#               quote = F, row.names = T, sep = "\t")
# }
# 
# for (i in 1:ncol(cell_interest)){
#   cat("Analysis with ",colnames(cell_interest)[i],".\n")
#   X_name <- "Simplexvirus"
#   Y_name <- colnames(cell_interest)[i]
#   M1_name <- "HostM_cg04751089"
#   M2_name <- "HostT_IRS1"
#   M3_name <- "HostP_IRS1"
#   Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
#                                              meds=c(M1_name,M2_name,M3_name),
#                                              med.type="serial",
#                                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
#   summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
#   summary_tab$project <- project
#   # print(summary_tab[8,4] < 0.05)
#   summary_tab$Treat <- X_name
#   summary_tab$Y <- Y_name
#   summary_tab$Mediator <- paste(c(M1_name,M2_name,M3_name),collapse = ",")
#   summary_tab$Path_simple <- rownames(summary_tab)
#   summary_tab$Path <- rownames(summary_tab)
#   # replace
#   summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
#   summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
#   summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
#   summary_tab <- summary_tab %>%
#     mutate(Significance = case_when(
#       (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
#       TRUE ~ "ns"
#     ))
#   summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
#   ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
#   mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
#   write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
#               quote = F, row.names = T, sep = "\t")
# }
# 
# 
# for (i in 1:ncol(cell_interest)){
#   cat("Analysis with ",colnames(cell_interest)[i],".\n")
#   X_name <- "Simplexvirus"
#   Y_name <- colnames(cell_interest)[i]
#   M1_name <- "HostM_cg19411146"
#   M2_name <- "HostT_ESR1"
#   M3_name <- "HostP_ERALPHA_pS118"
#   Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
#                                              meds=c(M1_name,M2_name,M3_name),
#                                              med.type="serial",
#                                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
#   summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
#   summary_tab$project <- project
#   # print(summary_tab[8,4] < 0.05)
#   summary_tab$Treat <- X_name
#   summary_tab$Y <- Y_name
#   summary_tab$Mediator <- paste(c(M1_name,M2_name,M3_name),collapse = ",")
#   summary_tab$Path_simple <- rownames(summary_tab)
#   summary_tab$Path <- rownames(summary_tab)
#   # replace
#   summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
#   summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
#   summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
#   summary_tab <- summary_tab %>%
#     mutate(Significance = case_when(
#       (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
#       TRUE ~ "ns"
#     ))
#   summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
#   ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
#   mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
#   write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
#               quote = F, row.names = T, sep = "\t")
# }
# 
# for (i in 1:ncol(cell_interest)){
#   cat("Analysis with ",colnames(cell_interest)[i],".\n")
#   X_name <- "Simplexvirus"
#   Y_name <- colnames(cell_interest)[i]
#   M1_name <- "HostM_cg24900983"
#   M2_name <- "HostT_ESR1"
#   M3_name <- "HostP_ERALPHA_pS118"
#   Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
#                                              meds=c(M1_name,M2_name,M3_name),
#                                              med.type="serial",
#                                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
#   summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
#   summary_tab$project <- project
#   # print(summary_tab[8,4] < 0.05)
#   summary_tab$Treat <- X_name
#   summary_tab$Y <- Y_name
#   summary_tab$Mediator <- paste(c(M1_name,M2_name,M3_name),collapse = ",")
#   summary_tab$Path_simple <- rownames(summary_tab)
#   summary_tab$Path <- rownames(summary_tab)
#   # replace
#   summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
#   summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
#   summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
#   summary_tab <- summary_tab %>%
#     mutate(Significance = case_when(
#       (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
#       TRUE ~ "ns"
#     ))
#   summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
#   ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
#   mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
#   write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
#               quote = F, row.names = T, sep = "\t")
# }
# 
# 
# for (i in 1:ncol(cell_interest)){
#   cat("Analysis with ",colnames(cell_interest)[i],".\n")
#   X_name <- "Simplexvirus"
#   Y_name <- colnames(cell_interest)[i]
#   M1_name <- "HostM_cg00953256"
#   M2_name <- "HostT_CCND1"
#   M3_name <- "HostP_CYCLIND1"
#   Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
#                                              meds=c(M1_name,M2_name,M3_name),
#                                              med.type="serial",
#                                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
#   summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
#   summary_tab$project <- project
#   # print(summary_tab[8,4] < 0.05)
#   summary_tab$Treat <- X_name
#   summary_tab$Y <- Y_name
#   summary_tab$Mediator <- paste(c(M1_name,M2_name,M3_name),collapse = ",")
#   summary_tab$Path_simple <- rownames(summary_tab)
#   summary_tab$Path <- rownames(summary_tab)
#   # replace
#   summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
#   summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
#   summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
#   summary_tab <- summary_tab %>%
#     mutate(Significance = case_when(
#       (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
#       TRUE ~ "ns"
#     ))
#   summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
#   ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
#   mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
#   write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
#               quote = F, row.names = T, sep = "\t")
# }
# 
# for (i in 1:ncol(cell_interest)){
#   cat("Analysis with ",colnames(cell_interest)[i],".\n")
#   X_name <- "Simplexvirus"
#   Y_name <- colnames(cell_interest)[i]
#   M1_name <- "HostM_cg00953256"
#   M2_name <- "HostT_CCND1"
#   M3_name <- "HostP_CYCLIND1"
#   Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
#                                              meds=c(M1_name,M2_name,M3_name),
#                                              med.type="serial",
#                                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
#   summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
#   summary_tab$project <- project
#   # print(summary_tab[8,4] < 0.05)
#   summary_tab$Treat <- X_name
#   summary_tab$Y <- Y_name
#   summary_tab$Mediator <- paste(c(M1_name,M2_name,M3_name),collapse = ",")
#   summary_tab$Path_simple <- rownames(summary_tab)
#   summary_tab$Path <- rownames(summary_tab)
#   # replace
#   summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
#   summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
#   summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
#   summary_tab <- summary_tab %>%
#     mutate(Significance = case_when(
#       (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
#       TRUE ~ "ns"
#     ))
#   summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
#   ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
#   mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
#   write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
#               quote = F, row.names = T, sep = "\t")
# }
# 
# for (i in 1:ncol(cell_interest)){
#   cat("Analysis with ",colnames(cell_interest)[i],".\n")
#   X_name <- "Simplexvirus"
#   Y_name <- colnames(cell_interest)[i]
#   M1_name <- "HostM_cg07781399"
#   M2_name <- "HostT_CCND1"
#   M3_name <- "HostP_CYCLIND1"
#   Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
#                                              meds=c(M1_name,M2_name,M3_name),
#                                              med.type="serial",
#                                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
#   summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
#   summary_tab$project <- project
#   # print(summary_tab[8,4] < 0.05)
#   summary_tab$Treat <- X_name
#   summary_tab$Y <- Y_name
#   summary_tab$Mediator <- paste(c(M1_name,M2_name,M3_name),collapse = ",")
#   summary_tab$Path_simple <- rownames(summary_tab)
#   summary_tab$Path <- rownames(summary_tab)
#   # replace
#   summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
#   summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
#   summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
#   summary_tab <- summary_tab %>%
#     mutate(Significance = case_when(
#       (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
#       TRUE ~ "ns"
#     ))
#   summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
#   ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
#   mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
#   write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
#               quote = F, row.names = T, sep = "\t")
# }
# 
# 
# for (i in 1:ncol(cell_interest)){
#   cat("Analysis with ",colnames(cell_interest)[i],".\n")
#   X_name <- "Simplexvirus"
#   Y_name <- colnames(cell_interest)[i]
#   M1_name <- "HostM_cg14513281"
#   M2_name <- "HostT_CCND1"
#   M3_name <- "HostP_CYCLIND1"
#   Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
#                                              meds=c(M1_name,M2_name,M3_name),
#                                              med.type="serial",
#                                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
#   summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
#   summary_tab$project <- project
#   # print(summary_tab[8,4] < 0.05)
#   summary_tab$Treat <- X_name
#   summary_tab$Y <- Y_name
#   summary_tab$Mediator <- paste(c(M1_name,M2_name,M3_name),collapse = ",")
#   summary_tab$Path_simple <- rownames(summary_tab)
#   summary_tab$Path <- rownames(summary_tab)
#   # replace
#   summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
#   summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
#   summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
#   summary_tab <- summary_tab %>%
#     mutate(Significance = case_when(
#       (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
#       TRUE ~ "ns"
#     ))
#   summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
#   ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
#   mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
#   write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
#               quote = F, row.names = T, sep = "\t")
# }
# 
# for (i in 1:ncol(cell_interest)){
#   cat("Analysis with ",colnames(cell_interest)[i],".\n")
#   X_name <- "Simplexvirus"
#   Y_name <- colnames(cell_interest)[i]
#   M1_name <- "HostM_cg13608094"
#   M2_name <- "HostT_CCND1"
#   M3_name <- "HostP_CYCLIND1"
#   Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
#                                              meds=c(M1_name,M2_name,M3_name),
#                                              med.type="serial",
#                                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
#   summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
#   summary_tab$project <- project
#   # print(summary_tab[8,4] < 0.05)
#   summary_tab$Treat <- X_name
#   summary_tab$Y <- Y_name
#   summary_tab$Mediator <- paste(c(M1_name,M2_name,M3_name),collapse = ",")
#   summary_tab$Path_simple <- rownames(summary_tab)
#   summary_tab$Path <- rownames(summary_tab)
#   # replace
#   summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
#   summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
#   summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
#   summary_tab <- summary_tab %>%
#     mutate(Significance = case_when(
#       (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
#       TRUE ~ "ns"
#     ))
#   summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
#   ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
#   mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
#   write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
#               quote = F, row.names = T, sep = "\t")
# }
# 
# for (i in 1:ncol(cell_interest)){
#   cat("Analysis with ",colnames(cell_interest)[i],".\n")
#   X_name <- "Simplexvirus"
#   Y_name <- colnames(cell_interest)[i]
#   M1_name <- "HostM_cg19964454"
#   M2_name <- "HostT_CCND1"
#   M3_name <- "HostP_CYCLIND1"
#   Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
#                                              meds=c(M1_name,M2_name,M3_name),
#                                              med.type="serial",
#                                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
#   summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
#   summary_tab$project <- project
#   # print(summary_tab[8,4] < 0.05)
#   summary_tab$Treat <- X_name
#   summary_tab$Y <- Y_name
#   summary_tab$Mediator <- paste(c(M1_name,M2_name,M3_name),collapse = ",")
#   summary_tab$Path_simple <- rownames(summary_tab)
#   summary_tab$Path <- rownames(summary_tab)
#   # replace
#   summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
#   summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
#   summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
#   summary_tab <- summary_tab %>%
#     mutate(Significance = case_when(
#       (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
#       TRUE ~ "ns"
#     ))
#   summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
#   ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
#   mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
#   write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
#               quote = F, row.names = T, sep = "\t")
# }
# 





X_name <- "Ruegeria"
Y_name <- "HostP_BCL2"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg24740531"
M2_name <- "HostT_BCL2"

Ruegeria_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                              meds=c(M1_name,M2_name),
                              med.type="serial",
                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"


summary_tab <- Ruegeria_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")

summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)



X_name <- "Simplexvirus"
Y_name <- "HostP_IRS1"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg04751089"
M2_name <- "HostT_IRS1"

Simplexvirus_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                                  meds=c(M1_name,M2_name),
                                  med.type="serial",
                                  ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"

summary_tab <- Simplexvirus_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")

summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)

# 
# 
# Simplexvirus_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
#                                       meds=c(M1_name,M2_name),
#                                       med.type="serial",
#                                       ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
# 
# 
# summary_tab <- Simplexvirus_model$results[[1]]$lavaan.mediation
# # print(summary_tab[8,4] < 0.05)
# summary_tab$Path <- rownames(summary_tab)
# # replace
# summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
# summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
# summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
# summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
# summary_tab <- summary_tab %>%
#   mutate(Significance = case_when(
#     (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
#     (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
#     (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
#     (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
#     TRUE ~ "ns"
#   ))
# summary_tab$project <- project
# write.table(summary_tab, file = paste("../result/mediation/",project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
#             quote = F, row.names = T, sep = "\t")
###
X_name <- "Simplexvirus"
Y_name <- "HostP_ERALPHA_pS118"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg24900983"
M2_name <- "HostT_ESR1"

Simplexvirus_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                                      meds=c(M1_name,M2_name),
                                      med.type="serial",
                                      ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"

summary_tab <- Simplexvirus_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")

summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)

###

X_name <- "Simplexvirus"
Y_name <- "HostP_CYCLIND1"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg00953256"
M2_name <- "HostT_CCND1"

Simplexvirus_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                                      meds=c(M1_name,M2_name),
                                      med.type="serial",
                                      ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"

summary_tab <- Simplexvirus_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")

summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
###

X_name <- "Simplexvirus"
Y_name <- "HostP_CYCLIND1"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg07781399"
M2_name <- "HostT_CCND1"

Simplexvirus_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                                      meds=c(M1_name,M2_name),
                                      med.type="serial",
                                      ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"

summary_tab <- Simplexvirus_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")

summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)

####

X_name <- "Simplexvirus"
Y_name <- "HostP_CYCLIND1"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg13608094"
M2_name <- "HostT_CCND1"

Simplexvirus_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                                      meds=c(M1_name,M2_name),
                                      med.type="serial",
                                      ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"

summary_tab <- Simplexvirus_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")

summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)



mediation_all_tab$Y_name <- str_extract(mediation_all_tab$biolink, "(?<=Y:\\[)[^\\]]+")
mediation_all_tab$X_name <- str_extract(mediation_all_tab$biolink, "(?<=Treat:\\[)[^\\]]+")
mediation_all_tab$M1_name <- str_extract(mediation_all_tab$biolink, "(?<=Mediator:\\[)[^|]+")
mediation_all_tab$M2_name <- str_extract(mediation_all_tab$biolink, "(?<=\\|)[^\\]]+")


#########
X_name <- "Simplexvirus"
Y_name <- "HostP_CYCLIND1"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg19964454"
M2_name <- "HostT_CCND1"

Simplexvirus_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                                      meds=c(M1_name,M2_name),
                                      med.type="serial",
                                      ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"

summary_tab <- Simplexvirus_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")
summary_tab$Y_name <- Y_name
summary_tab$X_name <- X_name
summary_tab$M1_name <- M1_name
summary_tab$M2_name <- M2_name
summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)




###
X_name <- "Simplexvirus"
Y_name <- "HostP_CYCLIND1"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg14513281"
M2_name <- "HostT_CCND1"

Simplexvirus_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                                      meds=c(M1_name,M2_name),
                                      med.type="serial",
                                      ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
summary_tab <- Simplexvirus_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")
summary_tab$Y_name <- Y_name
summary_tab$X_name <- X_name
summary_tab$M1_name <- M1_name
summary_tab$M2_name <- M2_name
summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)

#####





############
### LGG

# rm(list = ls()) ## check to make sure no precomputed dataframe before clearing
file_list<- ls()
rm(list=file_list[which(file_list != 'mediation_all_tab')])


setwd("C:/Lab/project/pan-cancer-host-microbiome-associations/src")
project <- "TCGA-LGG"
result.loc <- paste0("../result/SparseCCA_xena/",project,"/")



input_dirname <- paste0("../data/Preprocessed/xena_tcga/",project,"/case/")
filenames <- list.files(input_dirname)
cat("Read data from ",input_dirname," .....")


feat.rppa <- as.matrix(fread(paste0(input_dirname,grep("RPPA",filenames,value = TRUE))),rownames=1)
cat('Feature microbiome otu dimensions:', dim(feat.rppa), '\n')

# rppa_interest_list <-c("BCL2","IRS1","ERALPHA_pS118","ERALPHA") 
rppa_interest_list <-c("HER3","YAP_pS127","SYK")

feat.rppa.interest <- as.data.frame(feat.rppa[,rppa_interest_list])
# feat.rppa.interest <- sapply(feat.rppa.interest, as.numeric)
rm(feat.rppa)

feat.otu.relt <- read.table(paste0(input_dirname,grep("Kraken",filenames,value = TRUE)))
colnames(feat.otu.relt) <- str_split_fixed(colnames(feat.otu.relt),'g__',2)[,2]
taxa_interest <- c("Lachnoclostridium","Blautia")
feat.otu.relt.interest <- as.data.frame(feat.otu.relt[,taxa_interest])



cat('Feature microbiome otu dimensions:', dim(feat.otu.relt.interest), '\n')
rm(feat.otu.relt)

feat.gene.expr <- read.table(paste0(input_dirname,grep("Geneexpr",filenames,value = TRUE)))
gene_interest <- c("ESR1","SYK","YAP1")
feat.gene.expr.interest <- as.data.frame(feat.gene.expr[,gene_interest])
rm(feat.gene.expr)

feat.meth <- as.matrix(fread(paste0(input_dirname,grep("CpG",filenames,value = TRUE))),rownames=1)
meth_interest <- c("cg01442799","cg15894431","cg26369667")
feat.meth.interest <- as.data.frame(feat.meth[,meth_interest])
rm(feat.meth)
idx <- Reduce(intersect, list(rownames(feat.otu.relt.interest),
                              rownames(feat.gene.expr.interest),
                              rownames(feat.rppa.interest),
                              rownames(feat.meth.interest)))

feat.otu.relt.interest <- feat.otu.relt.interest[idx,]
feat.gene.expr.interest <- feat.gene.expr.interest[idx,]      
feat.rppa.interest <- feat.rppa.interest[idx,]
feat.meth.interest <- feat.meth.interest[idx,]

colnames(feat.rppa.interest) <- paste0("HostP_",colnames(feat.rppa.interest))
colnames(feat.gene.expr.interest) <- paste0("HostT_",colnames(feat.gene.expr.interest))
colnames(feat.meth.interest) <- paste0("HostM_",colnames(feat.meth.interest))

LGG_BG_deconvolution <- read.table(file="C:/Lab/project/pan-cancer-host-microbiome-associations/data/cell_type/LGG_BG_deconvolution.t.tsv",sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)


cell_interest <- LGG_BG_deconvolution[idx,grep("_std$", colnames(LGG_BG_deconvolution), invert = TRUE, value = TRUE)]




Meddata <- data.frame(feat.otu.relt.interest,feat.meth.interest,feat.gene.expr.interest,feat.rppa.interest,cell_interest)

Meddata_scale <- apply(Meddata,2,scale)


LGG_model <- bruceR::PROCESS(Meddata_scale, y="HostP_SYK", x="Blautia",
                                  meds=c("HostM_cg01442799","HostT_SYK"),
                                  med.type="serial",
                                  ci="boot", nsim=999, seed=42) # or omit "mod.type", default is "2-way"


X_name <- "Blautia"
Y_name <- "HostP_SYK"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg01442799"
M2_name <- "HostT_SYK"

summary_tab <- LGG_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")
summary_tab$Y_name <- Y_name
summary_tab$X_name <- X_name
summary_tab$M1_name <- M1_name
summary_tab$M2_name <- M2_name
summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)


#############

X_name <- "Blautia"
Y_name <- "HostP_YAP_pS127"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg15894431"
M2_name <- "HostT_YAP1"

LGG_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                             meds=c(M1_name,M2_name),
                             med.type="serial",
                             ci="boot", nsim=999, seed=2023)# or omit "mod.type", default is "2-way"
# LGG_model <- bruceR::PROCESS(Meddata_scale, y=M2_name, x=X_name,
#                              meds=c(M1_name),
#                              med.type="serial",
#                              ci="boot", nsim=999, seed=42) 

summary_tab <- LGG_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")
summary_tab$Y_name <- Y_name
summary_tab$X_name <- X_name
summary_tab$M1_name <- M1_name
summary_tab$M2_name <- M2_name
summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)

X_name <- "Blautia"
Y_name <- "HostP_YAP_pS127"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg26369667"
M2_name <- "HostT_YAP1"

LGG_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                             meds=c(M1_name,M2_name),
                             med.type="serial",
                             ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
# LGG_model <- bruceR::PROCESS(Meddata_scale, y=M2_name, x=X_name,
#                              meds=c(M1_name),
#                              med.type="serial",
#                              ci="boot", nsim=999, seed=42) 

summary_tab <- LGG_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")
summary_tab$Y_name <- Y_name
summary_tab$X_name <- X_name
summary_tab$M1_name <- M1_name
summary_tab$M2_name <- M2_name
summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)

###############
### COAD
###############
file_list<- ls()
rm(list=file_list[which(file_list != 'mediation_all_tab')])
setwd("C:/Lab/project/pan-cancer-host-microbiome-associations/src")
project <- "TCGA-COAD"
result.loc <- paste0("../result/SparseCCA_xena/",project,"/")



input_dirname <- paste0("../data/Preprocessed/xena_tcga/",project,"/case/")
filenames <- list.files(input_dirname)
cat("Read data from ",input_dirname," .....")


feat.rppa <- as.matrix(fread(paste0(input_dirname,grep("RPPA",filenames,value = TRUE))),rownames=1)
cat('Feature microbiome otu dimensions:', dim(feat.rppa), '\n')

# rppa_interest_list <-c("BCL2","IRS1","ERALPHA_pS118","ERALPHA") 
rppa_interest_list <-c("HER3","YAP_pS127","DUSP4")

feat.rppa.interest <- as.data.frame(feat.rppa[,rppa_interest_list])
# feat.rppa.interest <- sapply(feat.rppa.interest, as.numeric)
rm(feat.rppa)

feat.otu.relt <- read.table(paste0(input_dirname,grep("Kraken",filenames,value = TRUE)))
colnames(feat.otu.relt) <- str_split_fixed(colnames(feat.otu.relt),'g__',2)[,2]
taxa_interest <- c("Lachnoclostridium","Luteibacter")
feat.otu.relt.interest <- as.data.frame(feat.otu.relt[,taxa_interest])



cat('Feature microbiome otu dimensions:', dim(feat.otu.relt.interest), '\n')
rm(feat.otu.relt)

feat.gene.expr <- read.table(paste0(input_dirname,grep("Geneexpr",filenames,value = TRUE)))
gene_interest <- c("ESR1","SYK","DUSP4")
feat.gene.expr.interest <- as.data.frame(feat.gene.expr[,gene_interest])
rm(feat.gene.expr)

feat.meth <- as.matrix(fread(paste0(input_dirname,grep("CpG",filenames,value = TRUE))),rownames=1)
meth_interest <- c("cg18432895","cg03693434")
feat.meth.interest <- as.data.frame(feat.meth[,meth_interest])
rm(feat.meth)
idx <- Reduce(intersect, list(rownames(feat.otu.relt.interest),
                              rownames(feat.gene.expr.interest),
                              rownames(feat.rppa.interest),
                              rownames(feat.meth.interest)))

feat.otu.relt.interest <- feat.otu.relt.interest[idx,]
feat.gene.expr.interest <- feat.gene.expr.interest[idx,]      
feat.rppa.interest <- feat.rppa.interest[idx,]
feat.meth.interest <- feat.meth.interest[idx,]

colnames(feat.rppa.interest) <- paste0("HostP_",colnames(feat.rppa.interest))
colnames(feat.gene.expr.interest) <- paste0("HostT_",colnames(feat.gene.expr.interest))
colnames(feat.meth.interest) <- paste0("HostM_",colnames(feat.meth.interest))

Meddata <- data.frame(feat.otu.relt.interest,feat.meth.interest,feat.gene.expr.interest,feat.rppa.interest)

Meddata_scale <- apply(Meddata,2,scale)


X_name <- "Luteibacter"
Y_name <- "HostP_DUSP4"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg18432895"
M2_name <- "HostT_DUSP4"

COAD_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                             meds=c(M1_name,M2_name),
                             med.type="serial",
                             ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
# LGG_model <- bruceR::PROCESS(Meddata_scale, y=M2_name, x=X_name,
#                              meds=c(M1_name),
#                              med.type="serial",
#                              ci="boot", nsim=999, seed=42) 

summary_tab <- COAD_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")
summary_tab$Y_name <- Y_name
summary_tab$X_name <- X_name
summary_tab$M1_name <- M1_name
summary_tab$M2_name <- M2_name
summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)

###############
### STAD
###############
# rm(list = ls()) ## check to make sure no precomputed dataframe before clearing
file_list<- ls()
rm(list=file_list[which(file_list != 'mediation_all_tab')])

setwd("C:/Lab/project/pan-cancer-host-microbiome-associations/src")
project <- "TCGA-STAD"
result.loc <- paste0("../result/SparseCCA_xena/",project,"/")



input_dirname <- paste0("../data/Preprocessed/xena_tcga/",project,"/case/")
filenames <- list.files(input_dirname)
cat("Read data from ",input_dirname," .....")


feat.rppa <- as.matrix(fread(paste0(input_dirname,grep("RPPA",filenames,value = TRUE))),rownames=1)
cat('Feature microbiome otu dimensions:', dim(feat.rppa), '\n')

# rppa_interest_list <-c("BCL2","IRS1","ERALPHA_pS118","ERALPHA") 
rppa_interest_list <-c("HER3","TRANSGLUTAMINASE")

feat.rppa.interest <- as.data.frame(feat.rppa[,rppa_interest_list])
# feat.rppa.interest <- sapply(feat.rppa.interest, as.numeric)
rm(feat.rppa)

feat.otu.relt <- read.table(paste0(input_dirname,grep("Kraken",filenames,value = TRUE)))
colnames(feat.otu.relt) <- str_split_fixed(colnames(feat.otu.relt),'g__',2)[,2]
taxa_interest <- c("Lachnoclostridium","Lymphocryptovirus")
feat.otu.relt.interest <- as.data.frame(feat.otu.relt[,taxa_interest])



cat('Feature microbiome otu dimensions:', dim(feat.otu.relt.interest), '\n')
rm(feat.otu.relt)

feat.gene.expr <- read.table(paste0(input_dirname,grep("Geneexpr",filenames,value = TRUE)))
gene_interest <- c("ESR1","TGM2")
feat.gene.expr.interest <- as.data.frame(feat.gene.expr[,gene_interest])
rm(feat.gene.expr)

feat.meth <- as.matrix(fread(paste0(input_dirname,grep("CpG",filenames,value = TRUE))),rownames=1)
meth_interest <- c("cg18432895","cg08545268")
feat.meth.interest <- as.data.frame(feat.meth[,meth_interest])
rm(feat.meth)
idx <- Reduce(intersect, list(rownames(feat.otu.relt.interest),
                              rownames(feat.gene.expr.interest),
                              rownames(feat.rppa.interest),
                              rownames(feat.meth.interest)))

feat.otu.relt.interest <- feat.otu.relt.interest[idx,]
feat.gene.expr.interest <- feat.gene.expr.interest[idx,]      
feat.rppa.interest <- feat.rppa.interest[idx,]
feat.meth.interest <- feat.meth.interest[idx,]

colnames(feat.rppa.interest) <- paste0("HostP_",colnames(feat.rppa.interest))
colnames(feat.gene.expr.interest) <- paste0("HostT_",colnames(feat.gene.expr.interest))
colnames(feat.meth.interest) <- paste0("HostM_",colnames(feat.meth.interest))



STAD_BG_deconvolution <- read.table(file="C:/Lab/project/pan-cancer-host-microbiome-associations/data/cell_type/STAD_BG_deconvolution.t.tsv",sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)


cell_interest <- STAD_BG_deconvolution[idx,grep("_std$", colnames(STAD_BG_deconvolution), invert = TRUE, value = TRUE)]




Meddata <- data.frame(feat.otu.relt.interest,feat.meth.interest,feat.gene.expr.interest,feat.rppa.interest,cell_interest)

# 
# Meddata <- data.frame(feat.otu.relt.interest,feat.meth.interest,feat.gene.expr.interest,feat.rppa.interest)

Meddata_scale <- apply(Meddata,2,scale)
# 
# 
# for (i in 1:ncol(cell_interest)){
#   cat("Analysis with ",colnames(cell_interest)[i],".\n")
#   X_name <- "Lymphocryptovirus"
#   Y_name <- colnames(cell_interest)[i]
#   M1_name <- "HostM_cg08545268"
#   M2_name <- "HostT_TGM2"
#   M3_name <- "HostP_TRANSGLUTAMINASE"
#   Lachnoclostridium_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
#                                              meds=c(M1_name,M2_name,M3_name),
#                                              med.type="serial",
#                                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
#   summary_tab <- Lachnoclostridium_model$results[[1]]$lavaan.mediation
#   summary_tab$project <- project
#   # print(summary_tab[8,4] < 0.05)
#   summary_tab$Treat <- X_name
#   summary_tab$Y <- Y_name
#   summary_tab$Mediator <- paste(c(M1_name,M2_name,M3_name),collapse = ",")
#   summary_tab$Path_simple <- rownames(summary_tab)
#   summary_tab$Path <- rownames(summary_tab)
#   # replace
#   summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
#   summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
#   summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)
#   summary_tab$Path <- gsub("M3", M3_name, summary_tab$Path)
#   summary_tab <- summary_tab %>%
#     mutate(Significance = case_when(
#       (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
#       (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
#       TRUE ~ "ns"
#     ))
#   summary_tab_filepath <- paste0("../result/mediation/",project,"/cell_type/")
#   ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
#   mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
#   write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,"+",M3_name,".txt",sep = ""),
#               quote = F, row.names = T, sep = "\t")
# }

X_name <- "Lymphocryptovirus"
Y_name <- "HostP_TRANSGLUTAMINASE"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg08545268"
M2_name <- "HostT_TGM2"

STAD_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                              meds=c(M1_name,M2_name),
                              med.type="serial",
                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
# LGG_model <- bruceR::PROCESS(Meddata_scale, y=M2_name, x=X_name,
#                              meds=c(M1_name),
#                              med.type="serial",
#                              ci="boot", nsim=999, seed=42) 


summary_tab <- STAD_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")
summary_tab$Y_name <- Y_name
summary_tab$X_name <- X_name
summary_tab$M1_name <- M1_name
summary_tab$M2_name <- M2_name
summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)
###############
### PRAD
###############
file_list<- ls()
rm(list=file_list[which(file_list != 'mediation_all_tab')])

setwd("C:/Lab/project/pan-cancer-host-microbiome-associations/src")
project <- "TCGA-PRAD"
result.loc <- paste0("../result/SparseCCA_xena/",project,"/")



input_dirname <- paste0("../data/Preprocessed/xena_tcga/",project,"/case/")
filenames <- list.files(input_dirname)
cat("Read data from ",input_dirname," .....")


feat.rppa <- as.matrix(fread(paste0(input_dirname,grep("RPPA",filenames,value = TRUE))),rownames=1)
cat('Feature microbiome otu dimensions:', dim(feat.rppa), '\n')

# rppa_interest_list <-c("BCL2","IRS1","ERALPHA_pS118","ERALPHA") 
rppa_interest_list <-c("HER3","FASN")

feat.rppa.interest <- as.data.frame(feat.rppa[,rppa_interest_list])
# feat.rppa.interest <- sapply(feat.rppa.interest, as.numeric)
rm(feat.rppa)

feat.otu.relt <- read.table(paste0(input_dirname,grep("Kraken",filenames,value = TRUE)))
colnames(feat.otu.relt) <- str_split_fixed(colnames(feat.otu.relt),'g__',2)[,2]
taxa_interest <- c("Lachnoclostridium","Nitrospira")
feat.otu.relt.interest <- as.data.frame(feat.otu.relt[,taxa_interest])



cat('Feature microbiome otu dimensions:', dim(feat.otu.relt.interest), '\n')
rm(feat.otu.relt)

feat.gene.expr <- read.table(paste0(input_dirname,grep("Geneexpr",filenames,value = TRUE)))
gene_interest <- c("ESR1","FASN")
feat.gene.expr.interest <- as.data.frame(feat.gene.expr[,gene_interest])
rm(feat.gene.expr)

feat.meth <- as.matrix(fread(paste0(input_dirname,grep("CpG",filenames,value = TRUE))),rownames=1)
meth_interest <- c("cg03693434","cg04029738","cg16699850","cg02608019")
feat.meth.interest <- as.data.frame(feat.meth[,meth_interest])
rm(feat.meth)
idx <- Reduce(intersect, list(rownames(feat.otu.relt.interest),
                              rownames(feat.gene.expr.interest),
                              rownames(feat.rppa.interest),
                              rownames(feat.meth.interest)))

feat.otu.relt.interest <- feat.otu.relt.interest[idx,]
feat.gene.expr.interest <- feat.gene.expr.interest[idx,]      
feat.rppa.interest <- feat.rppa.interest[idx,]
feat.meth.interest <- feat.meth.interest[idx,]

colnames(feat.rppa.interest) <- paste0("HostP_",colnames(feat.rppa.interest))
colnames(feat.gene.expr.interest) <- paste0("HostT_",colnames(feat.gene.expr.interest))
colnames(feat.meth.interest) <- paste0("HostM_",colnames(feat.meth.interest))

Meddata <- data.frame(feat.otu.relt.interest,feat.meth.interest,feat.gene.expr.interest,feat.rppa.interest)

Meddata_scale <- apply(Meddata,2,scale)


X_name <- "Nitrospira"
Y_name <- "HostP_FASN"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg03693434"
M2_name <- "HostT_FASN"

PRAD_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                              meds=c(M1_name,M2_name),
                              med.type="serial",
                              ci="boot", nsim=999, seed=2023)# or omit "mod.type", default is "2-way"

summary_tab <- PRAD_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")
summary_tab$Y_name <- Y_name
summary_tab$X_name <- X_name
summary_tab$M1_name <- M1_name
summary_tab$M2_name <- M2_name
summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)


X_name <- "Nitrospira"
Y_name <- "HostP_FASN"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg04029738"
M2_name <- "HostT_FASN"

PRAD_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                              meds=c(M1_name,M2_name),
                              med.type="serial",
                              ci="boot", nsim=999, seed=2023)# or omit "mod.type", default is "2-way"

summary_tab <- PRAD_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")
summary_tab$Y_name <- Y_name
summary_tab$X_name <- X_name
summary_tab$M1_name <- M1_name
summary_tab$M2_name <- M2_name
summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)



X_name <- "Nitrospira"
Y_name <- "HostP_FASN"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg16699850"
M2_name <- "HostT_FASN"

PRAD_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                              meds=c(M1_name,M2_name),
                              med.type="serial",
                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"

summary_tab <- PRAD_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")
summary_tab$Y_name <- Y_name
summary_tab$X_name <- X_name
summary_tab$M1_name <- M1_name
summary_tab$M2_name <- M2_name
summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)

X_name <- "Nitrospira"
Y_name <- "HostP_FASN"
# mediator_name <- str_c(c("HostM_cg01442799","HostT_SYK"),collapse = "-")
M1_name <- "HostM_cg02608019"
M2_name <- "HostT_FASN"

PRAD_model <- bruceR::PROCESS(Meddata_scale, y=Y_name, x=X_name,
                              meds=c(M1_name,M2_name),
                              med.type="serial",
                              ci="boot", nsim=999, seed=42)# or omit "mod.type", default is "2-way"
# LGG_model <- bruceR::PROCESS(Meddata_scale, y=M2_name, x=X_name,
#                              meds=c(M1_name),
#                              med.type="serial",
#                              ci="boot", nsim=999, seed=42) 
summary_tab <- PRAD_model$results[[1]]$lavaan.mediation
# print(summary_tab[8,4] < 0.05)
summary_tab$Path_simple <- rownames(summary_tab)
summary_tab$Path <- rownames(summary_tab)
# replace
summary_tab$Path <- gsub("X", X_name,summary_tab$Path)
summary_tab$Path <- gsub("M2", M2_name, summary_tab$Path)
summary_tab$Path <- gsub("Y", Y_name, summary_tab$Path)
summary_tab$Path <- gsub("M1", M1_name, summary_tab$Path)

summary_tab <- summary_tab %>%
  mutate(Significance = case_when(
    (BootULCI == 0 | BootLLCI == 0) & pval > 0.05 ~ "ns",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.001 ~ "***",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.01 ~ "**",
    (BootULCI != 0 & BootLLCI != 0) & pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))
Total_effect <- summary_tab[6,1]
summary_tab$Prop <- summary_tab[,1]/Total_effect

summary_tab$project <- project
summary_tab$biolink <- paste0("Y:[",Y_name,"]~","Treat:[",X_name,"]+ Mediator:[",M1_name,"|",M2_name,"]")
summary_tab$Y_name <- Y_name
summary_tab$X_name <- X_name
summary_tab$M1_name <- M1_name
summary_tab$M2_name <- M2_name
summary_tab_filepath <- paste0("../result/mediation/",project,"/HostP/")
ifelse(!dir.exists(summary_tab_filepath), dir.create(summary_tab_filepath,recursive = TRUE), FALSE)
write.table(summary_tab, file = paste(summary_tab_filepath,project,"_",X_name,"_affects_",Y_name, "_through_",M1_name,"+",M2_name,".txt",sep = ""),
            quote = F, row.names = T, sep = "\t")

mediation_all_tab <- rbind(mediation_all_tab,summary_tab)



####################
write.table(mediation_all_tab, file = paste("../result/mediation/mediation_alls.txt",sep = ""),
            quote = F, row.names = T, sep = "\t")
write.table(mediation_all_tab, file = paste("../result/mediation/mediation_alls.csv",sep = ""),
            quote = F, row.names = F, sep = ",")

# write.table(mediation_all_tab, file = paste("../result/mediation/mediation_cells.txt",sep = ""),
#             quote = F, row.names = T, sep = "\t")
