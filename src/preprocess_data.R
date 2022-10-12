###################################################################
## @File        :   01_process_data.R
## @Time        :   2022/09/29 09:41:03
## @Desc        :   Process data before downstream analysis
## @Platform    :   All Linux Based Platform
## @Author      :   Su Changxing
## @Version     :   1.0
## @Contact     :   cheonghing.so@gmail.com                      
## @License     :   (C)Copyright 2017-2018, Liugroup-NLPR-CASIA  
###################################################################

.cran_packages <- c("tidyverse",#tidyverse数据处理和可视化
                    "data.table",
                    "yaml")#load yaml config

# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)) {
  #Here we use sustech cran
  options(repos=structure(c(CRAN="https://mirrors.sustech.edu.cn/CRAN/")))
  install.packages(.cran_packages[!.inst])
}

.bioc_packages <- c("TCGAutils",
                    "biomaRt"#To only keep data for protein-coding genes
                    )

# Install Bioconductor packages (if not already installed)
if (!require("BiocManager", quietly = TRUE)){
  options(repos=structure(c(CRAN="https://mirrors.sustech.edu.cn/CRAN/")))
  install.packages("BiocManager")
}


.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  BiocManager::install(.bioc_packages[!.inst])
}


## Import the packages
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
cat("Susscessfully loaded packages")

## This script relies heavily on the R package 'biomaRt',check its version
cat("BiomaRt version is ")
packageVersion('biomaRt')

cat('Starting data preprocee script\n')
start.time <- proc.time()[1]


## load parameters
parameters <- yaml.load_file('../parameters.yaml')
cat("Loaded yaml config\n")
set.seed(parameters$public$seed)
data.loc <- parameters$public$data.loc
tag <- parameters$preprocess$tag
qt <- parameters$preprocess$qt
cat("Preprocessed Data will be save in ",data.loc,"\n")
# ----------------------------------------------------------------------
#                              Function 
# ----------------------------------------------------------------------

load_gene_expr <- function(filename){
  genes <- data.frame(data.table::fread(filename,sep="\t",head=T), row.names = 1, check.names = F)
  genes <- as.matrix(genes)
  
}

load_microbiome_abnd <- function(filename){
  microbes <- read.table(filename,sep=",",head=T, row.names = 1, check.names =F, stringsAsFactors = F)
  microbes <- as.matrix(microbes)
  
}

filter_genes <- function(feat.gene, qt){
  ## feat.gene rownames is sample
  feat.gene.t <- t(feat.gene)
  genes.sd <- transform(as.data.frame(feat.gene.t), SD=apply(as.data.frame(feat.gene.t),1, sd, na.rm = TRUE))
  ## select top genes with high SD (~ variability) across samples
  SD_quantile <- quantile(genes.sd$SD) ## identical to summary(genes.sd$SD)
  SD_cutoff <- SD_quantile[qt] ## 2nd quantile -- 25th quantile.
  genes.sd <- genes.sd[order(genes.sd$SD, decreasing = T),]
  top.variable.genes <- rownames(genes.sd[genes.sd$SD > SD_cutoff,])
  ## subset these genes from gene table
  select <- which(colnames(feat.gene) %in% top.variable.genes)
  feat.gene <- feat.gene[,select]
}

# ----------------------------------------------------------------------
#                              Main 
# ----------------------------------------------------------------------

cat("Read raw data\n")

feat.gene.raw <- load_gene_expr(paste0(data.loc,"raw/PanCanAtlas/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"))
dim(feat.gene.raw)


meta.all <- load_microbiome_abnd(paste0(data.loc,"raw/Knight_2020/Metadata-TCGA-All-18116-Samples.csv"))
dim(meta.all)

## load microbiome data
feat.otu.relt <- load_microbiome_abnd(paste0(data.loc,"raw/Knight_2020/Kraken-TCGA-Voom-SNM-All-Putative-Contaminants-Removed-Data.csv"))
dim(feat.otu.relt)

## change meta rownames as aliquot_ids
meta.all<- meta.all %>%
  as.data.frame()
meta.all <- meta.all[!duplicated(meta.all$aliquot_uuid, fromLast = F),]
meta.all <-  tibble::rownames_to_column(meta.all,"ID")
rownames(meta.all) <- meta.all[,"aliquot_uuid"]

# Transform TCGAbarcode to UUID
TCGAbarcode <- colnames(feat.gene.raw)
# UUID.idx <- TCGAutils::barcodeToUUID(TCGAbarcode)
##If offline
UUID.idx <- readRDS("../data/TCGAbarcode2UUID.RDS")
UUID.idx$aliquot_ids<- toupper(UUID.idx$aliquot_ids)#Toupper
rownames(UUID.idx) <- UUID.idx$aliquot_ids
##get the submitter_aliquot_ids
meta.all[,"submitter_aliquot_ids"] <- UUID.idx[rownames(meta.all),"submitter_aliquot_ids"]

sum(!is.na(meta.all$submitter_aliquot_ids))
#11059

#get valid idx between gene expr and metadata
meta.all.valid <- meta.all[!is.na(meta.all$submitter_aliquot_ids),]

#get final valid idx 
sample.idx.valid <- intersect(rownames(feat.otu.relt),meta.all.valid[,"ID"])
length(sample.idx.valid)
# 10769
feat.otu.relt.valid <- feat.otu.relt[sample.idx.valid,]
rownames(meta.all ) <-meta.all$ID
rownames(meta.all.valid) <-meta.all.valid$ID
meta.all.valid <- meta.all[sample.idx.valid,]
rownames(feat.otu.relt.valid) <- meta.all.valid[rownames(feat.otu.relt.valid),"submitter_aliquot_ids"]
feat.gene.raw.valid <- feat.gene.raw[,rownames(feat.otu.relt.valid)]
rownames(meta.all.valid) <-meta.all.valid$submitter_aliquot_ids


gene.name <-str_split_fixed(rownames(feat.gene.raw.valid),"\\|",2)[,1]#get the gene name
rownames(feat.gene.raw.valid) <- gene.name
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",host="https://www.ensembl.org")
##If offline


#used the ‘biomaRt’ R package to only keep data for protein-coding genes
#ref:https://www.biostars.org/p/9517934/#9517940
# transcript.biotype <- biomaRt::getBM(
#   attributes=c("hgnc_symbol","transcript_biotype"),
#   #filters = c("transcript_biotype"), 
#   values=gene.name, mart=ensembl)
transcript.biotype <- readRDS("../data/transcript_biotype.RDS")

gene.proteincoding.keep <- transcript.biotype[which(transcript.biotype$"transcript_biotype"=="protein_coding"),]

gene.proteincoding.keep <- intersect(gene.proteincoding.keep$"hgnc_symbol",gene.name)
feat.gene.protein.valid <- feat.gene.raw.valid[gene.proteincoding.keep,]
feat.gene.protein.valid <- t(feat.gene.protein.valid)
dim(feat.gene.protein.valid)
# 
# ## 
# rownames(meta.all) <- meta.all$ID
# sample.idx <- intersect(rownames(feat.gene.protein),rownames(meta.all))
# sample.idx <- intersect(sample.idx,rownames(feat.otu.relt))
# cat('Feature genes expression dimensions:', length(sample.idx), '\n')
# 
# ## get data from idx
# feat.gene.protein.valid <- feat.gene.protein[sample.idx,]
# feat.otu.relt.valid <- feat.otu.relt[sample.idx,]
# meta.all.valid <- meta.all[sample.idx,]

## Ensure same sampleIDs in both genes and microbes data
stopifnot(all(rownames(feat.gene.protein.valid) == rownames(feat.otu.relt.valid)))
stopifnot(all(rownames(feat.gene.protein.valid) == rownames(meta.all.valid)))
stopifnot(all(rownames(feat.otu.relt.valid) == rownames(meta.all.valid)))


output_dirname <- paste0(data.loc, tag, "/")
ifelse(!dir.exists(output_dirname), dir.create(output_dirname), FALSE)

feat_relt_all.loc <- paste0(data.loc, tag, "/","Kraken-TCGA-Voom-SNM-",dim(feat.otu.relt.valid)[1],'.tsv')


write.table(feat.otu.relt.valid,
            sep="\t",
            file =paste0(output_dirname,"Kraken-TCGA-Voom-SNM-All-",dim(feat.otu.relt.valid)[1],'.tsv'))


write.table(feat.gene.protein.valid,
            sep="\t",
            file = paste0(output_dirname,"Geneexpr-TCGA-ProteinCoding-All-",dim(feat.otu.relt.valid)[1],'.tsv'))

write.table(meta.all.valid,
            sep="\t",
            file = paste0(output_dirname,"Metadata-TCGA-Valid-All-",dim(meta.all.valid)[1],'.tsv'))

cat("Successfully preprocessed of whole datasets \n")

cat("Start dividing the data for each TCGA project \n")
## Splitting data by project name
meta.all.valid$investigation <- as.factor(meta.all.valid$investigation)
projects <- levels(meta.all.valid$investigation)

for (project in projects) {
  meta.project <- meta.all.valid[which(meta.all.valid$investigation == project),]
  meta.project <- meta.all.valid %>% 
                    subset(investigation == project) %>%
                    subset(sample_type %in% c("Primary Tumor","Solid Tissue Normal"))
  
  feat.otu.relt.projcet <- feat.otu.relt.valid[rownames(meta.project),]
  feat.gene.protein.project <- feat.gene.protein.valid[rownames(meta.project),]
  
  stopifnot(all(rownames(feat.gene.protein.project) == rownames(feat.otu.relt.projcet)))
  feat.gene.protein.project <- filter_genes(feat.gene.protein.project,qt)
  
  output_dirname <- paste0(data.loc, tag, "/",project,"/")
  ifelse(!dir.exists(output_dirname), dir.create(output_dirname), FALSE)
  
  write.table(feat.otu.relt.projcet,
              sep="\t",
              file =paste0(output_dirname,"Kraken-TCGA-Voom-SNM-",project,"-",dim(feat.otu.relt.projcet)[1],'.tsv'))
  
  
  write.table(feat.gene.protein.project,
              sep="\t",
              file = paste0(output_dirname,"Geneexpr-TCGA-ProteinCoding-QtFilt-",project,"-",dim(feat.gene.protein.project)[1],'.tsv'))
  
  write.table(meta.project,
              sep="\t",
              file = paste0(output_dirname,"Metadata-TCGA-Valid-",project,"-",dim(meta.project)[1],'.tsv'))
  cat("Successfully preprocessed of ",project," whole datasets\n")
  
  ###
  ### Case
  
  case_output_dirname <- paste0(data.loc, tag, "/",project,"/case/")
  ifelse(!dir.exists(case_output_dirname), dir.create(case_output_dirname), FALSE)
  
  case.idx <-  rownames(meta.project[which(meta.project[,"sample_type"]=="Primary Tumor"),])#get case id
  feat.gene.protein.project.case <- feat.gene.protein.project[case.idx,]
  cat('Feature genes expression case dimensions:', dim(feat.gene.protein.project.case), '\n')
  write.table(feat.gene.protein.project.case,
              sep="\t",
              file = paste0(case_output_dirname,"Geneexpr-TCGA-ProteinCoding-QtFilt-",project,"-case.tsv"))
  
  feat.otu.relt.projcet.case <- feat.otu.relt.projcet[case.idx,]
  cat('Feature microbiome otu case dimensions:', dim(feat.otu.relt.projcet.case), '\n')
  
  write.table(feat.otu.relt.projcet.case,
              sep="\t",
              file =paste0(case_output_dirname,"Kraken-TCGA-Voom-SNM-",project,"-case.tsv"))
  
  meta.project.case <- meta.project[case.idx,]
  write.table(meta.project.case,
              sep="\t",
              file = paste0(case_output_dirname,"Metadata-TCGA-Valid-",project,"-case.tsv"))
  
  cat('Metadata case dimensions:', dim(meta.project.case), '\n')
  cat("Successfully segmented ",project," case datasets\n")

  ###
  ### Control
  ctr_output_dirname <- paste0(data.loc, tag, "/",project,"/control/")
  ifelse(!dir.exists(ctr_output_dirname), dir.create(ctr_output_dirname), FALSE)
  
  ctr.idx <-  rownames(meta.project[which(meta.project[,"sample_type"]=="Solid Tissue Normal"),])#get control id
  feat.gene.protein.project.ctr <- feat.gene.protein.project[ctr.idx,]
  cat('Feature genes expression control dimensions:', dim(feat.gene.protein.project.ctr), '\n')
  write.table(feat.gene.protein.project.ctr,
              sep="\t",
              file = paste0(ctr_output_dirname,"Geneexpr-TCGA-ProteinCoding-QtFilt-",project,"-ctr.tsv"))
  
  feat.otu.relt.projcet.ctr <- feat.otu.relt.projcet[ctr.idx,]
  cat('Feature microbiome otu case dimensions:', dim(feat.otu.relt.projcet.ctr), '\n')
  
  write.table(feat.otu.relt.projcet.ctr,
              sep="\t",
              file =paste0(ctr_output_dirname,"Kraken-TCGA-Voom-SNM-",project,"-ctr.tsv"))
  
  meta.project.ctr <- meta.project[ctr.idx,]
  write.table(meta.project.ctr,
              sep="\t",
              file = paste0(ctr_output_dirname,"Metadata-TCGA-Valid-",project,"-ctr.tsv"))
  
  cat('Metadata case dimensions:', dim(meta.project.ctr), '\n')
  cat("Successfully segmented ",project," control datasets\n")
}

# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",host="https://www.ensembl.org")
# test <- biomaRt::getBM(attributes = c("hgnc_symbol", "entrezgene_id"), 
#                        #filters = "hgnc_symbol", 
#                        values = gene.entrez.id,
#                        bmHeader = T, 
#                        mart = ensembl)



cat('Successfully complete preprocess script in',
    proc.time()[1]-start.time, 'second...\n')

