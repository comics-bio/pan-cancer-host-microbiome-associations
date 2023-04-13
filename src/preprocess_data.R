###################################################################
## @File        :   01_process_data.R
## @Time        :   2022/11/13 09:41:03
## @Desc        :   Process data before downstream analysis for xena data
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
pr.cutoff <- 0.05
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

filter_genes <- function(feat.gene, qt,pr.cutoff){
  if(dim(feat.gene)[1]==0){
    print("control sample is empty")
  }else {
    ## feat.gene rownames is sample
    feat.gene.t <- t(feat.gene)
    feat.gene.t <- feat.gene.t[which(rowSums(as.matrix(feat.gene.t)!=min((feat.gene.t)))>pr.cutoff*dim(feat.gene.t)[2]),]
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
  
}

filter_feature <- function(feat.gene, qt){
  ## feat.gene rownames is sample
  feat.name <- colnames(feat.gene)
  feat.gene.t <- as.data.frame(t(feat.gene))
  genes.sd <- transform(feat.gene.t, SD=apply(feat.gene.t,1, sd, na.rm = TRUE))
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
# data loc has change as xena 
feat.gene.raw <- load_gene_expr(paste0(data.loc,"raw/xena-gene-expression-tpm/tcga_RSEM_gene_tpm"))
cat("Read  gene expression RNAseq from",paste0(data.loc,"raw/xena-gene-expression-tpm/tcga_RSEM_gene_tpm"))

cat("\n dim(feat.gene.raw) is",dim(feat.gene.raw))


meta.all <- load_microbiome_abnd(paste0(data.loc,"raw/Knight_2020/Metadata-TCGA-All-18116-Samples.csv"))

cat("\n dim(meta.all) is",dim(meta.all))


## load microbiome data
feat.otu.relt <- load_microbiome_abnd(paste0(data.loc,"raw/Knight_2020/Kraken-TCGA-Voom-SNM-All-Putative-Contaminants-Removed-Data.csv"))

cat("\n dim(feat.otu.relt)is",dim(feat.otu.relt))

feat.pro.rppa <- read.table(file=paste0(data.loc,"raw/xena-RPPA-pancan-clean/TCGA-RPPA-pancan-clean.xena"),sep="",head=T, row.names = 1, check.names =F, stringsAsFactors = F)
feat.pro.rppa <- as.matrix(t(feat.pro.rppa))
feat.pro.rppa[is.na(feat.pro.rppa)] <- 0
cat("\n dim(feat.pro.rppa)is",dim(feat.pro.rppa))


# feat.dna.meth <-
meth_annotation <-  read.table(file=paste0(data.loc,"raw/xena-dna-methylation/humanmethylation450_15017482_v1-2.csv"),sep=",",head=T, row.names = 1, check.names =F, stringsAsFactors = F)
# feat.meth <- load_gene_expr(paste0(data.loc,"raw/xena-dna-methylation/HumanMethylation27"))
table(meth_annotation$Relation_to_UCSC_CpG_Island)
# Island N_Shelf N_Shore S_Shelf S_Shore
# 89338   75269   12418   31702   11073   24445
meth_island_annotation <- subset(meth_annotation,Relation_to_UCSC_CpG_Island == "Island")
meth_island_annotation$Gene_Name <- meth_island_annotation$UCSC_RefGene_Name
meth_island_annotation$Gene_CpG_Islands_Name <- paste0(meth_island_annotation$UCSC_CpG_Islands_Name,"_",meth_island_annotation$UCSC_RefGene_Name)
## remove probes that match chromosomes X and Y
meth_island_annotation <- subset(meth_island_annotation,!(CHR %in% c("Y")))
## remove SNPs overlapped probe
meth_island_annotation <- subset(meth_island_annotation,Probe_SNPs=="")

# ## Removing probes that have been demonstrated to map to multiple places in the genome.
# # list adapted from https://www.tandfonline.com/doi/full/10.4161/epi.23470

crs.reac <- read.csv(file= paste0(data.loc,"raw/xena-dna-methylation/cross_reactive_probe.chen2013.csv"))
crs.reac <- crs.reac$TargetID[-1]

meth_island_annotation <- subset(meth_island_annotation,!(Name %in% crs.reac))

feat.meth.raw <- load_gene_expr(paste0(data.loc,"raw/xena-dna-methylation/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena"))

# ## remove probes with NA
probe.na <- rowSums(is.na(feat.meth.raw))

table(probe.na == 0)
#FALSE   TRUE
#103553 382024
# chose those has no NA values in rows
probe <- probe.na[probe.na == 0]
feat.meth.raw <- feat.meth.raw[row.names(feat.meth.raw) %in% names(probe), ]
meth_cpg_idx <- intersect(rownames(feat.meth.raw),rownames(meth_island_annotation))
feat.meth.cpg <- feat.meth.raw[meth_cpg_idx,]
# rownames(feat.meth.cpg) <- meth_island_annotation[rownames(feat.meth.cpg),"Gene_CpG_Islands_Name"]

feat.meth.cpg.mval <- t(apply(feat.meth.cpg, 1, function(x) log2((x+0.00001)/((1-x)+0.00001))))

cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
cat("\n++++++++++++++++++++Metadata handle+++++++++++++++++++++++\n")

## change meta rownames as aliquot_ids
meta.all<- meta.all %>%
  as.data.frame()

# meta.all$case_uuid <- tolower(kraken_meta$case_uuid)
# meta.all$sample_uuid <- tolower(kraken_meta$sample_uuid)
meta.all$aliquot_uuid <- tolower(meta.all$aliquot_uuid)
meta.all <- meta.all[order(
  meta.all$investigation, meta.all$experimental_strategy,
  meta.all$case_uuid, meta.all$sample_uuid,
  meta.all$aliquot_uuid, meta.all$data_submitting_center_label
), ]




kraken_msg_pad <- 19
sample_msg_pad <- 5
cat(
  "[", str_pad("Kraken Metadata", kraken_msg_pad, side="right"), "]",
  nrow(meta.all), "samples", length(unique(meta.all$case_uuid)),
  "unique cases\n"
)
cat("\nTranslate Aliquot UUID to TCGA Barcode.\n")
# UUID2barcode <- UUIDtoBarcode(unique(meta.all$aliquot_uuid), from_type = "aliquot_ids")
# if offline
UUID2barcode  <- readRDS("../data/UUID2barcode.RDS")
rownames(UUID2barcode) <- UUID2barcode$portions.analytes.aliquots.aliquot_id
head(UUID2barcode)
#                                case_id submitter_id
# 1 0304b12d-7640-4150-a581-2eea2b1f2ad5 TCGA-OR-A5LL
# 2 075dbfd0-9cf4-4877-884f-ae858902c79e TCGA-OR-A5J7
# 3 0824b246-9fa2-4a8b-ad4c-1ffc7731bf7d TCGA-P6-A5OG
# 4 08e0d412-d4d8-4d13-b792-a4dd0bd9ec2b TCGA-OR-A5J8
# 5 0acc28bc-d348-41df-8f49-bc0332eb9be4 TCGA-OR-A5JL
# 6 0b73cbba-5520-4610-b649-912d76114033 TCGA-OR-A5JP

# meta.all <- meta.all[!duplicated(meta.all$case_uuid, fromLast = F),]
cat("\nFilter duplicated data.")
cat(
  "[", str_pad("Kraken Filter duplicated", kraken_msg_pad, side="right"), "]",
  nrow(meta.all), "samples", length(unique(meta.all$case_uuid)),
  "unique cases\n"
)

meta.all <-  tibble::rownames_to_column(meta.all,"Knight_ID")

meta.all$TCGA_barcode <-UUID2barcode[meta.all$aliquot_uuid,"portions.analytes.aliquots.submitter_id"]
# 
# meta.all[which(meta.all$sample_type == "Primary Tumor"),"TCGA_barcode"] <- paste0(meta.all[which(meta.all$sample_type == "Primary Tumor"),"TCGA_barcode"],"-01")
# meta.all[which(meta.all$sample_type == "Solid Tissue Normal"),"TCGA_barcode"] <- paste0(meta.all[which(meta.all$sample_type == "Solid Tissue Normal"),"TCGA_barcode"],"-11")
# meta.all[which(meta.all$sample_type == "Blood Derived Normal"),"TCGA_barcode"] <- paste0(meta.all[which(meta.all$sample_type == "Blood Derived Normal"),"TCGA_barcode"],"-10")

meta.all$TCGA_barcode_simple <- substr(meta.all$TCGA_barcode,start=1,stop=15) 

# RNA-seq filters
meta.all.rna <- meta.all[meta.all$experimental_strategy == "RNA-Seq",]

cat(
  "[", str_pad("Kraken RNA-seq ", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(meta.all.rna), sample_msg_pad, side="left"), "samples",
  str_pad(length(unique(meta.all.rna$case_uuid)), sample_msg_pad, side="left"),
  "unique cases\n"
)
meta.all.rna <- meta.all.rna[!duplicated(meta.all.rna$TCGA_barcode_simple, fromLast = F),]
rownames(meta.all.rna) <- meta.all.rna$TCGA_barcode_simple
cat(
  "[", str_pad("Kraken RNA-seq After deduplicated", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(meta.all.rna), sample_msg_pad, side="left"), "samples",
  str_pad(length(unique(meta.all.rna$TCGA_barcode_simple)), sample_msg_pad, side="left"),
  "unique cases\n"
)
###
meta.all.WGS <- meta.all[meta.all$experimental_strategy == "WGS", ]

cat(
  "[", str_pad("Kraken WGS", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(meta.all.WGS), sample_msg_pad, side="left"), "samples",
  str_pad(length(unique(meta.all.WGS$case_uuid)), sample_msg_pad, side="left"),"unique cases",
  str_pad(length(unique(meta.all.WGS$TCGA_barcode_simple)), sample_msg_pad, side="left"),
  "unique cases sample\n"
)

meta.all.WGS <- meta.all.WGS[!duplicated(meta.all.WGS$TCGA_barcode_simple, fromLast = F),]
rownames(meta.all.WGS) <- meta.all.WGS$TCGA_barcode_simple
cat(
  "[", str_pad("Kraken WGS After deduplicated", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(meta.all.WGS), sample_msg_pad, side="left"), "samples",
  str_pad(length(unique(meta.all.WGS$case_uuid)), sample_msg_pad, side="left"),
  "unique cases\n"
)

WGS.sample.idx.valid <- intersect(rownames(feat.otu.relt),meta.all.WGS[,"Knight_ID"])
cat("\nNo.Valid WGS sample both in Metadata and microbiome data is",length(WGS.sample.idx.valid))
feat.otu.relt.WGS.valid <- feat.otu.relt[WGS.sample.idx.valid,]
rownames(meta.all.WGS) <- meta.all.WGS$Knight_ID
rownames(feat.otu.relt.WGS.valid) <- meta.all.WGS[rownames(feat.otu.relt.WGS.valid),"TCGA_barcode_simple"]
WGS.sample.idx.valid <- intersect(colnames(feat.gene.raw),rownames(feat.otu.relt.WGS.valid))
feat.gene.raw.WGS <- feat.gene.raw[,WGS.sample.idx.valid]
feat.otu.relt.WGS.valid <- feat.otu.relt.WGS.valid[WGS.sample.idx.valid,]
rownames(meta.all.WGS) <- meta.all.WGS$TCGA_barcode_simple
meta.all.WGS.valid <- meta.all.WGS[WGS.sample.idx.valid,]


## Gene handle
# convert name
gene_convert_tab <- read.table(paste0(data.loc,"raw/xena-gene-expression-tpm/probeMap_gencode.v23.annotation.gene.probemap"),sep="\t",head=T, row.names = 1, check.names =F, stringsAsFactors = F)
rownames(feat.gene.raw.WGS) <- gene_convert_tab[rownames(feat.gene.raw.WGS),"gene"]
gene.name <- rownames(feat.gene.raw.WGS)

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
feat.gene.protein.WGS.valid <- feat.gene.raw.WGS[gene.proteincoding.keep,]
feat.gene.protein.WGS.valid <- t(feat.gene.protein.WGS.valid)
cat("\n Dim(feat.gene.protein.valid) is",dim(feat.gene.protein.WGS.valid),"\n")
## protein
feat.pro.rppa.WGS.valid <- feat.pro.rppa[intersect(rownames(feat.pro.rppa),rownames(meta.all.WGS)),]

cat("\n Dim(feat.pro.rppa.WGS.valid) is",dim(feat.pro.rppa.WGS.valid),"\n")
feat.meth.cpg.mval.WGS.valid <- t(feat.meth.cpg.mval)
feat.meth.cpg.mval.WGS.valid <-feat.meth.cpg.mval.WGS.valid[intersect(rownames(feat.meth.cpg.mval.WGS.valid),rownames(meta.all.WGS)),]
# Ensure same sampleIDs in both genes and microbes data
cat("\nCheck whether data have same sample number.")
stopifnot(all(rownames(feat.gene.protein.WGS.valid) == rownames(feat.otu.relt.WGS.valid)))
stopifnot(all(rownames(feat.gene.protein.WGS.valid) == rownames(meta.all.WGS.valid)))
stopifnot(all(rownames(feat.otu.relt.WGS.valid) == rownames(meta.all.WGS.valid)))


output_dirname <- paste0(data.loc, tag, "/xena_tcga_WGS/")
ifelse(!dir.exists(output_dirname), dir.create(output_dirname,recursive = TRUE), FALSE)


write.table(meth_island_annotation,
            sep="\t",
            file =paste0(output_dirname,"Meth-TCGA-Xena-CpG-All-Mvalue-NonSNP-NonY-Annotation.tsv"))



write.table(feat.otu.relt.WGS.valid,
            sep="\t",
            file =paste0(output_dirname,"Kraken-TCGA-Voom-SNM-All-",dim(feat.otu.relt.WGS.valid)[1],'.tsv'))


write.table(feat.gene.protein.WGS.valid,
            sep="\t",
            file = paste0(output_dirname,"Geneexpr-TCGA-RSEM-TPM-ProteinCoding-All-",dim(feat.otu.relt.WGS.valid)[1],'.tsv'))

write.table(meta.all.WGS,
            sep="\t",
            file = paste0(output_dirname,"Metadata-TCGA-Valid-WGS-All-",dim(meta.all.WGS)[1],'.tsv'))
write.table(feat.pro.rppa.WGS.valid,
            sep="\t",
            file = paste0(output_dirname,"Protein-TCGA-Xena-RPPA-All-",dim(feat.otu.relt.WGS.valid)[1],'.tsv'))

write.table(feat.meth.cpg.mval.WGS.valid,
            sep="\t",
            file = paste0(output_dirname,"Meth-TCGA-Xena-CpG-All-Mvalue-NonSNP-NonY-",dim(feat.meth.cpg.mval.WGS.valid)[1],'.tsv'))

# feat.meth.cpg.mval.WGS.valid<-load_gene_expr("../data/Preprocessed/xena_tcga_WGS/Meth-TCGA-Xena-CpG-All-Mvalue-NonSNP-NonY-1884.tsv")

#get final valid idx 
rna.sample.idx.valid <- intersect(rownames(feat.otu.relt),meta.all.rna[,"Knight_ID"])

cat("\nNo.Valid RNA sample both in Metadata and microbiome data is",length(rna.sample.idx.valid))
# 10125
feat.otu.relt.rna.valid <- feat.otu.relt[rna.sample.idx.valid,]
rownames(meta.all.rna) <- meta.all.rna$Knight_ID
rownames(feat.otu.relt.rna.valid) <- meta.all.rna[rownames(feat.otu.relt.rna.valid),"TCGA_barcode_simple"]

rna.sample.idx.valid <- intersect(colnames(feat.gene.raw),rownames(feat.otu.relt.rna.valid))
cat("\nNo.Valid RNA sample both in Gene and microbiome data is",length(rna.sample.idx.valid))
feat.gene.raw.valid <- feat.gene.raw[,rna.sample.idx.valid]
feat.otu.relt.rna.valid <- feat.otu.relt.rna.valid[rna.sample.idx.valid,]
# rownames(meta.all ) <-meta.all$ID
# rownames(meta.all.valid) <-meta.all.valid$ID
# rownames(meta.all.valid) <-meta.all.valid$TCGA_barcode_simple

rownames(meta.all.rna) <- meta.all.rna$TCGA_barcode_simple
meta.all.rna.valid <- meta.all.rna[rna.sample.idx.valid,]




## Gene handle
# convert name
gene_convert_tab <- read.table(paste0(data.loc,"raw/xena-gene-expression-tpm/probeMap_gencode.v23.annotation.gene.probemap"),sep="\t",head=T, row.names = 1, check.names =F, stringsAsFactors = F)
rownames(feat.gene.raw.valid) <- gene_convert_tab[rownames(feat.gene.raw.valid),"gene"]
gene.name <- rownames(feat.gene.raw.valid)

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
cat("\n Dim(feat.gene.protein.valid) is",dim(feat.gene.protein.valid))
## protein
feat.pro.rppa.rna.valid <- feat.pro.rppa[intersect(rownames(feat.pro.rppa),rownames(meta.all.rna)),]

feat.meth.cpg.mval.rna.valid <- t(feat.meth.cpg.mval)
feat.meth.cpg.mval.rna.valid<- feat.meth.cpg.mval.rna.valid[intersect(rownames(feat.meth.cpg.mval.rna.valid),rownames(meta.all.rna)),]
## Ensure same sampleIDs in both genes and microbes data
cat("\nCheck whether data have same sample number.")
stopifnot(all(rownames(feat.gene.protein.valid) == rownames(feat.otu.relt.rna.valid)))
stopifnot(all(rownames(feat.gene.protein.valid) == rownames(meta.all.rna.valid)))
stopifnot(all(rownames(feat.otu.relt.rna.valid) == rownames(meta.all.rna.valid)))

cat("\nCheck Pass\n")

output_dirname <- paste0(data.loc, tag, "/xena_tcga/")
ifelse(!dir.exists(output_dirname), dir.create(output_dirname,recursive = TRUE), FALSE)


write.table(feat.otu.relt.rna.valid,
            sep="\t",
            file =paste0(output_dirname,"Kraken-TCGA-Voom-SNM-All-",dim(feat.otu.relt.rna.valid)[1],'.tsv'))


write.table(feat.gene.protein.valid,
            sep="\t",
            file = paste0(output_dirname,"Geneexpr-TCGA-RSEM-TPM-ProteinCoding-All-",dim(feat.otu.relt.rna.valid)[1],'.tsv'))

write.table(meta.all.rna,
            sep="\t",
            file = paste0(output_dirname,"Metadata-TCGA-Valid-RNA-All-",dim(meta.all.rna)[1],'.tsv'))
write.table(feat.pro.rppa.rna.valid,
            sep="\t",
            file = paste0(output_dirname,"Protein-TCGA-Xena-RPPA-All-",dim(feat.pro.rppa.rna.valid)[1],'.tsv'))


write.table(feat.meth.cpg.mval.rna.valid,
            sep="\t",
            file = paste0(output_dirname,"Meth-TCGA-Xena-CpG-All-Mvalue-NonSNP-NonY-",dim(feat.meth.cpg.mval.rna.valid)[1],'.tsv'))

write.table(meth_island_annotation,
            sep="\t",
            file =paste0(output_dirname,"Meth-TCGA-Xena-CpG-All-Mvalue-NonSNP-NonY-Annotation.tsv"))

# feat.meth.cpg.mval.RNA.valid <-load_gene_expr("../data/Preprocessed/xena_tcga/Meth-TCGA-Xena-CpG-All-Mvalue-NonSNP-NonY-10245.tsv")

cat("Successfully preprocessed of whole datasets \n")

cat("Start dividing the data for each TCGA project \n")

## Print summary
fn_summary <- paste0("../log/Data_summary.txt")
sink(fn_summary)

cat(
  "[", str_pad("Kraken RNA-seq After deduplicated", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(meta.all.rna.valid), sample_msg_pad, side="left"), "samples",
  str_pad(length(unique(meta.all.rna.valid$TCGA_barcode_simple)), sample_msg_pad, side="left"),
  "unique cases\n"
)
cat(
  "[", str_pad("Gene expression RNA-seq After deduplicated and protein coding without filter", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(feat.gene.protein.valid), sample_msg_pad, side="left"), "samples",
  str_pad(ncol(feat.gene.protein.valid), sample_msg_pad, side="left"), "Gene feature\n"
)
cat(
  "[", str_pad("Protein RPPARNA-seq After deduplicated without filter", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(feat.pro.rppa.rna.valid), sample_msg_pad, side="left"), "samples",
  str_pad(ncol(feat.pro.rppa.rna.valid), sample_msg_pad, side="left"), "Protein feature\n"
)
cat(
  "[", str_pad("DNA methlylation Mvalue RNA-seq After deduplicated without filter", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(feat.meth.cpg.mval.rna.valid), sample_msg_pad, side="left"), "samples",
  str_pad(ncol(feat.meth.cpg.mval.rna.valid), sample_msg_pad, side="left"), "CpG island feature\n"
)

cat("\n\n++++++++++++++++++++++++++++++++++++\n\n")
cat(
  "[", str_pad("Kraken WGS After deduplicated", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(meta.all.WGS), sample_msg_pad, side="left"), "samples",
  str_pad(length(unique(meta.all.WGS$case_uuid)), sample_msg_pad, side="left"),
  "unique cases\n"
)
cat(
  "[", str_pad("Kraken WGS After deduplicated", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(meta.all.WGS.valid), sample_msg_pad, side="left"), "samples",
  str_pad(length(unique(meta.all.WGS.valid$TCGA_barcode_simple)), sample_msg_pad, side="left"),
  "unique cases\n"
)
cat(
  "[", str_pad("Taxa WGS After deduplicated  without filter", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(feat.otu.relt.WGS.valid), sample_msg_pad, side="left"), "samples",
  str_pad(ncol(feat.otu.relt.WGS.valid), sample_msg_pad, side="left"), "Taxa feature\n"
)
cat(
  "[", str_pad("Gene expression WGS After deduplicated and protein coding without filter", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(feat.gene.protein.WGS.valid), sample_msg_pad, side="left"), "samples",
  str_pad(ncol(feat.gene.protein.WGS.valid), sample_msg_pad, side="left"), "Gene feature\n"
)
cat(
  "[", str_pad("Protein RPPA WGS After deduplicated without filter", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(feat.pro.rppa.WGS.valid), sample_msg_pad, side="left"), "samples",
  str_pad(ncol(feat.pro.rppa.WGS.valid), sample_msg_pad, side="left"), "Protein feature\n"
)
cat(
  "[", str_pad("DNA methlylation Mvalue WGS After deduplicated without filter", kraken_msg_pad, side="right"), "]",
  str_pad(nrow(feat.meth.cpg.mval.WGS.valid), sample_msg_pad, side="left"), "samples",
  str_pad(ncol(feat.meth.cpg.mval.WGS.valid), sample_msg_pad, side="left"), "CpG island feature\n"
)
cat("\n\n++++++++++++++++++++++++++++++++++++\n\n")
sink()



##########################################################################
##########################################################################
## Splitting data by project name
meta.all.rna.valid$investigation <- as.factor(meta.all.rna.valid$investigation)
projects <- levels(meta.all.rna.valid$investigation)

for (project in projects) {
  #debug
  # project <- "TCGA-ACC"
  meta.project.rna <- meta.all.rna  %>% 
    subset(investigation == project) %>%
    subset(sample_type %in% c("Primary Tumor","Solid Tissue Normal","Metastatic"))
  
  feat.otu.relt.rna.projcet <- feat.otu.relt.rna.valid[intersect(rownames(feat.otu.relt.rna.valid),rownames(meta.project.rna)),]
  feat.gene.protein.rna.project <- feat.gene.protein.valid[intersect(rownames(feat.gene.protein.valid),rownames(meta.project.rna)),]
  
  stopifnot(all(rownames(feat.gene.protein.rna.project) == rownames(feat.otu.relt.rna.projcet)))
  #feat.gene.protein.rna.project <- filter_genes(feat.gene.protein.rna.project,qt,pr.cutoff)
  output_dirname <- paste0(data.loc, tag, "/xena_tcga")
  output_dirname <- paste0(output_dirname, "/",project,"/")
  
  ifelse(!dir.exists(output_dirname), dir.create(output_dirname,recursive = TRUE), FALSE)
  
  
  feat.pro.rppa.rna.project <-  feat.pro.rppa.rna.valid[intersect(rownames(feat.pro.rppa.rna.valid),rownames(meta.project.rna)),]
  feat.meth.cpg.mval.rna.project <-  feat.meth.cpg.mval.rna.valid[intersect(rownames(feat.meth.cpg.mval.rna.valid),rownames(meta.project.rna)),]
  
  write.table(feat.otu.relt.rna.projcet,
              sep="\t",
              file =paste0(output_dirname,"Kraken-TCGA-Voom-SNM-",project,"-",dim(feat.otu.relt.rna.projcet)[1],'.tsv'))
  
  write.table(feat.gene.protein.rna.project,
              sep="\t",
              file = paste0(output_dirname,"Geneexpr-TCGA-RSEM-TPM-ProteinCoding-QtFilt-",project,"-",dim(feat.gene.protein.rna.project)[1],'.tsv'))
  
  write.table(meta.project.rna,
              sep="\t",
              file = paste0(output_dirname,"Metadata-TCGA-Valid-",project,"-",dim(meta.project.rna)[1],'.tsv'))
  write.table(feat.pro.rppa.rna.project,
              sep="\t",
              file = paste0(output_dirname,"Protein-TCGA-Xena-RPPA-",project,"-",dim(feat.pro.rppa.rna.project)[1],'.tsv'))
  
  write.table(feat.meth.cpg.mval.rna.project,
              sep="\t",
              file = paste0(output_dirname,"Meth-TCGA-Xena-CpG-Mvalue-NonSNP-NonY-",project,"-",dim(feat.meth.cpg.mval.rna.project)[1],'.tsv'))
  
  cat("Successfully preprocessed of ",project," whole datasets\n")
  
  ###
  ### Case
  
  case_output_dirname <- paste0(output_dirname,"/case/")
  ifelse(!dir.exists(case_output_dirname), dir.create(case_output_dirname), FALSE)
  meta.project.rna.case <- meta.project.rna %>% subset(sample_type == "Primary Tumor")
  
  write.table(meta.project.rna.case,
              sep="\t",
              file = paste0(case_output_dirname,"Metadata-TCGA-Valid-",project,"-case-",dim(meta.project.rna.case)[1],".tsv"))
  
  cat('Metadata case dimensions:', dim(meta.project.rna.case), '\n')
  
  
  feat.gene.protein.rna.project.case <- feat.gene.protein.rna.project[intersect(rownames(meta.project.rna.case),rownames(feat.gene.protein.rna.project)),]
  feat.gene.protein.rna.project.case <- filter_genes(feat.gene.protein.rna.project.case,qt,pr.cutoff)
  
  cat('Feature genes expression case dimensions:', dim(feat.gene.protein.rna.project.case), '\n')
  write.table(feat.gene.protein.rna.project.case,
              sep="\t",
              file = paste0(case_output_dirname,"Geneexpr-TCGA-RSEM-TPM-ProteinCoding-QtFilt-",project,"-case-",dim(feat.gene.protein.rna.project.case)[1],".tsv"))
  
  
  feat.otu.relt.rna.projcet.case <- feat.otu.relt.rna.projcet[intersect(rownames(meta.project.rna.case),rownames(feat.otu.relt.rna.projcet)),]
  #otu case
  if (is.null(feat.otu.relt.rna.projcet.case) | dim(feat.otu.relt.rna.projcet.case)[1]==0){
    cat("feat.otu.relt.rna.projcet.case is empty")
  }else{
    feat.otu.relt.rna.projcet.case <- filter_feature(feat.otu.relt.rna.projcet.case,qt)
    cat('Feature microbiome otu case dimensions:', dim(feat.otu.relt.rna.projcet.case), '\n')
    write.table(feat.otu.relt.rna.projcet.case,
                sep="\t",
                file =paste0(case_output_dirname,"Kraken-TCGA-Voom-SNM-",project,"-case-",dim(feat.otu.relt.rna.projcet.case)[1],".tsv"))
  }
  
  #protein case 
  feat.pro.rppa.rna.projcet.case <- feat.pro.rppa.rna.project[intersect(rownames(meta.project.rna.case),rownames(feat.pro.rppa.rna.project)),]
  if (is.null(feat.pro.rppa.rna.projcet.case) | dim(feat.pro.rppa.rna.projcet.case)[1]==0){
    cat("feat.pro.rppa.rna.projcet.case is empty")
  }else{
    feat.pro.rppa.rna.projcet.case[is.na(feat.pro.rppa.rna.projcet.case)] <- 0
    # na equal 0
    feat.pro.rppa.rna.projcet.case <- filter_feature(feat.pro.rppa.rna.projcet.case,qt)
    cat('Feature Protein rppa  case dimensions:', dim(feat.pro.rppa.rna.projcet.case), '\n')
    write.table(feat.pro.rppa.rna.projcet.case,
                sep="\t",
                file = paste0(case_output_dirname,"Protein-TCGA-Xena-RPPA-",project,"-case-",dim(feat.pro.rppa.rna.projcet.case)[1],".tsv"))
    
  }
  
  #meth case
  feat.meth.cpg.mval.rna.project.case <- feat.meth.cpg.mval.rna.project[intersect(rownames(meta.project.rna.case),rownames(feat.meth.cpg.mval.rna.project)),]
  
  if (is.null(feat.meth.cpg.mval.rna.project.case) | dim(feat.meth.cpg.mval.rna.project.case)[1]==0){
    cat("feat.meth.cpg.mval.rna.project.case is empty")
  }else{
    feat.meth.cpg.mval.rna.project.case <- filter_feature(feat.meth.cpg.mval.rna.project.case,qt)
    cat('feat.meth.cpg.mval.rna.project.case case dimensions:', dim(feat.meth.cpg.mval.rna.project.case), '\n')
    write.table(feat.meth.cpg.mval.rna.project.case,
                sep="\t",
                file = paste0(case_output_dirname,"Meth-TCGA-Xena-CpG-Mvalue-NonSNP-NonY-",project,"-case-",dim(feat.meth.cpg.mval.rna.project.case)[1],".tsv"))
    
  }
  
  cat("Successfully segmented ",project," case datasets\n")
  

  ###
  ### Control
  ctr_output_dirname <- paste0(output_dirname,"/control/")
  ctr.idx <-  rownames(meta.project.rna[which(meta.project.rna[,"sample_type"]=="Solid Tissue Normal"),])#get control id
  
  #metadata ctr
  meta.project.rna.ctr <- meta.project.rna %>% subset(sample_type == "Solid Tissue Normal")
  
  if (is.null(meta.project.rna.ctr ) | length(ctr.idx)<5){
    cat("meta.project.rna.ctr is too small")
    # next;
  }else{
    ifelse(!dir.exists(ctr_output_dirname), dir.create(ctr_output_dirname), FALSE)
    cat("Dim of meta.project.rna.ctr is ",dim(meta.project.rna.ctr),"\n")
    write.table(meta.project.rna.ctr,
                sep="\t",
                file = paste0(ctr_output_dirname,"Metadata-TCGA-Valid-",project,"-ctr-",dim(meta.project.rna.ctr)[1],".tsv"))
    
  }
  
  # gene ctr
  ctr.rna.gene.idx <- intersect(rownames(meta.project.rna.ctr),rownames(feat.gene.protein.rna.project))
  feat.gene.protein.rna.project.ctr <- feat.gene.protein.rna.project[ctr.rna.gene.idx,]
  
  
  if(is.null(feat.gene.protein.rna.project.ctr) |length(ctr.rna.gene.idx)<5){
    cat("feat.gene.protein.rna.project.ctr is too small")
  }else if(ncol(feat.gene.protein.rna.project.ctr) ==0){
    cat("feat.gene.protein.rna.project.ctr col is zero")
  }else{
  feat.gene.protein.rna.project.ctr <- filter_genes(feat.gene.protein.rna.project.ctr,qt,pr.cutoff)
  cat('Feature genes expression control dimensions:', dim(feat.gene.protein.rna.project.ctr), '\n')
  
  write.table(feat.gene.protein.rna.project.ctr,
              sep="\t",
              file = paste0(ctr_output_dirname,"Geneexpr-TCGA-RSEM-TPM-ProteinCoding-QtFilt-",project,"-ctr-",dim(feat.gene.protein.rna.project.ctr)[1],".tsv"))
  }
  
  # gene ctr
  ctr.rna.taxa.idx <- intersect(rownames(meta.project.rna.ctr),rownames(feat.otu.relt.rna.projcet))
  feat.otu.relt.rna.projcet.ctr <- feat.otu.relt.rna.projcet[ctr.rna.taxa.idx,]
  
  
  if(is.null(feat.otu.relt.rna.projcet.ctr) |length(ctr.rna.taxa.idx)<5){
    cat("feat.otu.relt.rna.projcet.ctr is too small")
  }else if(ncol(feat.otu.relt.rna.projcet.ctr) ==0){
    cat("feat.otu.relt.rna.projcet.ctr col is zero")
  }else{
    feat.otu.relt.rna.projcet.ctr <- filter_feature(feat.otu.relt.rna.projcet.ctr,qt)
    cat('Feature microbiome otu case dimensions:', dim(feat.otu.relt.rna.projcet.ctr), '\n')
    write.table(feat.otu.relt.rna.projcet.ctr,
                sep="\t",
                file =paste0(ctr_output_dirname,"Kraken-TCGA-Voom-SNM-",project,"-ctr-",dim(feat.otu.relt.rna.projcet.ctr)[1],".tsv"))
  }

  # protein ctr
  ctr.rna.rppa.idx <- intersect(rownames(meta.project.rna.ctr),rownames(feat.pro.rppa.rna.project))
  feat.pro.rppa.rna.projcet.ctr <- feat.pro.rppa.rna.project[ctr.rna.rppa.idx,]
  
  
  if(is.null(feat.pro.rppa.rna.projcet.ctr)  |length(ctr.rna.rppa.idx)<5){
    cat("feat.meth.cpg.mval.rna.project.ctr is too small")
  }else if(ncol(feat.pro.rppa.rna.projcet.ctr) ==0){
    cat("feat.meth.cpg.mval.rna.project.ctr col is zero")
  }else{
    feat.pro.rppa.rna.projcet.ctr[is.na(feat.pro.rppa.rna.projcet.ctr)] <- 0
    feat.pro.rppa.rna.projcet.ctr <- filter_feature(feat.pro.rppa.rna.projcet.ctr,qt)
    cat('Feature Protein rppa  case dimensions:', dim(feat.pro.rppa.rna.projcet.case), '\n')
    write.table(feat.pro.rppa.rna.projcet.ctr,
                sep="\t",
                file = paste0(ctr_output_dirname,"Protein-TCGA-Xena-RPPA-",project,"-ctr-",dim(feat.pro.rppa.rna.projcet.ctr)[1],".tsv"))
  }
  
  
  #meth ctr
  ctr.rna.meth.idx <- intersect(rownames(feat.meth.cpg.mval.rna.project),rownames(meta.project.rna.ctr))
  feat.meth.cpg.mval.rna.project.ctr <- feat.meth.cpg.mval.rna.project[ctr.rna.meth.idx,]
  
  if(is.null(ctr.rna.meth.idx) |length(ctr.rna.meth.idx)<5){
    cat("feat.meth.cpg.mval.rna.project.ctr is too small")
  }else if(ncol(feat.meth.cpg.mval.rna.project.ctr) ==0){
    cat("feat.meth.cpg.mval.rna.project.ctr col is zero")
  }else{
    feat.meth.cpg.mval.rna.project.ctr[is.na(feat.meth.cpg.mval.rna.project.ctr)] <- 0
    feat.meth.cpg.mval.rna.project.ctr <- filter_feature(feat.meth.cpg.mval.rna.project.ctr,qt)
    cat('feat.meth.cpg.mval.rna.project.ctr dimensions:', dim(feat.meth.cpg.mval.rna.project.ctr), '\n')
    write.table(feat.meth.cpg.mval.rna.project.ctr,
                sep="\t",
                file = paste0(ctr_output_dirname,"Meth-TCGA-Xena-CpG-Mvalue-NonSNP-NonY-",project,"-ctr-",dim(feat.meth.cpg.mval.rna.project.ctr)[1],".tsv"))
  }
  
  
  cat("\nSuccessfully segmented ",project," control datasets\n")

  
  ###
  ### metastatic
  metastatic_output_dirname <- paste0(output_dirname,"/metastatic/")
  mets.idx <-  rownames(meta.project.rna[which(meta.project.rna[,"sample_type"]=="Metastatic"),])#get metastatic id
  
  #metadata mets
  meta.project.rna.mets <- meta.project.rna %>% subset(sample_type == "Metastatic")
  
  if (is.null(meta.project.rna.mets) | length(mets.idx)<5){
    cat("meta.project.rna.mets is too small")
    next;
  }else{
    ifelse(!dir.exists(metastatic_output_dirname), dir.create(metastatic_output_dirname), FALSE)
    cat("Dim of meta.project.rna.mets is ",dim(meta.project.rna.mets),"\n")
    write.table(meta.project.rna.mets,
                sep="\t",
                file = paste0(metastatic_output_dirname,"Metadata-TCGA-Valid-",project,"-mets","-",dim(meta.project.rna.mets)[1],".tsv"))
    
  }
  #gene ctr
  feat.gene.protein.rna.project.mets <- feat.gene.protein.rna.project[intersect(rownames(meta.project.rna.mets),rownames(feat.gene.protein.rna.project)),]
  feat.gene.protein.rna.project.mets <- filter_genes(feat.gene.protein.rna.project.mets,qt,pr.cutoff)
  cat('Feature genes expression mets dimensions:', dim(feat.gene.protein.rna.project.mets), '\n')
  
  write.table(feat.gene.protein.rna.project.mets,
              sep="\t",
              file = paste0(metastatic_output_dirname,"Geneexpr-TCGA-RSEM-TPM-ProteinCoding-QtFilt-",project,"-mets-",dim(feat.gene.protein.rna.project.mets)[1],".tsv"))
  
  
  # otu ctr
  feat.otu.relt.rna.projcet.mets <- feat.otu.relt.rna.projcet[intersect(rownames(meta.project.rna.mets),rownames(feat.otu.relt.rna.projcet)),]
  
  if (is.null(feat.otu.relt.rna.projcet.mets) | dim(feat.otu.relt.rna.projcet.mets)[1]==0){
    cat("feat.otu.relt.rna.projcet.mets is empty")
  }else{
    feat.otu.relt.rna.projcet.mets <- filter_feature(feat.otu.relt.rna.projcet.mets,qt)
    cat('Feature microbiome otu case dimensions:', dim(feat.otu.relt.rna.projcet.mets), '\n')
    write.table(feat.otu.relt.rna.projcet.mets,
                sep="\t",
                file =paste0(metastatic_output_dirname,"Kraken-TCGA-Voom-SNM-",project,"-mets-",dim(feat.otu.relt.rna.projcet.mets)[1],".tsv"))
  }
  
  # protein ctr
  mets.rna.rppa.idx <- intersect(rownames(meta.project.rna.mets),rownames(feat.pro.rppa.rna.project))
  feat.pro.rppa.rna.projcet.mets <- feat.pro.rppa.rna.project[mets.rna.rppa.idx,]
  
  
  if(is.null(feat.pro.rppa.rna.projcet.mets)  |length(mets.rna.rppa.idx)<5){
    cat("feat.meth.cpg.mval.rna.project.mets is too small")
  }else if(ncol(feat.pro.rppa.rna.projcet.mets) ==0){
    cat("feat.meth.cpg.mval.rna.project.metscol is zero")
  }else{
    feat.pro.rppa.rna.projcet.mets[is.na(feat.pro.rppa.rna.projcet.mets)] <- 0
    feat.pro.rppa.rna.projcet.mets <- filter_feature(feat.pro.rppa.rna.projcet.mets,qt)
    cat('Feature Protein rppa  mets dimensions:', dim(feat.pro.rppa.rna.projcet.mets), '\n')
    write.table(feat.pro.rppa.rna.projcet.mets,
                sep="\t",
                file = paste0(metastatic_output_dirname,"Protein-TCGA-Xena-RPPA-",project,"-mets-",dim(feat.pro.rppa.rna.projcet.mets)[1],".tsv"))
  }
  
  
  #meth ctr
  mets.rna.meth.idx <- intersect(rownames(feat.meth.cpg.mval.rna.project),rownames(meta.project.rna.mets))
  feat.meth.cpg.mval.rna.project.mets <- feat.meth.cpg.mval.rna.project[mets.rna.meth.idx,]
  
  if(is.null(mets.rna.meth.idx) |length(mets.rna.meth.idx)<5){
    cat("feat.meth.cpg.mval.rna.project.mets is too small")
  }else if(ncol(feat.meth.cpg.mval.rna.project.mets) ==0){
    cat("feat.meth.cpg.mval.rna.project.mets col is zero")
  }else{
    feat.meth.cpg.mval.rna.project.mets[is.na(feat.meth.cpg.mval.rna.project.mets)] <- 0
    feat.meth.cpg.mval.rna.project.mets <- filter_feature(feat.meth.cpg.mval.rna.project.mets,qt)

    cat('feat.meth.cpg.mval.rna.project.mets dimensions:', dim(feat.meth.cpg.mval.rna.project.mets), '\n')
    write.table(feat.meth.cpg.mval.rna.project.mets,
                sep="\t",
                file = paste0(metastatic_output_dirname,"Meth-TCGA-Xena-CpG-Mvalue-NonSNP-NonY-",project,"-mets-",dim(feat.meth.cpg.mval.rna.project.mets)[1],".tsv"))
  }
  
  
  cat("\nSuccessfully segmented ",project," metastatic datasets\n")
}

###########################################
########## WGS
################################################
cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
cat("\nBegin to handle WGS dataset\n")
cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

for (project in projects) {
  #debug
  # project <- "TCGA-ACC"
  project.idx <- rownames(meta.all.WGS [which(meta.all.WGS[,"investigation"]==project),])#get  id
  
  meta.project.WGS <- meta.all.WGS  %>% 
    subset(investigation == project) %>%
    subset(sample_type %in% c("Primary Tumor","Solid Tissue Normal","Metastatic"))

  if (is.null(meta.project.WGS) | length(project.idx)<5){
    cat("\nmeta.project.WGS is too small\n")
    next;
  }
  feat.otu.relt.WGS.projcet <- feat.otu.relt.WGS.valid[intersect(rownames(feat.otu.relt.WGS.valid),rownames(meta.project.WGS)),]
  feat.gene.protein.WGS.project <- feat.gene.protein.WGS.valid[intersect(rownames(feat.gene.protein.WGS.valid),rownames(meta.project.WGS)),]
  
  stopifnot(all(rownames(feat.gene.protein.WGS.project) == rownames(feat.otu.relt.WGS.projcet)))
  # feat.gene.protein.WGS.project <- filter_genes(feat.gene.protein.WGS.project,qt,pr.cutoff)
  output_dirname <- paste0(data.loc, tag, "/xena_tcga_WGS")
  output_dirname <- paste0(output_dirname, "/",project,"/")
  
  ifelse(!dir.exists(output_dirname), dir.create(output_dirname,recursive = TRUE), FALSE)
  
  
  feat.pro.rppa.WGS.project <-  feat.pro.rppa.WGS.valid[intersect(rownames(feat.pro.rppa.WGS.valid),rownames(meta.project.WGS)),]
  feat.meth.cpg.mval.WGS.project <-  feat.meth.cpg.mval.WGS.valid[intersect(rownames(feat.meth.cpg.mval.WGS.valid),rownames(meta.project.WGS)),]
  
  write.table(feat.otu.relt.WGS.projcet,
              sep="\t",
              file =paste0(output_dirname,"Kraken-TCGA-Voom-SNM-",project,"-",dim(feat.otu.relt.WGS.projcet)[1],'.tsv'))
  
  write.table(feat.gene.protein.WGS.project,
              sep="\t",
              file = paste0(output_dirname,"Geneexpr-TCGA-RSEM-TPM-ProteinCoding-QtFilt-",project,"-",dim(feat.gene.protein.WGS.project)[1],'.tsv'))
  
  write.table(meta.project.WGS,
              sep="\t",
              file = paste0(output_dirname,"Metadata-TCGA-Valid-",project,"-",dim(meta.project.WGS)[1],'.tsv'))
  write.table(feat.pro.rppa.WGS.project,
              sep="\t",
              file = paste0(output_dirname,"Protein-TCGA-Xena-RPPA-",project,"-",dim(feat.pro.rppa.WGS.project)[1],'.tsv'))
  
  write.table(feat.meth.cpg.mval.WGS.project,
              sep="\t",
              file = paste0(output_dirname,"Meth-TCGA-Xena-CpG-Mvalue-NonSNP-NonY-",project,"-",dim(feat.meth.cpg.mval.WGS.project)[1],'.tsv'))
  
  cat("Successfully preprocessed of ",project," whole datasets\n")
  
  ###
  ### Case
  
  case_output_dirname <- paste0(output_dirname,"/case/")
  ifelse(!dir.exists(case_output_dirname), dir.create(case_output_dirname), FALSE)
  meta.project.WGS.case <- meta.project.WGS %>% subset(sample_type == "Primary Tumor")
  
  #gene wgs
  feat.gene.protein.WGS.project.case <- feat.gene.protein.WGS.project[intersect(rownames(meta.project.WGS.case),rownames(feat.gene.protein.WGS.project)),]
  feat.gene.protein.WGS.project.case <- filter_genes(feat.gene.protein.WGS.project.case,qt,pr.cutoff)
  
  cat('Feature genes expression case dimensions:', dim(feat.gene.protein.WGS.project.case), '\n')
  write.table(feat.gene.protein.WGS.project.case,
              sep="\t",
              file = paste0(case_output_dirname,"Geneexpr-TCGA-RSEM-TPM-ProteinCoding-QtFilt-",project,"-case-",dim(feat.gene.protein.WGS.project.case)[1],".tsv"))
  #otu case
  feat.otu.relt.WGS.projcet.case <- feat.otu.relt.WGS.projcet[intersect(rownames(meta.project.WGS.case),rownames(feat.otu.relt.WGS.projcet)),]
  feat.otu.relt.WGS.projcet.case <- filter_feature(feat.otu.relt.WGS.projcet.case,qt)
  cat('Feature microbiome otu case dimensions:', dim(feat.otu.relt.WGS.projcet.case), '\n')
  
  #protein case
  case.WGS.rppa.idx <- intersect(rownames(meta.project.WGS.case),rownames(feat.pro.rppa.WGS.project))
  
  feat.pro.rppa.WGS.projcet.case <- feat.pro.rppa.WGS.project[case.WGS.rppa.idx,]
  if(is.null(feat.pro.rppa.WGS.projcet.case) |length(case.WGS.rppa.idx)<3){
    cat("feat.pro.rppa.WGS.projcet.case is too small")
  }else if(ncol(feat.pro.rppa.WGS.projcet.case) <4){
    cat("feat.pro.rppa.WGS.projcet.case feature is too small")
  }else{
  # na equal 0
  feat.pro.rppa.WGS.projcet.case[is.na(feat.pro.rppa.WGS.projcet.case)] <- 0
  feat.pro.rppa.WGS.projcet.case<- filter_feature(feat.pro.rppa.WGS.projcet.case,qt)
  cat('Feature Protein rppa  case dimensions:', dim(feat.pro.rppa.WGS.projcet.case), '\n')
  write.table(feat.pro.rppa.WGS.projcet.case,
              sep="\t",
              file = paste0(case_output_dirname,"Protein-TCGA-Xena-RPPA-",project,"-case-",dim(feat.pro.rppa.WGS.projcet.case)[1],".tsv"))
  }
  
  
  
  #meth case
  case.WGS.meth.idx <- intersect(rownames(meta.project.WGS.case),rownames(feat.meth.cpg.mval.WGS.project))
  feat.meth.cpg.mval.WGS.project.case <- feat.meth.cpg.mval.WGS.project[case.WGS.meth.idx,]
  
  if(is.null(feat.meth.cpg.mval.WGS.project.case) |length(case.WGS.meth.idx)<3){
    cat("feat.meth.cpg.mval.WGS.project.case is too small")
  }else if(ncol(feat.meth.cpg.mval.WGS.project.case) <4){
    cat("feat.meth.cpg.mval.WGS.project.case feature is too small")
  }else{
  feat.meth.cpg.mval.WGS.project.case[is.na(feat.meth.cpg.mval.WGS.project.case)] <- 0
  feat.meth.cpg.mval.WGS.project.case <- filter_feature(feat.meth.cpg.mval.WGS.project.case,qt)
  write.table(feat.meth.cpg.mval.WGS.project.case,
              sep="\t",
              file = paste0(case_output_dirname,"Meth-TCGA-Xena-CpG-Mvalue-NonSNP-NonY-",project,"-case-",dim(feat.meth.cpg.mval.WGS.project.case)[1],".tsv"))
  
  }
  write.table(feat.otu.relt.WGS.projcet.case,
              sep="\t",
              file =paste0(case_output_dirname,"Kraken-TCGA-Voom-SNM-",project,"-case-",dim(feat.otu.relt.WGS.projcet.case)[1],".tsv"))
  
  
  
  write.table(meta.project.WGS.case,
              sep="\t",
              file = paste0(case_output_dirname,"Metadata-TCGA-Valid-",project,"-case-",dim(meta.project.WGS.case)[1],".tsv"))
  
  
  
  cat('Metadata case dimensions:', dim(meta.project.WGS.case), '\n')
  cat("Successfully segmented ",project," case datasets\n")
  
  ###
  ### Control
  ctr_output_dirname <- paste0(output_dirname,"/control/")
  ifelse(!dir.exists(ctr_output_dirname), dir.create(ctr_output_dirname), FALSE)
  # ctr.idx <-  rownames(meta.project.WGS[which(meta.project.WGS[,"sample_type"]=="Solid Tissue Normal"),])#get control id
  # if(length(ctr.idx)<5){
  #   print("Well, control sample is too small")
  #   next;
  # }
  
  ctr.WGS.idx <- rownames(meta.project.WGS[which(meta.project.WGS[,"sample_type"]=="Solid Tissue Normal"),])#get control id
  meta.project.WGS.ctr <- meta.project.WGS %>% subset(sample_type == "Solid Tissue Normal")
  
  if (is.null(meta.project.WGS.ctr) | length(ctr.WGS.idx)<5){
    cat("meta.project.WGS.ctr is too small")
    # next;
  }else{

    cat("\nDim of meta.project.WGS.ctr is",dim(meta.project.WGS.ctr),"\n")
    write.table(meta.project.WGS.ctr,
                sep="\t",
                file = paste0(ctr_output_dirname,"Metadata-TCGA-Valid-",project,"-ctr-",dim(meta.project.WGS.ctr)[1],".tsv"))
  }
  
  #gene ctr
  ctr.WGS.gene.idx <- intersect(rownames(meta.project.WGS.ctr),rownames(feat.gene.protein.WGS.project))
  feat.gene.protein.WGS.project.ctr <- feat.gene.protein.WGS.project[ctr.WGS.gene.idx,]
  
    
    
  if(is.null(feat.gene.protein.WGS.project.ctr) |length(ctr.WGS.gene.idx)<3){
      cat("feat.gene.protein.WGS.project.ctr is too small")
    }else if(ncol(feat.gene.protein.WGS.project.ctr) <4){
      cat("feat.gene.protein.WGS.project.ctr feature is too small")
    }else{
    feat.gene.protein.WGS.project.ctr <- filter_genes(feat.gene.protein.WGS.project.ctr,qt,pr.cutoff)
    cat('Feature genes expression control dimensions:', dim(feat.gene.protein.WGS.project.ctr), '\n')
    
    write.table(feat.gene.protein.WGS.project.ctr,
                sep="\t",
                file = paste0(ctr_output_dirname,"Geneexpr-TCGA-RSEM-TPM-ProteinCoding-QtFilt-",project,"-ctr-",dim(feat.gene.protein.WGS.project.ctr)[1],".tsv"))
  }
  ##otu ctr
  ctr.WGS.otu.idx <- intersect(rownames(meta.project.WGS.ctr),rownames(feat.otu.relt.WGS.projcet))
  feat.otu.relt.WGS.projcet.ctr <- feat.otu.relt.WGS.projcet[ctr.WGS.otu.idx,]
  if(is.null(feat.otu.relt.WGS.projcet.ctr) |length(ctr.WGS.otu.idx)<3){
    cat("feat.otu.relt.WGS.projcet.ctr is too small")
  }else if(ncol(feat.otu.relt.WGS.projcet.ctr) <4){
    cat("feat.otu.relt.WGS.projcet.ctr feature is too small")
  }else{
    feat.otu.relt.WGS.projcet.ctr <- filter_feature(feat.otu.relt.WGS.projcet.ctr,qt)
    cat('Feature microbiome otu case dimensions:', dim(feat.otu.relt.WGS.projcet.ctr), '\n')
    
    write.table(feat.otu.relt.WGS.projcet.ctr,
                sep="\t",
                file =paste0(ctr_output_dirname,"Kraken-TCGA-Voom-SNM-",project,"-ctr-",dim(feat.otu.relt.WGS.projcet.ctr)[1],".tsv"))
  }
  
  
  # protein ctr
  ctr.WGS.rppa.idx <- intersect(rownames(feat.pro.rppa.WGS.project),rownames(meta.project.WGS.ctr))
  feat.pro.rppa.WGS.projcet.ctr <- feat.pro.rppa.WGS.project[ctr.WGS.rppa.idx,]

    
  if(is.null(feat.pro.rppa.WGS.projcet.ctr) |length(ctr.WGS.rppa.idx)<3){
      cat("feat.pro.rppa.WGS.projcet.ctr is too small")
    }else if(ncol(feat.pro.rppa.WGS.projcet.ctr) <4){
      cat("feat.pro.rppa.WGS.projcet.ctr feature is too small")
    }else{
    # na equal 0
    feat.pro.rppa.WGS.projcet.ctr[is.na(feat.pro.rppa.WGS.projcet.ctr)] <- 0
    feat.pro.rppa.WGS.projcet.ctr <- filter_feature(feat.pro.rppa.WGS.projcet.ctr,qt)
    cat('Feature Protein rppa  case dimensions:', dim(feat.pro.rppa.WGS.projcet.ctr), '\n')
    write.table(feat.pro.rppa.WGS.projcet.ctr,
                sep="\t",
                file = paste0(ctr_output_dirname,"Protein-TCGA-Xena-RPPA-",project,"-ctr-",dim(feat.pro.rppa.WGS.projcet.ctr)[1],".tsv"))
  }
  
  # meth ctr
  ctr.WGS.meth.idx <- intersect(rownames(feat.meth.cpg.mval.WGS.project),rownames(meta.project.WGS.ctr))
  feat.meth.cpg.mval.WGS.project.ctr <- feat.meth.cpg.mval.WGS.project[ctr.WGS.meth.idx,]
  
  if(is.null(feat.meth.cpg.mval.WGS.project.ctr) |length(ctr.WGS.meth.idx)<3){
    cat("feat.pro.rppa.WGS.projcet.ctr is too small")
  }else if(ncol(feat.meth.cpg.mval.WGS.project.ctr) <4){
    cat("feat.pro.rppa.WGS.projcet.ctr feature is too small")
  }else{
    feat.meth.cpg.mval.WGS.project.ctr[is.na(feat.meth.cpg.mval.WGS.project.ctr)] <- 0
    feat.meth.cpg.mval.WGS.project.ctr <- filter_feature(feat.meth.cpg.mval.WGS.project.ctr,qt)
    cat('feat.meth.cpg.mval.WGS.project.ctr dimensions:', dim(feat.meth.cpg.mval.WGS.project.ctr), '\n')
    write.table(feat.meth.cpg.mval.WGS.project.ctr,
                sep="\t",
                file = paste0(ctr_output_dirname,"Meth-TCGA-Xena-CpG-Mvalue-NonSNP-NonY-",project,"-ctr-",dim(feat.meth.cpg.mval.WGS.project.ctr)[1],".tsv"))
  }
  
  cat('Metadata control dimensions:', dim(meta.project.WGS.ctr), '\n')

 ##############################
  ### metastatic
  metastatic_output_dirname <- paste0(output_dirname,"/metastatic/")
  mets.idx <-  rownames(meta.project.WGS[which(meta.project.WGS[,"sample_type"]=="Metastatic"),])#get metastatic id
  ifelse(!dir.exists(metastatic_output_dirname), dir.create(metastatic_output_dirname), FALSE)
  #metadata mets
  meta.project.WGS.mets <- meta.project.WGS %>% subset(sample_type == "Metastatic")
  
  if (is.null(meta.project.WGS.mets) | length(mets.idx)<5){
    cat("meta.project.WGS.mets is too small")
    next;
  }else{

    cat("Dim of meta.project.WGS.mets is ",dim(meta.project.WGS.mets),"\n")
    write.table(meta.project.WGS.mets,
                sep="\t",
                file = paste0(metastatic_output_dirname,"Metadata-TCGA-Valid-",project,"-mets","-",dim(meta.project.WGS.mets)[1],".tsv"))
    
  }
  #gene mets
  mets.WGS.gene.idx <- intersect(rownames(meta.project.WGS.mets),rownames(feat.gene.protein.WGS.project))
  feat.gene.protein.WGS.project.mets <- feat.gene.protein.WGS.project[mets.WGS.gene.idx,]
  
    
    
  if(is.null(feat.gene.protein.WGS.project.mets ) |length(mets.WGS.gene.idx)<3){
      cat("feat.gene.protein.WGS.project.mets  is too small")
    }else if(ncol(feat.gene.protein.WGS.project.mets ) <4){
      cat("feat.gene.protein.WGS.project.mets  feature is too small")
    }else{
  feat.gene.protein.WGS.project.mets <- filter_genes(feat.gene.protein.WGS.project.mets,qt,pr.cutoff)
  cat('Feature genes expression mets dimensions:', dim(feat.gene.protein.WGS.project.mets), '\n')
  
  write.table(feat.gene.protein.WGS.project.mets,
              sep="\t",
              file = paste0(metastatic_output_dirname,"Geneexpr-TCGA-RSEM-TPM-ProteinCoding-QtFilt-",project,"-mets-",dim(feat.gene.protein.WGS.project.mets)[1],".tsv"))
  }

  

  
  ##otu mets
  mets.WGS.otu.idx <- intersect(rownames(meta.project.WGS.mets),rownames(feat.otu.relt.WGS.projcet))

  feat.otu.relt.WGS.projcet.mets <- feat.otu.relt.WGS.projcet[mets.WGS.otu.idx,]

  if(is.null(feat.otu.relt.WGS.projcet.mets) |length(mets.WGS.otu.idx)<3){
    cat("feat.otu.relt.WGS.projcet.mets is too small")
  }else if(ncol(feat.otu.relt.WGS.projcet.mets) <4){
    cat("feat.otu.relt.WGS.projcet.mets feature is too small")
  }else{
    feat.otu.relt.WGS.projcet.mets <- filter_feature(feat.otu.relt.WGS.projcet.mets,qt)
    cat('Feature microbiome otu case dimensions:', dim(feat.otu.relt.WGS.projcet.mets), '\n')
    write.table(feat.otu.relt.WGS.projcet.mets,
                sep="\t",
                file =paste0(metastatic_output_dirname,"Kraken-TCGA-Voom-SNM-",project,"-mets-",dim(feat.otu.relt.WGS.projcet.mets)[1],".tsv"))
  }


  # protein ctr
  mets.WGS.rppa.idx <- intersect(rownames(meta.project.WGS.mets),rownames(feat.pro.rppa.WGS.project))
  feat.pro.rppa.WGS.projcet.mets <- feat.pro.rppa.WGS.project[mets.WGS.rppa.idx,]
  
  
  if(is.null(feat.pro.rppa.WGS.projcet.mets)  |length(mets.WGS.rppa.idx)<5){
    cat("feat.meth.cpg.mval.WGS.project.mets is too small")
  }else if(ncol(feat.pro.rppa.WGS.projcet.mets) ==0){
    cat("feat.meth.cpg.mval.WGS.project.metscol is zero")
  }else{
    feat.pro.rppa.WGS.projcet.mets[is.na(feat.pro.rppa.WGS.projcet.mets)] <- 0
    feat.pro.rppa.WGS.projcet.mets <- filter_feature(feat.pro.rppa.WGS.projcet.mets,qt)
    cat('Feature Protein rppa  mets dimensions:', dim(feat.pro.rppa.WGS.projcet.mets), '\n')
    write.table(feat.pro.rppa.WGS.projcet.mets,
                sep="\t",
                file = paste0(metastatic_output_dirname,"Protein-TCGA-Xena-RPPA-",project,"-mets-",dim(feat.pro.rppa.WGS.projcet.mets)[1],".tsv"))
  }
  
  
  #meth ctr
  mets.WGS.meth.idx <- intersect(rownames(feat.meth.cpg.mval.WGS.project),rownames(meta.project.WGS.mets))
  feat.meth.cpg.mval.WGS.project.mets <- feat.meth.cpg.mval.WGS.project[mets.WGS.meth.idx,]
  
  if(is.null(mets.WGS.meth.idx) |length(mets.WGS.meth.idx)<5){
    cat("feat.meth.cpg.mval.WGS.project.mets is too small")
  }else if(ncol(feat.meth.cpg.mval.WGS.project.mets) ==0){
    cat("feat.meth.cpg.mval.WGS.project.mets col is zero")
  }else{
    feat.meth.cpg.mval.WGS.project.mets[is.na(feat.meth.cpg.mval.WGS.project.mets)] <- 0
    feat.meth.cpg.mval.WGS.project.mets <- filter_feature(feat.meth.cpg.mval.WGS.project.mets,qt)

    cat('feat.meth.cpg.mval.WGS.project.mets dimensions:', dim(feat.meth.cpg.mval.WGS.project.mets), '\n')
    write.table(feat.meth.cpg.mval.WGS.project.mets,
                sep="\t",
                file = paste0(metastatic_output_dirname,"Meth-TCGA-Xena-CpG-Mvalue-NonSNP-NonY-",project,"-mets-",dim(feat.meth.cpg.mval.WGS.project.mets)[1],".tsv"))
  }
  
  
  cat("\nSuccessfully segmented ",project," metastatic datasets\n")


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

