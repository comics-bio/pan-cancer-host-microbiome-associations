###################################################################
## @File        :   grid_search_sparseCCA.R
## @Time        :   2022/09/29 11:03:48
## @Desc        :   Parallel grid search hyperparameters for sparse CCA analysis
## @Platform    :   All Linux Based Platform
## @Author      :   Su Changxing
## @Version     :   1.0
## @Contact     :   suchangxing@genomics.cn                      
## @License     :   (C)Copyright 2017-2018, Liugroup-NLPR-CASIA  
###################################################################
.cran_packages <- c("tidyverse",#tidyverse数据处理和可视化
                    "PMA",
                    "data.table",
                    "foreach",
                    "doParallel",
                    "itertools",
                    "iterators",
                    "optparse",
                    "yaml")#load yaml config

# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

.bioc_packages <- c()

# Install Bioconductor packages (if not already installed)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  BiocManager::install(.bioc_packages[!.inst])
}

## Import the packages
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
cat("Susscessfully loaded packages")

cat('Starting gird search sparseCCA analysis script\n')
start.time <- proc.time()[1]

## load parameters
parameters <- yaml.load_file('../parameters.yaml')
cat("Loaded yaml config\n")
set.seed(parameters$public$seed)
data.loc <- parameters$public$data.loc
tag <- parameters$SparseCCA$tag
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-p", "--project"), type="character", default= "TCGA-BRCA",
                help="TCGA cancer project [default: %default]"),
    make_option(c("-o","--output"), type="character",default= parameters$public$result.loc,
                help="The analysis result location [default: %default]"),
    make_option(c("-n","--num"), type="integer", default= parameters$SparseCCA$num_cores,
                help="Cluster cores numbers [default: %default]"),
    make_option(c("-f","--fdr"), type="numeric", default=parameters$SparseCCA$FDR.cutoff,
                help="significance cutoff[default: %default]"),
    make_option(c("-c","--componet"), type="integer", default= parameters$SparseCCA$num_componets,
                help="CCA componet[default: %default]")
    
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  result.loc <- paste0(opts$output,tag,"/",opts$project,"/")
  # 显示输入输出确认是否正确
  print(paste("The name of the analyzed data project is ", opts$project,  sep = ""))
  print(paste("Cluster cores numbers is ", opts$num,  sep = ""))
  print(paste("CCA componet is ", opts$componet,  sep = ""))
  print(paste("The significance cutoff is ", opts$fdr, sep = ""))
  print(paste("The analysis result location  is ", result.loc, sep = ""))
}
# 
# 
# # ----------------------------------------------------------------------
# #                              Function 
# # ----------------------------------------------------------------------

tune_params_grid_search_parallel <- function(X, Y,penaltyX,penaltyY){
  num_samples <- nrow(X)
  tune_start_time <- Sys.time()
  #   cl <- makeCluster(100)
  #   registerDoParallel(cl)
  cca_res <- foreach(i =icount(length(penaltyX)),.combine="rbind",.errorhandling = "pass") %:%
    foreach(j =icount(length(penaltyY)),.combine="rbind",.errorhandling = "pass") %:%
    foreach(k=icount(num_samples),.combine="rbind",.packages="PMA",.errorhandling = "pass") %dopar% {
      res <- CCA(X[-k,],Y[-k,], penaltyx = penaltyX[i], penaltyz = penaltyY[j], K=1, niter = 5, trace = F, standardize = T)
      # mat_u <- matrix(unlist(res[k,"u"]), ncol=1)#将列表转化成矩阵
      # mat_v <- matrix(unlist(res[k,"v"]), ncol=1)#将列表转化成矩阵
      score_demo_single <- matrix(ncol = 5 , nrow = 1)
      score_demo_single[,1] <- X[k,]%*%res$u
      score_demo_single[,2] <- Y[k,]%*%res$v
      score_demo_single[,3] <- i
      score_demo_single[,4] <- j
      score_demo_single[,5] <- k
      score_demo_single
    }
  tune_end_time <- Sys.time()
  time_elapsed <- tune_end_time - tune_start_time
  print(paste0("Time elapsed for param tuning = ", time_elapsed))
  return (cca_res)
}

test_significance_LOOCV_parallel <- function(X, Y, bestpenaltyX, bestpenaltyY, num_components){
  cca.k = num_components
  scoresXcv <- matrix(nrow = nrow(X), ncol = cca.k)
  scoresYcv <-  matrix(nrow = nrow(Y), ncol = cca.k)
  corr_pval <- c()
  corr_r <- c()
  score <- foreach(i = icount(nrow(X)),.packages = "PMA",.combine = "rbind") %dopar%{
    res <- CCA(X[-i,],Y[-i,], penaltyx=bestpenaltyX, penaltyz=bestpenaltyY, K=cca.k, trace = F) ## default niter = 15 which is spit out when trace = T (default)
    scoresXcv_single<- matrix(nrow = 1, ncol = cca.k)
    scoresYcv_single <-  matrix(nrow =1, ncol =cca.k)
    scoresXcv_single[,1] <- X[i,]%*%res$u[,1]
    scoresYcv_single[,1] <- Y[i,]%*%res$v[,1]
    scoresXcv_single[,2] <- X[i,]%*%res$u[,2]
    scoresYcv_single[,2] <- Y[i,]%*%res$v[,2]
    scoresXcv_single[,3] <- X[i,]%*%res$u[,3]
    scoresYcv_single[,3] <- Y[i,]%*%res$v[,3]
    scoresXcv_single[,4] <- X[i,]%*%res$u[,4]
    scoresYcv_single[,4] <- Y[i,]%*%res$v[,4]
    scoresXcv_single[,5] <- X[i,]%*%res$u[,5]
    scoresYcv_single[,5] <- Y[i,]%*%res$v[,5]
    scoresXcv_single[,6] <- X[i,]%*%res$u[,6]
    scoresYcv_single[,6] <- Y[i,]%*%res$v[,6]
    scoresXcv_single[,7] <- X[i,]%*%res$u[,7]
    scoresYcv_single[,7] <- Y[i,]%*%res$v[,7]
    scoresXcv_single[,8] <- X[i,]%*%res$u[,8]
    scoresYcv_single[,8] <- Y[i,]%*%res$v[,8]
    scoresXcv_single[,9] <- X[i,]%*%res$u[,9]
    scoresYcv_single[,9] <- Y[i,]%*%res$v[,9]
    scoresXcv_single[,10] <- X[i,]%*%res$u[,10]
    scoresYcv_single[,10] <- Y[i,]%*%res$v[,10]
    score_single <- cbind(scoresXcv_single,scoresYcv_single)
    score_single
  }
  
  ## Test for each components
  for(j in 1:cca.k){
    corr <- cor.test(score[,j],score[,j+10]) ## Pearson correlation.
    corr_pval[j] <- corr$p.value
  }
  corr_pval
}

run_sparseCCA <- function(X, Z, CCA.K, penaltyX, penaltyZ, vInit=NULL, outputFile=NULL){
  CCA.out <-  CCA(X,Z,typex="standard",typez="standard",K=CCA.K,
                  penaltyx=penaltyX,penaltyz=penaltyZ,
                  v=vInit)
  #Output CCA result
  if(!is.null(outputFile)){
    sink(outputFile)
    print(CCA.out)
    sink()
  }
  
  ifelse(!dir.exists(outputFile), dir.create(outputFile), FALSE)
  
  ## add rownames to output factors
  rownames(CCA.out$u) <- colnames(X)
  rownames(CCA.out$v) <- colnames(Z)
  ## compute contribution of selected features to each of the samples.
  CCA_var_genes <- X %*% CCA.out$u ## canonical variance for genes
  CCA_var_microbes <- Z %*% CCA.out$v ## canonical variance for microbes
  
  return(list(CCA.out, CCA_var_genes, CCA_var_microbes))
  
}

save_CCA_components <- function(CCA.out, CCA.K, dirname){
  ## Print canonical covariates in files
  for(i in CCA.K){
    print(paste0("Writing significant component = ", i))
    selected_X <- which(CCA.out$u[,i]!=0)
    selected_X <- rownames(CCA.out$u)[selected_X]
    coeff_X <- unname(CCA.out$u[selected_X,i])
    selected_Z <- which(CCA.out$v[,i]!=0)
    selected_Z <- rownames(CCA.out$v)[selected_Z]
    coeff_Z <- unname(CCA.out$v[selected_Z,i])
    ## Make all vectors of same length to avoid repetition of elements from shorter vectors.
    n <- max(length(selected_X), length(selected_Z))
    length(selected_X) <- n
    length(selected_Z) <- n
    length(coeff_X) <- n
    length(coeff_Z) <- n
    selected_XZ <- as.data.frame(cbind(gene = selected_X, gene_coeff = coeff_X,
                                       taxa = selected_Z, taxa_coeff = coeff_Z))
    write.table(selected_XZ, file=paste0(dirname,"gene_taxa_component_",i,".txt"), sep = "\t", col.names = NA)
  }
  
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


get_avg_features <- function(cca_cov, CCA.K){
  num_features <- 0
  for(k in 1:CCA.K){
    num_features <- num_features + length(which(cca_cov[,k]!=0))
  }
  avg_features <- num_features/CCA.K
}


# # ----------------------------------------------------------------------
# #                              Main
# # ----------------------------------------------------------------------
# ## now initiate parallel processors.
cl <- makeCluster(opts$num) # 8 workers per core -- itasca
registerDoParallel(cl)
print(paste0("Number of workers: ", getDoParWorkers()))

penaltyX <- seq(0.1,0.4,length=10)
penaltyY <- seq(0.15,0.4,length=10)
cca.k <- opts$componet
# ##############################################################################
# Case compute

## load gene expression and microbiome tables
input_dirname <- paste0(data.loc,parameters$preprocess$tag,"/by_type/",opts$project,"/case/")
filenames <- list.files(input_dirname)
cat("Read data from ",input_dirname," .....\n")

feat.gene.expr <- read.table(paste0(input_dirname,grep("Geneexpr",filenames,value = TRUE)))
feat.gene.expr <- as.matrix(feat.gene.expr)
feat.gene.expr[is.na(feat.gene.expr)] <- 0
feat.gene.expr <- filter_genes(feat.gene.expr,2)
feat.gene.expr <- as.matrix(feat.gene.expr)
cat('Feature genes expression dimensions:', dim(feat.gene.expr), '\n')

feat.otu.relt <- read.table(paste0(input_dirname,grep("Kraken",filenames,value = TRUE)))
feat.otu.relt <- as.matrix(feat.otu.relt)
cat('Feature microbiome otu dimensions:', dim(feat.otu.relt), '\n')

## Ensure same sampleIDs in both genes and microbes data before sparse CCA
stopifnot(all(rownames(feat.gene.expr) == rownames(feat.otu.relt)))

meta.project <- read.table(paste0(input_dirname,grep("Metadata",filenames,value = TRUE)))
cat('Metadata dimensions:', dim(meta.project), '\n')


## select tuning parameters using grid-search
tune_cca_mat <- tune_params_grid_search_parallel(feat.gene.expr,
                                                 feat.otu.relt,
                                                 penaltyX,
                                                 penaltyY)
tune_cca <- as.data.frame(tune_cca_mat)

cat('Successfully complete hyperparameters grid search analysis in',
    proc.time()[1]-start.time, 'second...\n')


## identify best penalty parameters
### find index with max absolute corr
colnames(tune_cca) <-  c("scoreXcv","scoreYcv","i","j","k")
corr_demo <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))
for( i in 1:length(penaltyX)){
  for(j in 1:length(penaltyY)){
    data <- tune_cca[which(tune_cca$i==i&tune_cca$j==j),]
    corr_demo[i,j] = cor(data[,"scoreXcv"],data[,"scoreYcv"])
  }
}

row.names(corr_demo) <- as.character(penaltyX)
colnames(corr_demo) <- as.character(penaltyY)
corr_demo_df <- as.data.frame(corr_demo)
rownames(corr_demo_df)
colnames(corr_demo_df)

bestpenalty <- which(abs(corr_demo) == max(abs(corr_demo)), arr.ind = TRUE)

bestpenalty

bestpenaltyX <- penaltyX[bestpenalty[1]]
bestpenaltyY <- penaltyY[bestpenalty[2]]
cat("\nTCGA PROJECT:",opts$project,"case bestpenaltyX is ",bestpenaltyX)
cat("\nTCGA PROJECT:",opts$project,"case bestpenaltyY is ",bestpenaltyY)

## Run sparse CCA using selected tuning param using permutation search

cat("\nRun sparse CCA using selected tuning param using permutation search\n")

## This will return FALSE if the directory already exists or is uncreatable,
## and TRUE if it didn't exist but was succesfully created.
output_dirname <- paste0(result.loc,"geneExp_taxa/case/grid_search_",bestpenaltyX,"_",bestpenaltyY,"/")
ifelse(!dir.exists(output_dirname), dir.create(output_dirname,recursive = TRUE), FALSE)

cca <- run_sparseCCA(feat.gene.expr, feat.otu.relt,
                     cca.k,
                     bestpenaltyX,
                     bestpenaltyY,
                     outputFile=paste0(output_dirname,"CCA_case_output_gridsearch_",bestpenaltyX,"_",bestpenaltyY,".txt"))

## average number of genes and microbes in resulting components
avg_genes <- get_avg_features(cca[[1]]$u, cca.k)
cat("case avg_genes",avg_genes)

avg.microbes <- get_avg_features(cca[[1]]$v, cca.k)
cat("case avg_microbiomes",avg.microbes)

## Test for each components
CCA_pval <- test_significance_LOOCV_parallel(feat.gene.expr,feat.otu.relt, bestpenaltyX, bestpenaltyY, cca.k)

length(which(CCA_pval < opts$fdr))
which(CCA_pval < opts$fdr)

CCA_padj <- p.adjust(CCA_pval, method = "BH")
CCA_padj

length(which(CCA_padj < opts$fdr))
which(CCA_padj < opts$fdr)

## Save significant components
sig_cutoff <- opts$fdr
sig <- which(CCA_padj < sig_cutoff)
sig_components_dirname <- paste0(output_dirname,"sig_components_case_",bestpenaltyX,"_",bestpenaltyY,"/")
ifelse(!dir.exists(sig_components_dirname), dir.create(sig_components_dirname), FALSE)
save_CCA_components(cca[[1]],sig,sig_components_dirname)

## Print summary
fn_summary <- paste0(sig_components_dirname,opts$project,"_CCA_case_summary_",bestpenaltyX,"_",bestpenaltyY,".txt")
sink(fn_summary)
cat(paste0(" bestpenaltyX = ", bestpenaltyX, ", bestpenaltyY = ", bestpenaltyY))
cat(paste0("\n cor(Xu,Yv): \n"))
cat(paste0(signif(cca[[1]]$cors, digits = 4)))
cat(paste0("\n Avg. no. of genes across components = ",avg_genes))
cat(paste0("\n Avg. no. of microbes across components= ", avg.microbes))
cat(paste0("\n P-value for components (LOOCV): \n"))
cat(paste0(signif(CCA_pval, digits = 4)))
# cat(paste0("\n LOOCV corr: \n"))
# cat(paste0(signif(corr_r, digits = 4)))
cat(paste0("\n No. of components with p-value < 0.1 = ", length(which(CCA_pval < 0.1))))
cat(paste0("\n No. of components with p-value < 0.05 = ", length(which(CCA_pval < 0.05))))
cat(paste0("\n No. of components with p-value < 0.005 = ", length(which(CCA_pval < 0.005))))
cat(paste0("\n No. of components with FDR < 0.1 = ", length(which(CCA_padj < 0.1))))
cat(paste0("\n No. of components with FDR < 0.05 = ", length(which(CCA_padj < 0.05))))
cat(paste0("\n No. of components with FDR < 0.005 = ", length(which(CCA_padj < 0.005))))
cat(paste0("\n Significant components with cutoff ",opts$fdr," : \n" ))
cat(paste0(which(CCA_padj < opts$fdr)))
cat(paste0("\n Significant components result location in  ",sig_components_dirname," : \n" ))
sink()


## Print total summary
fn_total <- paste0(opts$output,tag,"/","CCA_case_total_summary.txt")
sink(fn_total,append = TRUE)
cat(paste0("#########################################\n"))
cat(paste0("Analysis Project:",opts$project,"\n"))
cat(paste0(" bestpenaltyX = ", bestpenaltyX, ", bestpenaltyY = ", bestpenaltyY))
cat(paste0("\n cor(Xu,Yv): \n"))
cat(paste0(signif(cca[[1]]$cors, digits = 4)))
cat(paste0("\n Avg. no. of genes across components = ",avg_genes))
cat(paste0("\n Avg. no. of microbes across components= ", avg.microbes))
cat(paste0("\n P-value for components (LOOCV): \n"))
cat(paste0(signif(CCA_pval, digits = 4)))
# cat(paste0("\n LOOCV corr: \n"))
# cat(paste0(signif(corr_r, digits = 4)))
cat(paste0("\n No. of components with p-value < 0.1 = ", length(which(CCA_pval < 0.1))))
cat(paste0("\n No. of components with p-value < 0.05 = ", length(which(CCA_pval < 0.05))))
cat(paste0("\n No. of components with p-value < 0.005 = ", length(which(CCA_pval < 0.005))))
cat(paste0("\n No. of components with FDR < 0.1 = ", length(which(CCA_padj < 0.1))))
cat(paste0("\n No. of components with FDR < 0.05 = ", length(which(CCA_padj < 0.05))))
cat(paste0("\n No. of components with FDR < 0.005 = ", length(which(CCA_padj < 0.005))))
cat(paste0("\n Significant components: \n" ))
cat(paste0(which(CCA_padj < 0.1)))
cat(paste0("\n#########################################\n"))
sink()

cat('Successfully complete',opts$project,' CCA case pval analysis and save results in',
    proc.time()[1]-start.time, 'second...\n')
cat("Result saved in ",sig_components_dirname,".\n")




# ##############################################################################
# Control compute

## load gene expression and microbiome tables
input_dirname <- paste0(data.loc,parameters$preprocess$tag,"/by_type/",opts$project,"/control/")
filenames <- list.files(input_dirname)
cat("Read data from ",input_dirname," .....")

feat.gene.expr <- read.table(paste0(input_dirname,grep("Geneexpr",filenames,value = TRUE)))
feat.gene.expr <- as.matrix(feat.gene.expr)
feat.gene.expr[is.na(feat.gene.expr)] <- 0
feat.gene.expr <- filter_genes(feat.gene.expr,2)
feat.gene.expr <- as.matrix(feat.gene.expr)
if(dim(feat.gene.expr)[1]<10)
  stop("Control samples size is too small.\n")
cat('Feature genes expression dimensions:', dim(feat.gene.expr), '\n')

feat.otu.relt <- read.table(paste0(input_dirname,grep("Kraken",filenames,value = TRUE)))
feat.otu.relt <- as.matrix(feat.otu.relt)
cat('Feature microbiome otu dimensions:', dim(feat.otu.relt), '\n')

## Ensure same sampleIDs in both genes and microbes data before sparse CCA
stopifnot(all(rownames(feat.gene.expr) == rownames(feat.otu.relt)))

meta.project <- read.table(paste0(input_dirname,grep("Metadata",filenames,value = TRUE)))
cat('Metadata dimensions:', dim(meta.project), '\n')
dim(meta.project)

## select tuning parameters using grid-search
tune_cca_mat <- tune_params_grid_search_parallel(feat.gene.expr,
                                                 feat.otu.relt,
                                                 penaltyX,
                                                 penaltyY)
tune_cca <- as.data.frame(tune_cca_mat)

cat('Successfully complete hyperparameters grid search analysis in',
    proc.time()[1]-start.time, 'second...\n')


## identify best penalty parameters
### find index with max absolute corr
colnames(tune_cca) <-  c("scoreXcv","scoreYcv","i","j","k")
corr_demo <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))
for( i in 1:length(penaltyX)){
  for(j in 1:length(penaltyY)){
    data <- tune_cca[which(tune_cca$i==i&tune_cca$j==j),]
    corr_demo[i,j] = cor(data[,"scoreXcv"],data[,"scoreYcv"])
  }
}

row.names(corr_demo) <- as.character(penaltyX)
colnames(corr_demo) <- as.character(penaltyY)
corr_demo_df <- as.data.frame(corr_demo)
rownames(corr_demo_df)
colnames(corr_demo_df)

bestpenalty <- which(abs(corr_demo) == max(abs(corr_demo)), arr.ind = TRUE)

bestpenalty

bestpenaltyX <- penaltyX[bestpenalty[1]]
bestpenaltyY <- penaltyY[bestpenalty[2]]
cat("\nTCGA PROJECT:",opts$project,"control bestpenaltyX is ",bestpenaltyX)
cat("\nTCGA PROJECT:",opts$project,"control bestpenaltyY is ",bestpenaltyY)

## Run sparse CCA using selected tuning param using permutation search

cat("Run sparse CCA using selected tuning param using permutation search\n")

## This will return FALSE if the directory already exists or is uncreatable,
## and TRUE if it didn't exist but was succesfully created.
output_dirname <- paste0(result.loc,"geneExp_taxa/control/grid_search_",bestpenaltyX,"_",bestpenaltyY,"/")
ifelse(!dir.exists(output_dirname), dir.create(output_dirname,recursive = TRUE), FALSE)

cca <- run_sparseCCA(feat.gene.expr, feat.otu.relt,
                     cca.k,
                     bestpenaltyX,
                     bestpenaltyY,
                     outputFile=paste0(output_dirname,"CCA_control_output_gridsearch_",bestpenaltyX,"_",bestpenaltyY,".txt"))

## average number of genes and microbes in resulting components
avg_genes <- get_avg_features(cca[[1]]$u, cca.k)
cat("\ncontrol avg_genes",avg_genes)

avg.microbes <- get_avg_features(cca[[1]]$v, cca.k)
cat("\ncontrol avg_microbiomes",avg.microbes)

## Test for each components
CCA_pval <- test_significance_LOOCV_parallel(feat.gene.expr,feat.otu.relt, bestpenaltyX, bestpenaltyY, cca.k)

length(which(CCA_pval < opts$fdr))
which(CCA_pval < opts$fdr)

CCA_padj <- p.adjust(CCA_pval, method = "BH")
CCA_padj

length(which(CCA_padj < opts$fdropts$fdr))
which(CCA_padj < opts$fdr)

## Save significant components
sig_cutoff <- opts$fdr
sig <- which(CCA_padj < sig_cutoff)
sig_components_dirname <- paste0(output_dirname,"sig_components_control",bestpenaltyX,"_",bestpenaltyY,"/")
ifelse(!dir.exists(sig_components_dirname), dir.create(sig_components_dirname), FALSE)
save_CCA_components(cca[[1]],sig,sig_components_dirname)

## Print summary
fn_summary <- paste0(sig_components_dirname,opts$project,"_CCA_control_summary_",bestpenaltyX,"_",bestpenaltyY,".txt")
sink(fn_summary)
cat(paste0(" bestpenaltyX = ", bestpenaltyX, ", bestpenaltyY = ", bestpenaltyY))
cat(paste0("\n cor(Xu,Yv): \n"))
cat(paste0(signif(cca[[1]]$cors, digits = 4)))
cat(paste0("\n Avg. no. of genes across components = ",avg_genes))
cat(paste0("\n Avg. no. of microbes across components= ", avg.microbes))
cat(paste0("\n P-value for components (LOOCV): \n"))
cat(paste0(signif(CCA_pval, digits = 4)))
# cat(paste0("\n LOOCV corr: \n"))
# cat(paste0(signif(corr_r, digits = 4)))
cat(paste0("\n No. of components with p-value < 0.1 = ", length(which(CCA_pval < 0.1))))
cat(paste0("\n No. of components with p-value < 0.05 = ", length(which(CCA_pval < 0.05))))
cat(paste0("\n No. of components with p-value < 0.005 = ", length(which(CCA_pval < 0.005))))
cat(paste0("\n No. of components with FDR < 0.1 = ", length(which(CCA_padj < 0.1))))
cat(paste0("\n No. of components with FDR < 0.05 = ", length(which(CCA_padj < 0.05))))
cat(paste0("\n No. of components with FDR < 0.005 = ", length(which(CCA_padj < 0.005))))
cat(paste0("\n Significant components with cutoff ",opts$fdr," : \n" ))
cat(paste0(which(CCA_padj < opts$fdr)))
cat(paste0("\n Significant components result location in  ",sig_components_dirname," : \n" ))
sink()

## Print total summary
fn_total <- paste0(opts$output,tag,"/","CCA_control_total_summary.txt")
sink(fn_total,append = TRUE)
cat(paste0("#########################################\n"))
cat(paste0("Analysis Project:",opts$project,"\n"))
cat(paste0(" bestpenaltyX = ", bestpenaltyX, ", bestpenaltyY = ", bestpenaltyY))
cat(paste0("\n cor(Xu,Yv): \n"))
cat(paste0(signif(cca[[1]]$cors, digits = 4)))
cat(paste0("\n Avg. no. of genes across components = ",avg_genes))
cat(paste0("\n Avg. no. of microbes across components= ", avg.microbes))
cat(paste0("\n P-value for components (LOOCV): \n"))
cat(paste0(signif(CCA_pval, digits = 4)))
# cat(paste0("\n LOOCV corr: \n"))
# cat(paste0(signif(corr_r, digits = 4)))
cat(paste0("\n No. of components with p-value < 0.1 = ", length(which(CCA_pval < 0.1))))
cat(paste0("\n No. of components with p-value < 0.05 = ", length(which(CCA_pval < 0.05))))
cat(paste0("\n No. of components with p-value < 0.005 = ", length(which(CCA_pval < 0.005))))
cat(paste0("\n No. of components with FDR < 0.1 = ", length(which(CCA_padj < 0.1))))
cat(paste0("\n No. of components with FDR < 0.05 = ", length(which(CCA_padj < 0.05))))
cat(paste0("\n No. of components with FDR < 0.005 = ", length(which(CCA_padj < 0.005))))
cat(paste0("\n Significant components: \n" ))
cat(paste0(which(CCA_padj < 0.1)))
cat(paste0("\n#########################################\n"))
sink()

cat('Successfully complete',opts$project,' CCA control pval analysis and save results in',
    proc.time()[1]-start.time, 'second...\n')
cat("Result saved in ",sig_components_dirname,".\n")

cat('Successfully complete',opts$project,'whole sparseCCA analysis script in',
    proc.time()[1]-start.time, 'second...\n')