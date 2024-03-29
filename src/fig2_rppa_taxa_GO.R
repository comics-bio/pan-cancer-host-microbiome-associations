.cran_packages <- c("tidyverse",#
                    "data.table",
                    "stringr",
                    "ggraph",
                    "igraph",
                    "tidygraph",
                    "ggforce",
                    "ggpubr",
                    "ggplot2",
                    "patchwork",
                    "data.table",
                    "showtext",
                    "stringr",
                    "msigdbr",
                    "yaml")#load yaml config


# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)) {
  #Here we use sustech cran
  options(repos=structure(c(CRAN="https://mirrors.sustech.edu.cn/CRAN/")))
  install.packages(.cran_packages[!.inst])
}

.bioc_packages <- c("clusterProfiler",
                    "topGO",
                    "pathview",
                    "org.Hs.eg.db",
                    "Rgraphviz"
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
# showtext_auto()
setwd("C:/Lab/project/pan-cancer-host-microbiome-associations/src")
parameters <- yaml.load_file('../parameters.yaml')
cat("Susscessfully loaded yaml config\n")
set.seed(parameters$public$seed)
data.loc <- parameters$public$data.loc
fig.loc <- parameters$public$fig.loc

if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}

if (TRUE){
  option_list = list(
    make_option(c("-p", "--project"), type="character", default= "TCGA-BRCA",
                help="TCGA cancer project [default: %default]"),
    make_option(c("-o","--output"), type="character",default= parameters$public$fig.loc,
                help="The analysis figure location [default: %default]"),
    make_option(c("-n", "--num"), type="numeric", default= 0.15,
                help="Percent [default: %default]")
    
  )
  opts = parse_args(OptionParser(option_list=option_list))
  project <- opts$project
  percent <- opts$num
  fig.loc <- paste0(opts$output,parameters$SparseCCA$tag,"/componet_visualiazation/")

  print(paste("The name of the analyzed data project is ", opts$project,  sep = ""))
  print(paste("The analysis percent is ", opts$num,  sep = ""))
  print(paste("The analysis figure location  is ", fig.loc, sep = ""))
}
msigdb_enrichment <- function(pathway_DB, pathways, background_genes, genes_of_interest){
  
  enrichment_list <- list()
  for(i in 1:length(pathways)){
    
    pathway <- pathways[i]
    
    ## genes in this pathway
    pathway_gene_set <- pathway_DB[pathway_DB$gs_name == pathway,]$human_gene_symbol
    length(pathway_gene_set) 
    ## If the criteria for min and max #genes in a given pathway is not satified, 
    ## skip testing the current pathway
    if(length(pathway_gene_set) < min_genes || length(pathway_gene_set) > max_genes) next
    
    ## The contingency table
    ## Inspired by: https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
    ##                    genes_of_interest genes_NOT_of_interest  Total
    ## In_pathway               x                 m - x              m
    ## Not_in_pathway         k - x             n - (k - x)          n
    ## Total                    k                 m + n - k          m + n
    
    ## m, #overlapping genes in this pathway and background genes for CRC
    m <- length(intersect(background_genes,pathway_gene_set))
    m ## 11
    
    ## n, #genes in background but not in pathway 
    n <- length(setdiff(background_genes,pathway_gene_set))
    n #12502
    
    ## x, #genes of interest in pathway 
    x <- length(intersect(pathway_gene_set,genes_of_interest))
    x ## 1
    ## If the overlap between genes of interest and the pathway is less than the overlap cut-off, 
    ## then skip testing the current pathway
    if(x < overlap_genes) next 
    
    ## Extract list of genes in the genes of interest that are included in the pathway. 
    gene_names_in_pathway = noquote(paste(intersect(pathway_gene_set,genes_of_interest), collapse = ","))
    
    ## k, total #genes in genes list of interest
    k <- length(genes_of_interest)
    k
    
    ## Build the contigency table
    contingency_table <- matrix(c(x, k-x, m-x, n-(k-x)), #matrix is filled along the column by default
                                nrow = 2, ncol = 2, 
                                dimnames = list(c("In_pathway","Not_in_pathway"),
                                                c("Genes_of_interest","Genes_NOT_of_interest"))
    )
    
    contingency_table
    
    fisher_result <- fisher.test( contingency_table, alternative = "greater")
    
    ## save details in a dataframe
    enrichment_result_df <- data.frame( pathway = pathway,
                                        BG_genes = length(background_genes),
                                        genes_in_pathway = length(pathway_gene_set),
                                        genes_in_path_and_BG = m,
                                        genes_of_interest = k,
                                        genes_of_interest_in_pathway = x,
                                        gene_names_in_pathway = gene_names_in_pathway,
                                        ## fill in contingency table entries: x, k-x, m-x, n-(k-x) for z-score computation.
                                        cont_n1 = x,
                                        cont_n2 = k-x,
                                        cont_n3 = m-x,
                                        cont_n4 = n-(k-x),
                                        CI_95 = paste0(signif(fisher_result$conf.int[1],5),"-",signif(fisher_result$conf.int[2],5)),
                                        odds_ratio = unname(fisher_result$estimate),
                                        p_val = fisher_result$p.value
    )
    
    enrichment_list[[i]] <- enrichment_result_df
    
  }
  
  return(enrichment_list)
  
}

## Do enrichment per component
perform_enrichment <- function(input_dir, filenames, msigdb_collection_name, pathway_DB, pathways, background_genes, output_dir){
  enriched_list <- list()
  count <- 0
  
  ## perform enrichment for each component 
  for(i in filenames){
    
    ## debug
    # i <- "gene_taxa_component_1.txt"
    
    count <- count + 1
    
    print(paste0("Enrichment for component: ", i));flush.console()
    ## load genes list for a given component
    sparseCCA_genes <- read.table(paste0(input_dir,"/",i),sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
    sparseCCA_genes <- sparseCCA_genes$gene
    genes_of_interest <- sparseCCA_genes
    
    ## perform enrichment using pathways in msigdb collection
    msigdb_pathways <- msigdb_enrichment(pathway_DB, pathways, background_genes, genes_of_interest) 
    ## convert the list to dataframe
    msigdb_pathways_df <- do.call(rbind.data.frame, msigdb_pathways)
    ## Add current component as a column to enrichment result
    msigdb_pathways_df$Component <- rep(paste0(i),nrow(msigdb_pathways_df))
    ## Add collection name as a column to the enrichemnt result
    msigdb_pathways_df$Collection <- msigdb_collection_name
    
    ## Make last column(i.e., component) as first column
    msigdb_pathways_df <- msigdb_pathways_df[,c(ncol(msigdb_pathways_df)-1, ncol(msigdb_pathways_df),1:(ncol(msigdb_pathways_df)-2))]
    
    ## update column name to reflect component, except pathway colname
    # colnames(msigdb_pathways_df[,-3]) <- paste(colnames(msigdb_pathways_df[,-3]), count, sep = "_")
    
    ## sort pathways by pval
    msigdb_pathways_df <- msigdb_pathways_df[order(msigdb_pathways_df$p_val, decreasing = F),]
    length(msigdb_pathways_df$pathway)
    
    ## save dataframe of pathways for a component to file
    if(!is.null(output_dir)){
      filename <- strsplit(i,"\\.")[[1]][1]
      write.table(msigdb_pathways_df, file = paste0(output_dir,"/msigdb_",filename,".txt") , sep="\t", row.names = F)
    }
    ## add df to list
    enriched_list[[count]] <- msigdb_pathways_df
  }
  
  ## This is the list of dataframes, where each dataframe holds enrichment result for a component.  
  return(enriched_list)
}
plot_taxa_gene <- function(i,tmp_data,top.taxa.idx){
  p1 <- ggplot(tmp_data, mapping = aes(x= tmp_data[,i], y=tmp_data[,5])) +
    geom_point(shape = 21, color = "black", size = 4,fill="#9970AB")+
    geom_smooth(method = lm,color="#9970AB")+
    xlab(paste0("Taxa:",str_split_fixed(top.taxa.idx[i],'g__',2)[,2],"  (r=",signif(cor(tmp_data[,i],tmp_data[,5]),5),")"))+
    ylab(paste0("Gene:",colnames(tmp_data)[5]))+
    ggtitle(paste0(colnames(tmp_data)[5]," ~ ",str_split_fixed(top.taxa.idx[i],'g__',2)[,2]))+
    # annotate("text", x = 4, y = 25, label = paste0("rho:",cor(tmp_data[,1],tmp_data[,2])),
    #          color="#9970AB",size = 5, fontface="bold" )+
    theme_classic()
  
  
  p2 <- ggplot(tmp_data, mapping = aes(x= tmp_data[,i], y=tmp_data[,6])) +
    geom_point(shape = 21, color = "black", size = 4,fill="#9970AB")+
    geom_smooth(method = lm,color="#9970AB")+
    xlab(paste0("Taxa:",str_split_fixed(top.taxa.idx[i],'g__',2)[,2],"  (r=",signif(cor(tmp_data[,i],tmp_data[,6]),5),")"))+
    ylab(paste0("Gene:",colnames(tmp_data)[6]))+
    ggtitle(paste0(colnames(tmp_data)[6]," ~ ",str_split_fixed(top.taxa.idx[i],'g__',2)[,2]))+
    # annotate("text", x = 4, y = 25, label = paste0("rho:",cor(tmp_data[,1],tmp_data[,2])),
    #          color="#9970AB",size = 5, fontface="bold" )+
    theme_classic()
  
  
  p3 <- ggplot(tmp_data, mapping = aes(x= tmp_data[,1], y=tmp_data[,7])) +
    geom_point(shape = 21, color = "black", size = 4,fill="#9970AB")+
    geom_smooth(method = lm,color="#9970AB")+
    xlab(paste0("Taxa:",str_split_fixed(top.taxa.idx[i],'g__',2)[,2],"  (r=",signif(cor(tmp_data[,i],tmp_data[,7]),5),")"))+
    ylab(paste0("Gene:",colnames(tmp_data)[7]))+
    ggtitle(paste0(colnames(tmp_data)[7]," ~ ",str_split_fixed(top.taxa.idx[i],'g__',2)[,2]))+
    # annotate("text", x = 4, y = 25, label = paste0("rho:",cor(tmp_data[,1],tmp_data[,2])),
    #          color="#9970AB",size = 5, fontface="bold" )+
    theme_classic()
  
  p4 <- ggplot(tmp_data, mapping = aes(x= tmp_data[,i], y=tmp_data[,8])) +
    geom_point(shape = 21, color = "black", size = 4,fill="#9970AB")+
    geom_smooth(method = lm,color="#9970AB")+
    xlab(paste0("Taxa:",str_split_fixed(top.taxa.idx[i],'g__',2)[,2],"  (r=",signif(cor(tmp_data[,i],tmp_data[,8]),5),")"))+
    ylab(paste0("Gene:",colnames(tmp_data)[8]))+
    ggtitle(paste0(colnames(tmp_data)[8]," ~ ",str_split_fixed(top.taxa.idx[i],'g__',2)[,2]))+
    # annotate("text", x = 4, y = 25, label = paste0("rho:",cor(tmp_data[,1],tmp_data[,2])),
    #          color="#9970AB",size = 5, fontface="bold" )+
    theme_classic()
  
  plot <- p1/p2/p3/p4
  return(plot)
}

plot_componet <- function(sparseCCA_component,filename){

  cor.coef= sparseCCA_component
  
  cor.coef.gene=cor.coef[,c("gene","gene_coeff"),drop=F]
  rownames(cor.coef.gene)=cor.coef.gene$gene
  cor.coef.gene$gene=NULL
  cor.coef.gene=cor.coef.gene[order(cor.coef.gene$gene_coeff,decreasing = T),,drop=F]
  cor.coef.gene$name=rownames(cor.coef.gene)
  cor.coef.gene=cor.coef.gene[order(cor.coef.gene$gene_coeff),]
  cor.coef.gene$name=factor(cor.coef.gene$name,levels = cor.coef.gene$name)
  cor.coef.gene.tail <- tail(cor.coef.gene,n=10)
  cor.coef.gene.head <- head(cor.coef.gene,n=10)
  cor.coef.gene.top <- rbind(cor.coef.gene.head,cor.coef.gene.tail)
  # plot.gene <- ggplot(cor.coef.gene.top, aes(x=name, y=gene_coeff)) +
  #   geom_bar(stat="identity", color="#CCBB44", fill="white",width=0.05)+xlab("")+ylab("Correlation Coefficient")+coord_flip()+theme_classic()+
  #   theme(axis.text.y=element_blank(),
  #         axis.ticks.y=element_blank())
  
  plot.gene.top <- ggplot(cor.coef.gene.top, aes(x=name, y=gene_coeff)) +
    geom_bar(stat="identity", color="#CCBB44", fill="#CCBB44",width=0.8)+xlab("")+ylab("Gene Correlation Coefficient(Top 10)")+coord_flip()+theme_classic()+
    theme()
  # ggsave("OutputPlot/sparseCCA.plot8.1.pdf",width = 3,height = 5)
  
  cor.coef.taxa=cor.coef[,c("taxa","taxa_coeff"),drop=F]
  cor.coef.taxa=na.omit(cor.coef.taxa)
  cor.coef.taxa <- cor.coef.taxa[str_detect(cor.coef.taxa$taxa,'g__'),]
  rownames(cor.coef.taxa) <- str_split_fixed(cor.coef.taxa$taxa,'g__',2)[,2]
  #cor.coef.taxa$taxa=NULL
  cor.coef.taxa=cor.coef.taxa[order(cor.coef.taxa$taxa_coeff,decreasing = T),,drop=F]
  # cor.coef.taxa=as.data.frame(t(cor.coef.taxa))
  # pdf("./sparseCCA.plot5.2.pdf",width = 8,height = 3)
  # pheatmap(cor.coef.taxa,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = T,
  #          cellheight = 20,cellwidth = 20,color = colorRampPalette(c("#2166AC", "#D1E5F0"))(100),legend = F)
  # dev.off()
  
  # cor.coef.taxa=as.data.frame(t(cor.coef.taxa))
  cor.coef.taxa$name=rownames(cor.coef.taxa)
  # cor.coef.taxa=cor.coef.taxa[cor.coef.taxa$name!="Actinobacteria.1",]
  cor.coef.taxa=cor.coef.taxa[order(cor.coef.taxa$taxa_coeff),]
  cor.coef.taxa$name=factor(cor.coef.taxa$name,levels = cor.coef.taxa$name)
  cor.coef.taxa.tail <- tail(cor.coef.taxa,n=10)
  cor.coef.taxa.head <- head(cor.coef.taxa,n=10)
  cor.coef.taxa.top=rbind(cor.coef.taxa.head,cor.coef.taxa.tail)
  plot.taxa.top <- ggplot(cor.coef.taxa.top, aes(x=taxa_coeff, y=name)) +
    geom_bar(stat="identity", color="grey", fill="#0077BB",width=0.8)+xlab("Taxa Correlation Coefficient(Top 10)")+ylab("")+theme_classic()
  
  top.taxa.idx.top <- head(cor.coef.taxa.top$taxa,2)
  top.taxa.idx.tail <- tail(cor.coef.taxa.top$taxa,2)
  top.taxa.idx <- c(top.taxa.idx.tail,top.taxa.idx.top)
  
  top.gene.idx.head <- head(rownames(cor.coef.gene.top),2)
  top.gene.idx.tail <- tail(rownames(cor.coef.gene.top),2)
  top.gene.idx <- c(top.gene.idx.head,top.gene.idx.tail)
  # print(top.gene.idx)
  
  top.taxa <- as.data.frame(feat.otu.relt[,top.taxa.idx])
  top.gene <- as.data.frame(feat.gene.expr[,top.gene.idx])
  
  tmp_data=cbind(top.taxa,top.gene)
  # tmp_data <- as.data.frame(tmp_data)
  
  # print(colnames(tmp_data))
  
  a <- plot_taxa_gene(1,tmp_data = tmp_data,top.taxa.idx=top.taxa.idx)
  b <- plot_taxa_gene(2,tmp_data = tmp_data,top.taxa.idx=top.taxa.idx)
  c <- plot_taxa_gene(3,tmp_data = tmp_data,top.taxa.idx=top.taxa.idx)
  d <- plot_taxa_gene(4,tmp_data = tmp_data,top.taxa.idx=top.taxa.idx)
  patchwork <-(a|b|c|d)|(plot.taxa.top+plot.gene.top)
  patchwork <- patchwork + patchwork::plot_annotation(title = paste0('SparseCCA component ',filename,' for ',project))
  return(patchwork)
}






opts$project <- "TCGA-HNSC"
project <- opts$project
percent <- opts$num
fig.loc <- paste0(opts$output,parameters$SparseCCA$tag,"/componet_visualiazation/fig2/")
# 显示输入输出确认是否正确
print(paste("The name of the analyzed data project is ", opts$project,  sep = ""))
print(paste("The analysis percent is ", opts$num,  sep = ""))
print(paste("The analysis figure location  is ", fig.loc, sep = ""))

result.loc <- paste0(parameters$public$result.loc,"/SparseCCA_xena/",project,"/")
result.loc 
# Get Data

# result.loc  <- paste0("../figures/SparseCCA/componet_visualiazation/fig2/SparseCCA_xena/",project,"/geneExp_taxa")
input_dirname <- paste0(data.loc,parameters$preprocess$tag,"/xena_tcga/",project,"/case/")
filenames <- list.files(input_dirname)
cat("Read data from ",input_dirname," .....")



feat.rppa <- as.matrix(fread(paste0(input_dirname,grep("RPPA",filenames,value = TRUE))),rownames=1)
cat('Feature microbiome otu dimensions:', dim(feat.rppa), '\n')
# feat.meth.colname <- read.table(paste0(input_dirname,grep("CpG",filenames,value = TRUE)),sep="\t",nrow = 1)
# feat.meth <- fread(paste0(input_dirname,grep("CpG",filenames,value = TRUE)),skip=1,sep="\t")
# feat.meth <- as.matrix(feat.gene.expr)


feat.otu.relt <- read.table(paste0(input_dirname,grep("Kraken",filenames,value = TRUE)))
feat.otu.relt <- as.matrix(feat.otu.relt)
cat('Feature microbiome otu dimensions:', dim(feat.otu.relt), '\n')


idx <- intersect(rownames(feat.otu.relt),rownames(feat.rppa))
feat.otu.relt <- feat.otu.relt[idx,]
feat.rppa <- feat.rppa[idx,]

## Ensure same sampleIDs in both genes and microbes data before sparse CCA
stopifnot(all(rownames(feat.rppa) == rownames(feat.otu.relt)))

meta.project <- read.table(paste0(input_dirname,grep("Metadata",filenames,value = TRUE)))
cat('Metadata dimensions:', dim(meta.project), '\n')

## read methlathin 
sparseCCA_result_dirname <- paste0(result.loc,"taxa_rppa/case")
cat("\nRead sparseCCA componet from",sparseCCA_result_dirname)
dirnames <- list.dirs(sparseCCA_result_dirname)
sig_components_dirname <- grep("sig_components",dirnames,value = TRUE)
# Get the components filename
filenames <- list.files(sig_components_dirname)
filenames <- grep("component",filenames,value = TRUE)
stopifnot(length(filenames) > 0) # make sure the input dir path is set correctly
cat("File contains:",filenames)
filename <- filenames[9]

# 
# for(filename in filenames){

cat(paste0("Enrichment for component: ", filename,"\n"));
fig_component_dir <- paste0(fig.loc,project,"/rppa_taxa/")
ifelse(!dir.exists(fig_component_dir), dir.create(fig_component_dir,recursive = TRUE), FALSE)

sparseCCA_component <- read.table(paste0(sig_components_dirname,"/",filename),sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
#correlatarion_plot <- plot_componet(sparseCCA_component,filename)

sparseCCA_genes <- sparseCCA_component[,c("rppa","rppa_coeff")]
sparseCCA_genes <- sparseCCA_genes %>% arrange(desc(abs(rppa_coeff)))
sparseCCA_genes <- na.omit(sparseCCA_genes)
# idx<- percent * length(sparseCCA_genes$gene_coeff)
# sparseCCA_genes <- sparseCCA_genes[1:idx,]
cat("Dim of sparseCCA genes:",dim(sparseCCA_genes),"\n")
genes_of_interest <- sparseCCA_genes$rppa

## annotation of meth data
annotation <- read.table(file="C:/Lab/project/pan-cancer-host-microbiome-associations/data/protein-meta.txt",sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
annotation$name <- rownames(annotation)
annotation_split <- annotation %>%
  separate_rows(gene, sep = ",")
annotation_split_unique <- distinct(annotation_split, gene, name, .keep_all = TRUE)

annotation_split_unique_interest <- subset(annotation_split_unique,name %in% sparseCCA_genes$rppa)
genes_of_interest <- unique(annotation_split_unique_interest$gene)


KEGG_enrich_symbol <- unique(genes_of_interest)
# cor.coef <- subset(sparseCCA_genes,abs(sparseCCA_genes$gene_coeff) >0.02)
cor.coef <- sparseCCA_component
cor.coef.taxa=cor.coef[,c("taxa","taxa_coeff"),drop=F]
cor.coef.taxa=na.omit(cor.coef.taxa)

cor.coef.gene=cor.coef[,c("rppa","rppa_coeff"),drop=F]

cor.coef.gene <- na.omit(cor.coef.gene)
rownames(cor.coef.gene)=cor.coef.gene$rppa

df_gene <- (feat.rppa)[,cor.coef.gene$rppa]

## load correaltion data
cor.spearman <- as.matrix(fread("C:/Lab/project/pan-cancer-host-microbiome-associations/data/spearman/HNSC.bacteria.protein.spearman.correlation.test.result.tsv"))

colnames(cor.spearman)[1] <- "rppa"
sparseCCA_taxas <- sparseCCA_component[,c("taxa","taxa_coeff")]
sparseCCA_taxas <- na.omit(sparseCCA_taxas)

intersect(test$rppa,sparseCCA_genes$rppa)

## get cca result 
cor_spearman_CCA <- subset(as.data.frame(cor.spearman), rppa %in% sparseCCA_genes$rppa)
cor_spearman_CCA <- subset(as.data.frame(cor_spearman_CCA), taxa %in% sparseCCA_taxas$taxa)
rm(cor.spearman)

# cor_spearman_gene_interest <- subset(as.data.frame(cor_spearman_CCA), rppa %in% KEGG_enrich_symbol)
cor_spearman_gene_interest <- cor_spearman_CCA

cor_spearman_gene_interest$rho <- as.numeric(cor_spearman_gene_interest$rho)
cor_spearman_gene_interest$rho_padj <-as.numeric( cor_spearman_gene_interest$rho_padj)



#regexpr("ABCC1",annotation$gene)
interest_whole_taxa <- c("Campylobacter","Terrabacter","Flammeovirga","Lachnoclostridium","Sutterella","Collimonas")

rho_cutoff <- 0.1
cor.data.single.spearman = cor_spearman_gene_interest 


cor.data <- cor.data.single.spearman
colnames(cor.data)[1:2] <-c("from","to") 
cor.data <- cor.data %>% filter(abs(rho)>rho_cutoff)
# if (is.null(cor.data )| dim(cor.data)[1]==0){
#   cat("Could not find correaltio >0.3")
#   pdf(paste0(fig_component_dir,project,"_",filename,"_network.pdf"),width = 24,height = 12)
#   print(correlatarion_plot)
#   dev.off()
#   next;
# }else{
#   head(cor.data)
# }
interest_whole_taxa <- c("Campylobacter","Terrabacter","Flammeovirga","Lachnoclostridium","Sutterella","Collimonas")
cor.data$to <- str_split_fixed(cor.data$to,'g__',2)[,2]

cor.data <- subset(cor.data,to %in% interest_whole_taxa)


if (is.null(cor.data )| dim(cor.data)[1]==0){
  cat("Could not find correaltio >0.3")
  pdf(paste0(fig_component_dir,project,"_",filename,"_network.pdf"),width = 24,height = 12)
  print(correlatarion_plot)
  dev.off()
  next;
}else{
  head(cor.data)
}

## 3.1.1 计算每个节点具有的链接数(degree)
c(as.character(cor.data$from),as.character(cor.data$to)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> cor.nodes
#加上coef

colnames(cor.nodes) <- c("name", "edges.num")
rownames(cor.nodes) <- cor.nodes$name


# 改变coef中细菌的名字
cor.coef.taxa.curated <- cor.coef.taxa
rownames(cor.coef.taxa.curated) <- cor.coef.taxa.curated$taxa
rownames(cor.coef.taxa.curated) <- str_split_fixed(rownames(cor.coef.taxa.curated),'g__',2)[,2]
colnames(cor.coef.taxa.curated) <- c("item","coef")
cor.coef.gene.curated <- cor.coef.gene
colnames(cor.coef.gene.curated) <- c("item","coef")

cor.coef.curated <- rbind(cor.coef.gene.curated,cor.coef.taxa.curated)
cor.nodes[,"coef"] <-cor.coef.curated[rownames(cor.nodes),"coef"]
cor.nodes[,"node_size"] <-abs(cor.nodes[,"coef"])
cor.nodes <- cor.nodes %>%
  mutate(
    subgroup = case_when(
      name %in% rownames(cor.coef.gene.curated) ~"rppa",
      name %in% rownames(cor.coef.taxa.curated) ~"taxa")
  )


cor.nodes.test <- cor.nodes
head(cor.nodes.test)
cor.nodes.test$"key_label" <- NA
#cor.nodes.test[which(abs(cor.nodes.test$coef)>0.04),"key_label"] <- cor.nodes.test[which(abs(cor.nodes.test$coef)>0.04),"name"]


cor.nodes.test[which(cor.nodes.test$name %in% interest_whole_taxa),"key_label"] <- cor.nodes.test[which(cor.nodes.test$name %in% interest_whole_taxa),"name"]

cat("\nCreate graph\n")
graph <- tbl_graph(nodes = cor.nodes.test, edges = cor.data, directed = FALSE)
graph
cat("Filter")
graph = graph %>%
  activate("edges") %>%
  filter(!edge_is_multiple()) %>%
  filter(!edge_is_loop())
graph


cat("\n Ploting+++++++++++++++++++++++++++++++++\n\n\n")
color = c("positive" ="#f96f71","negative" ="#9fc0e5")
Correlation <- guide_legend(title="Correlation type",
                            direction="vertical",
                            order=2,# 图例排列序号
                            ncol=1,
                            byrow=FALSE,
                            title.theme = element_text(
                              size = 14,
                              face = "bold",
                              colour = "black"),
                            label.theme = element_text(
                              size = 12,colour = "black")
)
width_legend <- guide_legend(title="|Spearman rho|",
                             direction="vertical",
                             order=1,# 图例排列序号
                             ncol=1,
                             byrow=FALSE,
                             title.theme = element_text(
                               size = 14,
                               face = "bold",
                               colour = "black"),
                             label.theme = element_text(
                               size = 12,colour = "black")
)

p_rsize <- ggraph(graph,layout ='dh') +
  # geom_edge_link(aes(edge_color =r
  #                    
  # ),edge_width=1.5) +
  #geom_edge_link(aes(edge_color = linecolor),edge_width=0.2)+
  geom_edge_link(edge_color="#969896",edge_width=0.2)+
  geom_node_point(aes(
                      shape=subgroup,
                      color=subgroup
                      ),size=7)+
  scale_shape_manual(values=c(18,17))+#which could both use fill and colour 21-25
  scale_color_manual(
    values = c("#efc883","#c8ea93")
  )+
  scale_size_continuous(name = "sparse CCA coefficients",breaks = seq(0.1,1,0.2),range = c(3,10))+
geom_node_text(aes(label = name), repel = TRUE,check_overlap=TRUE)+
  theme_graph()+
  labs(title = paste0(project,"_",filename))
# ggsave(p_rsize,file="test.pdf")

# 
# font_import()
# # This tries to autodetect the directory containing the TrueType fonts.
# # If it fails on your system, please let me know.
# # Vector of font family names
# fonts()
# 
# # Show entire table
# 
# fonttable()
# 
# loadfonts()

pdf(paste0(fig_component_dir,project,"_",filename,"_network",rho_cutoff,".pdf"),width = 12,height = 12)
# print(correlatarion_plot)
print(p_rsize)
dev.off()

cat("\nSucessfully save plot in ",paste0(fig_component_dir,project,"_",filename,"_network.pdf"),"\n")


rho_cutoff <- "Non"
write.table(cor.data,file=paste0(fig_component_dir,project,"_",filename,"_correlation_data_",rho_cutoff,".txt"),sep="\t", row.names = F)



## save nodes
cor_spearman_CCA$taxa <-  str_split_fixed((cor_spearman_CCA$taxa) ,'g__',2)[,2]
c(as.character(cor_spearman_CCA$rppa),as.character(cor_spearman_CCA$taxa)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> cor.nodes.all
cor.nodes.all
colnames(cor.nodes.all) <- c("name", "edges.num")
rownames(cor.nodes.all) <- cor.nodes.all$name


colnames(sparseCCA_genes) <- c("name","sCCA_coeff")
colnames(sparseCCA_taxas) <- c("name","sCCA_coeff")
sparseCCA_taxas$name <- str_split_fixed((sparseCCA_taxas$name) ,'g__',2)[,2]
### combinde data
sparseCCA_nodes <- rbind(sparseCCA_genes,sparseCCA_taxas)

sparseCCA_nodes <- na.omit(sparseCCA_nodes)

cor.nodes.all <- left_join(cor.nodes.all,sparseCCA_nodes,by="name")
cor.nodes.all  <- cor.nodes.all  %>%
  mutate(
    feature_type = case_when(
      name %in% sparseCCA_genes$name ~"protein",
      name %in% sparseCCA_taxas$name ~"taxa")
  )
cor.nodes.all.annotation <- cor.nodes.all

cor.nodes.all.annotation$edges.num <- NULL
cor.nodes.all.annotation$SparseCCA_Componet <- filename
cor.nodes.all.annotation$Cancer_type <- project

write.table(cor.nodes.all.annotation,file=paste0(fig_component_dir,project,"_",filename,"_nodes_data.csv"),sep=",", row.names = F)

