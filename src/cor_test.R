library(data.table)

args = commandArgs(trailingOnly=TRUE)

project_code <- args[1]
cell <- read.table(args[2], sep='\t', row.names = 1, header = T, check.names = F)
microbe <- read.table(args[3], sep='\t', row.names = 1, header = T, check.names = F)

cell <- t(cell)
short_row_names <- substr(rownames(microbe), 1, 15)
microbe <- microbe[!duplicated(short_row_names), ]  # remove duplication randomly
rownames(microbe) <- substr(rownames(microbe), 1, 15)
keep <- intersect(rownames(microbe), rownames(cell))

microbe <- microbe[keep,]
cell <- cell[keep,]

pair_wise_res <- as.data.table(expand.grid(colnames(cell), colnames(microbe)))
colnames(pair_wise_res) <- c('cell_type', 'taxa')
pair_wise_res[, c("rho", "rho_p", "tau", "tau_p"):=list(numeric(), numeric(), numeric(), numeric())]
setkey(pair_wise_res, cell_type, taxa)

for (microbe_idx in 1:ncol(microbe)) {
  for (cell_idx in 1:ncol(cell)) {
    cell_name <- colnames(cell)[cell_idx]
    microbe_type <- colnames(microbe)[microbe_idx]
    rho_test_res <- cor.test(microbe[,microbe_idx], cell[,cell_idx],
                             method='spearman', continuity=T)
    tau_test_res <- cor.test(microbe[,microbe_idx], cell[,cell_idx],
                             method='kendall', continuity=T)
    pair_wise_res[.(cell_name, microbe_type), rho_p:=rho_test_res$p.value]
    pair_wise_res[.(cell_name, microbe_type), rho:=rho_test_res$estimate]
    pair_wise_res[.(cell_name, microbe_type), tau_p:=tau_test_res$p.value]
    pair_wise_res[.(cell_name, microbe_type), tau:=tau_test_res$estimate]
  }
}

pair_wise_res[, rho_padj:=p.adjust(rho_p, 'BH')]
pair_wise_res[, tau_padj:=p.adjust(tau_p, 'BH')]
write.table(pair_wise_res, paste0(project_code, '.correlation.test.result.tsv'), sep='\t', row.names = F, quote=F)

