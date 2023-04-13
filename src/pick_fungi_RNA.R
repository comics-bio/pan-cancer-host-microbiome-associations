library(data.table)

args = commandArgs(trailingOnly=TRUE)

project_code <- args[1]
microbe <- read.table(args[2], sep='\t', row.names = 1, header = T, check.names = F)

short_row_names <- substr(rownames(microbe), 1, 20)
microbe <- microbe[!duplicated(short_row_names), ]  # remove duplication randomly
rownames.microbe <- substr(rownames(microbe), 1, 20)

RNA.microbe.name <- c()
DNA.microbe.name <- c()

for (rownames in rownames.microbe) {
  last.cha <- substr(rownames,20,20)
  if (last.cha == "R")
    RNA.microbe.name <- append(RNA.microbe.name,rownames)
  if (last.cha == "D")
    DNA.microbe.name <- append(DNA.microbe.name,rownames)
    
}

RNA.microbe <- microbe[RNA.microbe.name,]
DNA.microbe <- microbe[DNA.microbe.name,]

write.table(RNA.microbe, paste0(project_code, ".fungi.RNA.tsv"), sep='\t')
write.table(DNA.microbe, paste0(project_code, ".fungi.DNA.tsv"), sep='\t')