library(survival)
library(survminer)

args = commandArgs(trailingOnly=TRUE)

cancer <- args[1] # cancer name
data <- read.table(args[2], sep=',', row.names = 1, header = T, check.names = F) # survival group.csv
name <- read.table(args[3]) # header.txt


for (g in 1418:2828){
    fit <- survfit(Surv(OS_time, OS) ~ data[,g], data=data)
    p <- ggsurvplot(fit, data=data,
                    pval = TRUE, conf.int = TRUE,
                    risk.table = TRUE, 
                    risk.table.col = "strata", 
                    linetype = "strata", 
                    surv.median.line = "hv", 
                    ggtheme = theme_bw(),
                    palette = c("#E7B800", "#2E9FDF"))
    filename=(paste0('OS',"-",cancer,"-", name[g,], '.pdf'))
    pdf(file=filename, width=6, height=6, onefile=FALSE)
    print(p)
    dev.off()

    fit.2 <- survfit(Surv(PFS_time, PFS) ~ data[,g], data=data)
    p.2 <- ggsurvplot(fit.2, data=data,
                    pval = TRUE, conf.int = TRUE,
                    risk.table = TRUE, 
                    risk.table.col = "strata", 
                    linetype = "strata", 
                    surv.median.line = "hv", 
                    ggtheme = theme_bw(),
                    palette = c("#E7B800", "#2E9FDF"))
    filename=(paste0('PFS',"-",cancer,"-", name[g,], '.pdf'))
    pdf(file=filename, width=6, height=6, onefile=FALSE)
    print(p.2)
    dev.off()
  }
dev.off()