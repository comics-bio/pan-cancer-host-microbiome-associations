library(glmnet)
library(hdi)
library(stabs)
library(doMC)
library(car)

## Input
args = commandArgs(trailingOnly=TRUE)

project_code <- args[1]
cell <- read.table(args[2], sep='\t', row.names = 1, header = T, check.names = F)
microbe <- read.table(args[3], sep='\t', row.names = 1, header = T, check.names = F)
nthread <- args[4]

registerDoMC(cores = nthread)

cell <- t(cell)
short_row_names <- substr(rownames(microbe), 1, 15)
microbe <- microbe[!duplicated(short_row_names), ]  # remove duplication randomly
rownames(microbe) <- substr(rownames(microbe), 1, 15)
keep <- intersect(rownames(microbe), rownames(cell))

x <- as.matrix(microbe[keep,]) # sample x taxa
y <- as.matrix(cell[keep,]) # sample x cell_type

## lasso
fit.cv.lasso <- function(x, y_i, kfold, repeats){
  
  lambdas = NULL
  r.sqr.final <- numeric()
  r.sqr.final.adj <- numeric()
  r.sqr.CV.test <- numeric()
  lasso.cv.list <- list()
  
  for (i in 1:repeats){
    
    ## glmnet CV
    fit <- cv.glmnet(x, y_i, alpha=1, nfolds=kfold, type.measure = "mse",
                     keep = TRUE, grouped=FALSE, parallel = TRUE)  
    errors = data.frame(fit$lambda,fit$cvm)
    lambdas <- rbind(lambdas,errors)
    
    ## Get R^2 of final model
    r.sqr.final[i] <- r_squared(as.vector(y_i), 
                                as.vector(predict(fit$glmnet.fit, 
                                                  newx = x, s = fit$lambda.min)))
    ## Get adjusted R^2
    r.sqr.final.adj[i] <- adj_r_squared(r.sqr.final[i], n = nrow(x), 
                                        p = sum(as.vector(coef(fit$glmnet.fit, 
                                                               s = fit$lambda.min)) > 0))
    
  }
  
  # take mean cvm for each lambda
  lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)
  # dim(lambdas)
  # select the best one
  bestindex = which(lambdas[2]==min(lambdas[2]))
  bestlambda = lambdas[bestindex,1]
  
  return(list(bestlambda = bestlambda, r.sqr = median(r.sqr.final), 
              r.sqr.adj = median(r.sqr.final.adj)
  ))
}
## functions to compute R2
## Adapted from https://rpubs.com/kaz_yos/alasso
r_squared <- function(y, yhat) {
  ybar <- mean(y)
  ## Total SS
  ss_tot <- sum((y - ybar)^2)
  ## Residual SS
  ss_res <- sum((y - yhat)^2)
  ## R^2 = 1 - ss_res/ ss_tot
  1 - (ss_res / ss_tot)
}
## Function for Adjusted R^2
## n sample size, p number of prameters
adj_r_squared <- function(r_squared, n, p) {
  1 - (1 - r_squared) * (n - 1) / (n - p - 1)
}

## Inspired from the details for estimateSigma() function in selectiveInference package: https://cran.r-project.org/web/packages/selectiveInference/selectiveInference.pdf
estimate.sigma.loocv <- function(x, y_i, bestlambda, tol) {
  
  lasso.fit = glmnet(x,y_i,alpha = 1)
  # beta <- coef(lasso.fit, s = bestlambda)[-1] ## This gives coefficients of fitted model, not predicted coeff. 
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
    sigma.flag = 0
  } else{
    sigma = 1 ## conservative option
    # sigma = ss_res ## lenient option
    sigma.flag = 2
  }
  
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
}


## Start
Z = NULL
lasso_res <- foreach (i = 1:ncol(y), .combine=rbind) %do% {
  y_i <- as.vector(logit(y[,i])) # transform percentage y

  ## Find best lambda with loocv
  fit.model <- fit.cv.lasso(x, y_i,  kfold = length(y_i), repeats = 1)
  bestlambda <- fit.model$bestlambda
  r.sqr <- fit.model$r.sqr
  r.sqr.adj <- fit.model$r.sqr.adj
  
  ## Estimate sigma using the estimated lambda param
  sigma.myfun <- estimate.sigma.loocv(x, y_i, bestlambda, tol=1e-4)
  sigma <- sigma.myfun$sigmahat
  beta <- as.vector(sigma.myfun$betahat)[-1]
  sigma.flag <- sigma.myfun$sigmaflag
  
  ## Inference
  lasso.proj.fit <- lasso.proj(x, y_i, multiplecorr.method = "BH",
                               betainit = beta,
                               parallel = T,
                               ncores = nthread,
                               sigma = sigma,
                               suppress.grouptesting = T, 
                               Z=Z,
                               return.Z = T)
  if (is.null(Z)) Z <- lasso.proj.fit$Z
  
  lasso.ci <- as.data.frame(confint(lasso.proj.fit, level = 0.95))
  
  ## perform stability selection using glmnet lasso
  stab.glmnet <- stabsel(x = x, y = y_i,
                         fitfun = glmnet.lasso, cutoff = 0.6,
                         PFER = 1)
  
  data.frame(y = rep(colnames(y)[i], length(lasso.proj.fit$pval)), 
             x = names(lasso.proj.fit$pval.corr), 
             r.sqr = r.sqr, r.sqr.adj = r.sqr.adj,
             pval = lasso.proj.fit$pval, padj = lasso.proj.fit$pval.corr, 
             ci.lower = lasso.ci$lower, ci.upper = lasso.ci$upper,
             sigma = sigma, sigma.flag = sigma.flag,
             stabsel_prob = stab.glmnet$max,
             stabsel_pass = names(stab.glmnet$max) %in% names(stab.glmnet$selected),
             row.names=NULL)
}

write.table(lasso_res, paste0(project_code, '.lasso.projection.test.result.tsv'), sep='\t', row.names = F, quote=F)

