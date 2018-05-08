## load used packages
library(cvTools)
library(ggplot2)
library(boot)
library(glmnet)
library(MASS)

## load DATA, get N, SE, p, weights, !r2, and dummies? Yase, OK

lightart <- read.csv("Lightart.csv", header = TRUE)

N <- nrow(lightart)




##
slcs <- c(
    "Imputed",
    "MAF",
    "FEMALE",
    "AGE",
    "PLAT",
    "SOFT"
    )


catSLCs <- c(
    "FEMALE",
    "PLAT",
    "SOFT"
    )

sum(unlist(lapply(lightart[,catSLCs], nlevels))) - length(catSLCs)
 

# SEs


lightartSEs <- vector(length = N)


lightartWTs <- (1/lightart$SE^2)/sum(1/lightart$SE^2)



##############
hist(lightartWTs)




simple_reg <- function(X, Y, wts){
    P <- ncol(X)
    fixselect <- vector(mode = "numeric", length = P)
    coefs <- vector(mode = "numeric", length = P)
    r2 <- vector(mode = "numeric", length = P)
    
    for (j in 1: P){
        templm <- lm(Y~X[,j], weights = wts)
        wmp<- which.min(summary(templm)$coefficients[-1,"Pr(>|t|)"]) + 1

        if(is.numeric(wmp)&&length(wmp) ==0L){
            fixselect[j] <- 0
            coefs[j] <- NA
            r2[j] <- NA
            next
        }
        tempp <- summary(templm)$coefficients[wmp,"Pr(>|t|)"]
        fixselect[j] <- as.numeric(tempp <= 0.05)
        coefs[j] <- summary(templm)$coefficients[wmp,"Estimate"]
        r2[j] <- summary(templm)$r.squared
    }

    output <- list(select = fixselect, coefs = coefs, r2 = r2)
    return(output)
}



lasso_reg <- function(X, Y, wts){

    lassomodel <- glmnet(x = X, y = Y, weights = wts)
    return(output)

}









boot_rep_simp <- function(data, idx, wts){

    replicate <- simple_reg(X = data[idx,-1], Y = data[idx,1], wts = wts[idx])
    
    return(replicate$coefs)
}


boot_rep_las <- function(data, idx, wts){

    replicate <- lasso_reg(X = data[idx,-1], Y = data[idx,1], wts = wts[idx])
    
    return(replicate$coefs)

}




lightart <- data.matrix(lightart)



light1 <- simple_reg(X = lightart[,slcs], Y = lightart[,'Beta'], wts = lightartWTs)


llight1 <- lasso_reg(X = lightart[,slcs], Y = lightart[,"Beta"],wts = lightartWTs)

bdataM1 <- matrix(cbind(lightart[,'Beta'], lightart[,slcs[-2]]), nrow = N, ncol = length(slcs) )




bootsimpM1 <- boot(bdataM1, boot_rep_simp, sim = btype, R = 1000, wts = lightartWTs)
bsimpmat1 <- matrix(nrow = (ncol(bdataM1) - 1), ncol = 3)
for (j in 1:nrow(bsimpmat1)){
    bsimpmat1[j,1:2] <- boot.ci(bootsimpM1,index = j, type = "perc")$percent[4:5]
    bsimpmat1[j,3] <- (sign(bsimpmat1[j,1])*sign(bsimpmat1[j,2]) > 0)
}



bootlas1 <- boot(bldataM1, boot_rep_las, R = 1000, wts = lightartWTs[CC1])
blasmat1 <- matrix( nrow = (ncol(bldataM1) -1 ), ncol = 3)
for (j in 1:nrow(blasmat1)){
    blasmat1[j,1:2] <- boot.ci(bootlas1, index = j, type = "perc")$percent[4:5]
    blasmat1[j,3] <- (sign(blasmat1[j,1])*sign(blasmat1[j,2]) > 0)
}


