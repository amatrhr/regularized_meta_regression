## load used packages
library(cvTools)
library(ggplot2)
library(boot)
library(glmnet)
library(MASS)

## load DATA, get N, SE, p, weights, !r2, and dummies? Yase, OK

murphy <- read.csv("Murphy.csv", header = TRUE)

N <- nrow(murphy)

murphy$time <- murphy$Year + murphy$Month/12


##
slcs <- c(
    "pubs",
    "time",
    "n",
    "Stimuli",
    "Emotion",
    "Presentation",
    "Task",
    "Contrast",
    "At.Rest",
    "Hemisphere",
    "Method",
    "Sample",
    "Ancestry",
    "Age",
    "MF",
    "HWE"
    )


catSLCs <- c(
    "Stimuli",
    "Emotion",
    "Presentation",
    "Task",
    "Contrast",
    "At.Rest",
    "Hemisphere",
    "Method",
    "Sample",
    "Ancestry",
    "HWE"
    )

sum(unlist(lapply(murphy[,catSLCs], nlevels))) - length(catSLCs)
 
murphy$pubs <- (!is.na(murphy$time))
MPS <-  murphy[which(murphy$pubs == TRUE),]

# SEs

gSE <- function(m, effect){
    cm <- gamma(m/2)/(sqrt(m/2)*gamma((m-1)/2))
    ene <- m/2
    a <- m*cm^2/(m-2) 
    var <- (a/ene)*(1 + ene*effect^2) - effect^2
    se <- sqrt(var)

 
    return(se)
}

murphySEs <- vector(length = N)

for(j in 1:N){
    murphySEs[j] <- gSE(m = murphy$n[j], effect = murphy$Effect[j])
}

MPSEs <- murphySEs[murphy$pubs]

murphyWTs <- (1/murphySEs)/sum(1/murphySEs)
MPWTs <- (1/MPSEs)/sum(1/MPSEs)


##############
hist(murphyWTs)
hist(MPWTs)



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

    lassomodel <- cv.glmnet(x = X, y = Y, weights = wts, nfolds =4)
    coefs <- coefficients(lassomodel, s = lassomodel$lambda.min)[-1]
    lasselect <- as.numeric(coefs != 0)
    output <- list(select = lasselect, coefs = coefs)
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




murphy <- data.matrix(murphy)
MPS <- data.matrix(MPS)


murph1 <- simple_reg(X = murphy[,slcs], Y = murphy[,'Effect'], wts = murphyWTs)
murph2 <- simple_reg(X = MPS[,slcs[-1]], Y = MPS[,'Effect'], wts = MPWTs)

CC0 <- complete.cases(murphy[,slcs])

lmurph1 <- lasso_reg(X = murphy[CC0,slcs], Y = murphy[CC0,"Effect"],wts = murphyWTs[CC0])
CC00 <- complete.cases(MPS[,slcs[-1]])
lmurph2 <- lasso_reg(X = MPS[CC00,slcs[-1]], Y = MPS[CC00,"Effect"],wts = MPWTs[CC00])


        bdataM1 <- matrix(cbind(murphy[,'Effect'], murphy[,slcs[-2]]), nrow = N, ncol = length(slcs) )


        btype <- "ordinary"

        bootsimpM1 <- boot(bdataM1, boot_rep_simp, sim = btype, R = 1000, wts = murphyWTs)
        bsimpmat1 <- matrix(nrow = (ncol(bdataM1) - 1), ncol = 3)
        for (j in 1:nrow(bsimpmat1)){
            bsimpmat1[j,1:2] <- boot.ci(bootsimpM1,index = j, type = "perc")$percent[4:5]
            bsimpmat1[j,3] <- (sign(bsimpmat1[j,1])*sign(bsimpmat1[j,2]) > 0)
        }

CC1 <- complete.cases(bdataM1)
bldataM1 <- bdataM1[CC1,]
bootlas1 <- boot(bldataM1, boot_rep_las, R = 1000, wts = murphyWTs[CC1])
        blasmat1 <- matrix( nrow = (ncol(bldataM1) -1 ), ncol = 3)
        for (j in 1:nrow(blasmat1)){
            blasmat1[j,1:2] <- boot.ci(bootlas1, index = j, type = "perc")$percent[4:5]
            blasmat1[j,3] <- (sign(blasmat1[j,1])*sign(blasmat1[j,2]) > 0)
        }


bdataM2 <- matrix(cbind(MPS[,'Effect'], MPS[,slcs[-1]]), nrow = 29, ncol = length(slcs)  )


bootsimpM2 <- boot(bdataM2, boot_rep_simp, sim = btype, R = 1000, wts = MPWTs)
bsimpmat2 <- matrix(nrow = (ncol(bdataM2) - 1), ncol = 3)
for (j in 1:nrow(bsimpmat2)){
    bsimpmat2[j,1:2] <- boot.ci(bootsimpM2,index = j, type = "perc")$percent[4:5]
    bsimpmat2[j,3] <- (sign(bsimpmat2[j,1])*sign(bsimpmat2[j,2]) > 0)
}

CC2 <- complete.cases(bdataM2)
bldataM2 <- bdataM2[CC2,]
bootlas2 <- boot(bldataM2, boot_rep_las, R = 1000, wts = MPWTs[CC2])
blasmat2 <- matrix( nrow = (ncol(bldataM2) -1 ), ncol = 3)
for (j in 1:nrow(blasmat2)){
    blasmat2[j,1:2] <- boot.ci(bootlas2, index = j, type = "perc")$percent[4:5]
    blasmat2[j,3] <- (sign(blasmat2[j,1])*sign(blasmat2[j,2]) > 0)
}



