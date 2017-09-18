# aggregation script for dissertation simulations

# Joins all tables resulting from the run script, produces the 108 x 13 total table of results
# And analyses on each margin

# 15 sets of results, each with 10800 rows
# merge together into a giant data frame
# calculate summary and marginal stats
# present summary and marginal stats in nice .tex form  

Ns <- c(6, 18, 54)
ps <- c(5, 10, 20)
rs <- c(0.2, 0.5)
wts <- c('Eq', 'Skw')
model <- c('Fix', 'Rand', 'Hetero')
method <- c("simp", "las", "scad", "bsimp", "blas")
cells_table <- expand.grid(Ns, ps, rs, wts, model)
resultsmatrix <- matrix(ncol = 14)


for (i in 1:nrow(cells_table)){
       tempnamelist <- system("ls -1 ds_t*_s*.txt", intern = TRUE)

       for (name in tempnamelist){
           tempresults <- read.table(name, header = FALSE, colClasses= c(rep("numeric",3),"character","character", "numeric","numeric", "character", rep("numeric", 6)))
           resultsmatrix <- rbind(resultsmatrix, tempresults)
       }
     
}
colnames(resultsmatrix) <- c("N", "p", "R2", "wts", "model", "seed","achR","meth","vsp", "tpr", "fpr", "hvsp", "htpr", "hfpr")
 
analysis_table <- expand.grid(Ns, ps, rs, wts, model, method)
analyses_matrix <- cbind(analysis_table, matrix(nrow = nrow(analysis_table), ncol = 14))
for (q in 1:nrow(analysis_table)){
  W1 <- which(resultsmatrix[,'N'] == analysis_table[q,1])
  W2 <- which(resultsmatrix[,'p'] == analysis_table[q,2])
  W3 <- which(resultsmatrix[,'R2'] == analysis_table[q,3])
  W4 <- which(resultsmatrix[,'wts'] == as.numeric(analysis_table[q,4]))
  W5 <- which(resultsmatrix[,'model'] == as.numeric(analysis_table[q,5]))
  W6 <- which(resultsmatrix[,'meth'] == analysis_table[q,6])

  indices <- intersect(W1, intersect( W2, intersect( W3, intersect(W4, intersect(W5, W6)))))
  meanvec <- apply(resultsmatrix[indices, c(7,9:14)], 2, mean, na.rm = TRUE)
  sdvec <- apply(resultsmatrix[indices, c(7,9:14)], 2, sd, na.rm = TRUE)
  analyses_matrix[q, 7:20] <- c(unlist(meanvec), unlist(sdvec))
}
colnames(analyses_matrix) <- c("N", "p", "R2", "wts", "model","meth", "achR","vsp", "tpr", "fpr", "hvsp", "htpr", "hfpr","sachR","svsp", "stpr", "sfpr", "shvsp", "shtpr", "shfpr")
write.table(analyses_matrix,"analyses_matrix.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)"
