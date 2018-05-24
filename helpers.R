# get active set from lasso fit ----
activeset_names <- function(fit, lambda_pos, X){
  colnames(X)[
    unlist(predict(fit, type = "nonzero", s = fit$lambda[lambda_pos]))
      ]
 }


# lasso fitted values function ----
glmnet_fitted <- function(fit, lambda_pos, X){
  predict(fit, newx = X, s = fit$lambda[lambda_pos])
}

# lasso residuals function ----
glmnet_resids <- function(fit, lambda_pos, X,Y){
  Y - glmnet_fitted(fit, lambda_pos, X)
}


# Make summary of variables, add column of correlation with current residual ----
glmnet_summextra <- function(fit, lambda_pos, X,Y){
  active_set <- activeset_names(fit, lambda_pos, X)

  current_residuals <- glmnet_resids(fit, lambda_pos, X,Y)

  active_corrs <- cor(current_residuals, X[,active_set])

  other_corrs <- cor(current_residuals, X[, setdiff(colnames(X),active_set)])

  next_to_join <- which.max(abs(other_corrs))

  next_name <- colnames(other_corrs)[[next_to_join]]
  extra_summ <- data.frame(
    psych::describe(X[,c(active_set, next_name)]),
    resid.cor = c(active_corrs,other_corrs[[next_to_join]])
    )

  extra_summ$vars <- c(active_set, next_name)
  # Parallel coordinates plot of active set plus next variable ----
  pcp_next <- MASS::parcoord(X[,c(active_set, next_name)], col =  1 + (0:nrow(X))%/%50)

return(list(extra_summary = extra_summ, pcp_extra = pcp_next))
  }


# Active set PC plot ----
active_pcs <- function(X, active_set, study_names){
if(length(active_set) < 2){
  # plot first 2 pcs
  plot(prcomp(X, retx = TRUE)$x[,1:2] )
  text(prcomp(X, retx = TRUE)$x[,1:2] , labels=study_names)

} else {
  # plot
  plot(prcomp(X[,active_set], retx = TRUE)$x[,1:2] )
  text(prcomp(X[, active_set], retx = TRUE)$x[,1:2] , labels=study_names)
}
}




