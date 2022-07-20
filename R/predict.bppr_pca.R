predict.bppr_pca <-function(object, newdata, idx_use = NULL, n_cores = 1, par_type = 'fork', ...){

  n_keep <- object$fit_list[[1]]$n_keep
  if(is.null(idx_use)){
    idx_use <- 1:n_keep
  }else if(max(idx_use) > n_keep){
    stop("invalid 'idx_use'")
  }
  n_use <- length(idx_use)

  n <- nrow(newdata)
  D <- length(object$pca_Y$mn_Y)

  n_pc <- object$pca_Y$n_pc

  run_predict <- parse(text = 'predict(object$mod.list[[i]], newdata, idx_use = idx_use, ...)')
  run_pca_reverse <- parse(text = 'pca_reverse(matrix(preds_Y_new[i, , ],
                                                      ncol = object$pca_Y$n_pc, nrow = nrow(newdata)),
                                               object$dat)')

  if(n.cores == 1){
    preds_Y_new <- array(
      unlist(lapply(1:n_pc, function(i) eval(run_predict))),
      dim = c(n_use, n, n_pc))
    out <- array(
      unlist(lapply(1:n_use, function(i) eval(run_pca_reverse))),
      dim = c(D, n, n_use))
  }else if(parType == "socket"){
    cl <- parallel::makeCluster(min(n_cores, n_pc, parallel::detectCores()),
                                setup_strategy = "sequential")
    parallel::clusterExport(cl, varlist = c("newdata"), envir = environment())
    preds_Y_new <- array(
      unlist(parallel::parLapply(cl, 1:n_pc, function(i) eval(run_predict))),
      dim = c(n_use, n, n_pc))
    out <- array(
      unlist(parallel::parLapply(cl, 1:n_use, function(i) eval(run_pca_reverse))),
      dim = c(D, n, n_use))
    parallel::stopCluster(cl)
  }
  else if (parType == "fork") {
    preds_Y_new <- array(
      unlist(parallel::mclapply(1:n_pc, function(i) eval(run_predict),
                                mc.cores = n_cores)),
      dim = c(n_use, n, n_pc))
    out <- array(
      unlist(parallel::mclapply(1:n_use, function(i) eval(run_pca_reverse),
             mc.cores = n_cores)),
      dim = c(D, n, n_use))
  }

  out <- aperm(out, c(3, 2, 1))
  return(out)
}
