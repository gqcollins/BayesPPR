bppr_pca <- function(X, Y, n_pc = NULL, prop_var = 0.99, n_cores = 1, par_type = 'fork', n_ridge_mean = 10, n_ridge_max = NULL, n_act_max = NULL, df_spline = 4, prob_relu = 2/3, shape_var_coefs = NULL, rate_var_coefs = NULL, n_dat_min = NULL, scale_proj_dir_prop = NULL, w_n_act_init = NULL, w_feat_init = NULL, n_post = 1000, n_burn = 9000, n_thin = 1, print_every = 0, model = 'bppr'){

  pca_Y <- pca_setup(X, Y, n_pc = n_pc, prop_var = prop_var)

  n_cores_max <- parallel::detectCores()

  if(n_cores > n_cores_max){
    warning(paste0("Specified n_cores = ", n_cores, ". Proceeding with n_cores = min(n_cores, n_pc, detectCores()) = ",
                   min(n_cores, n_pc, n_cores_max)))
  }
  n_cores <- min(n_cores, n_pc, n_cores_max)

  run_bppr <- parse(text =
  "bppr(X, pca_Y$Y_new[, i], n_ridge_mean = n_ridge_mean, n_ridge_max = n_ridge_max,
           n_act_max = n_act_max, df_spline = df_spline, prob_relu = prob_relu,
           shape_var_coefs = shape_var_coefs, rate_var_coefs = rate_var_coefs,
           n_dat_min = n_dat_min, scale_proj_dir_prop = scale_proj_dir_prop,
           w_n_act_init = w_n_act_init, w_feat_init = w_feat_init,
           n_post = n_post, n_burn = n_burn, n_thin = n_thin, print_every = print_every,
           model = model)"
  )

  if (n_cores == 1){
    fit_list <- lapply(1:n_pc, function(i) eval(run_bppr))
  }else if(par_type == "socket"){
    cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential")
    fit_list <- parallel::parLapply(cl, 1:n_pc, function(i) eval(run_bppr))
    parallel::stopCluster(cl)
  }else if(par_type == "fork"){
    fit_list <- parallel::mclapply(1:pca_Y$n_pc, function(i) eval(run_bppr),
      mc.cores = n_cores, mc.preschedule = FALSE)
  }

  out <- list(pca_Y = pca_Y, fit_list = fit_list)
  class(out) <- "bppr_pca"
  return(out)
}
