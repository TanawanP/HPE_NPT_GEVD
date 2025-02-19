# Title: R for Distribution-Based Estimation of Hyperparameters for the Negative Power Transformation of the GEV Distribution
# Modified version: 2025-02-19

getReturnLevels_npt <- function(f, a, rp=c(25,50,100), conf=0.95, sd_scale=NULL, trans=TRUE){
  # make empty list
  dn <- list(NULL,c('period','2.5% ','Estimate','97.5%','se'))
  rl <- matrix(NA,nrow=length(rp),ncol=5,dimnames=dn)
  rl[,1] <- rp
  
  # Extract the estimated parameters (location, scale, shape)
  params <- f$mle
  names(params) <- c("location", "scale", "shape")
  
  # Manually compute the covariance matrix using the Fisher Information
  Cov_mat <- f$cov
  if(any(is.na(Cov_mat))){
    resam_x <- sample(f$data, replace = TRUE)
    refit_x <- gev.fit(resam_x, show=FALSE)
    cov_matrix <- refit_x$cov
  }else{
    cov_matrix <- f$cov
  }
  
  # Define the return period and calculate the return level for the transformed data
  rl_transformed <- 
    params["location"] - (params["scale"] / params["shape"]) * (1 - (-log(1 - 1/rp))^(-params["shape"]))
  
  
  if(trans){
    rl_transformed_scaled <- (rl_transformed^(-1/a) - 1) * (sd(f$min.data) / sd_scale) + mean(f$min.data)
    rl[,3] <-  rl_transformed_scaled
  }else{
    rl[,3] <-  rl_transformed^(-1/a)
  }
  
  # Delta method: Calculate the gradient of the rl with respect to the parameters (for maxima data)
  gradient <- NULL
  for (l in 1:length(rp)) {
    location_grad <- 1
    scale_grad <- -((1 - (-log(1 - 1/rp[l]))^(-params["shape"])) / params["shape"])
    shape_grad <- (params["scale"] / params["shape"]^2) * 
      (1 - (-log(1 - 1/rp[l]))^(-params["shape"])) - 
      (params["scale"] / params["shape"]) * 
      (-log(1 - 1/rp[l]))^(-params["shape"]) * log(-log(1 - 1/rp[l]))
    
    # Combine the gradients into a numeric vector
    gradient <- rbind(gradient, c(location_grad, scale_grad, shape_grad)) 
  }
  
  # Compute the standard error of the return level for the transformed data (maxima)
  v <- gradient%*%cov_matrix%*%t(gradient)
  
  # Conf with rl
  rl[,2] <- rl[,3]-qnorm(1-(1-conf)/2)*sqrt(diag(v))
  rl[,4] <- rl[,3]+qnorm(1-(1-conf)/2)*sqrt(diag(v))
  rl[,5] <- sqrt(diag(v))
  
  return(rl)
}

boots.lme_par <- function(data, nboot=1000, method="l-moment"){
  # Initialize a matrix to store bootstrap parameters
  boot_params <- matrix(NA, nrow = nboot, ncol = 3)
  
  # Bootstrap loop
  for (i in 1:nboot) {
    sample_data <- sample_fit  <- NULL
    # Resample data with replacement
    sample_data <- sample(data, replace = TRUE, size = length(data))
    
    # Refit distribution
    if(method=="l-moment"){
      sample_fit <- tryCatch({
        fevd(sample_data, method="Lmoments", type="GEV")$results
      }, error = function(e) rep(NA, 3))
    }else if(method=="mle"){
      sample_fit <- tryCatch({
        gev.fit(sample_data, show=FALSE)$mle
      }, error = function(e) rep(NA, 3))
    }
    # Store fitted parameters
    boot_params[i, ] <- sample_fit
  }
  # Calculate standard errors for each parameter
  param_se <- apply(boot_params, 2, sd, na.rm = TRUE)
  return(param_se)
}

boots.lme_rl_npt <- function(data, a, min.data=NULL, nboot=1000, rp=c(25,50,100), conf=0.95, sd_scale=NULL, method0="Lmoments", trans=TRUE){
  # Set results matrix
  dn <- list(NULL,c('period','2.5% ','Estimate','97.5%','se'))
  rl <- matrix(NA,nrow=length(rp),ncol=5,dimnames=dn)
  rl[,1] <- rp
  
  if(trans){
    # True RL
    true_fit <- fevd(data, method=method0, type="GEV")
    true_rl <- tryCatch({ 
      ((return.level(true_fit, return.period=rp)^(-1/a) - 1) * (sd(min.data) / sd_scale) + mean(min.data)) %>% as.numeric() 
    }, error = function(e) rep(NA, length(rp)))
    
    # Initialize a matrix to store bootstrap parameters
    boot_params <- matrix(NA, nrow = nboot, ncol = length(rp))
    
    # Bootstrap loop
    for (i in 1:nboot) { 
      sample_data <- sample_fit <- sample_rl <- NULL
      # Resample data with replacement
      sample_data <- sample(data, replace = TRUE, size = length(data))
      
      # Refit distribution
      sample_fit <- fevd(sample_data, method=method0, type="GEV")
      
      # Recal RL
      sample_rl <- tryCatch({ 
        ((return.level(sample_fit, return.period=rp)^(-1/a) - 1) * (sd(min.data) / sd_scale) + mean(min.data)) %>% as.numeric() 
      }, error = function(e) rep(NA, length(rp)))
      
      # Store fitted parameters
      boot_params[i, ] <- sample_rl
    }
    
  }else{
    # True RL
    true_fit <- fevd(data, method=method0, type="GEV")
    true_rl <- tryCatch({ return.level(true_fit, return.period=rp)^(-1/a) %>% as.numeric() 
    }, error = function(e) rep(NA, length(rp)))
    
    # Initialize a matrix to store bootstrap parameters
    boot_params <- matrix(NA, nrow = nboot, ncol = length(rp))
    
    # Bootstrap loop
    for (i in 1:nboot) { 
      sample_data <- sample_fit <- sample_rl <- NULL
      # Resample data with replacement
      sample_data <- sample(data, replace = TRUE, size = length(data))
      
      # Refit distribution
      sample_fit <- fevd(sample_data, method=method0, type="GEV")
      
      # Recal RL
      sample_rl <- tryCatch({ return.level(sample_fit, return.period=rp)^(-1/a) %>% as.numeric() 
      }, error = function(e) rep(NA, length(rp)))
      
      # Store fitted parameters
      boot_params[i, ] <- sample_rl
    }
  }
  
  # Calculate standard errors for each parameter
  param_se <- apply(boot_params, 2, sd, na.rm = TRUE)
  rl[,3] <- true_rl
  rl[,2] <- rl[,3]-qnorm(1-(1-conf)/2)*param_se
  rl[,4] <- rl[,3]+qnorm(1-(1-conf)/2)*param_se
  rl[,5] <- param_se
  return(rl)
}

cal.f.cvm.lm <- function(a, sam.min.data, scale_min.data){
  # transform data
  trans.data <- sam.min.data^(-a)
  
  # gev distribution fit
  gev.data <- fevd(trans.data, method="Lmoments", type="GEV")
  params <- gev.data$results %>% as.numeric()
  
  # test gof
  gof.test <- goftest::cvm.test(scale_min.data,'pgev', params[1], params[2], params[3])
  return(gof.test$statistic)
}

cal.f.cvm.mle <- function(a, sam.min.data, scale_min.data){
  # transform data
  trans.data <- sam.min.data^(-a)
  
  # gev distribution fit
  # gev.data <- fevd(trans.data, method="MLE", type="GEV")
  # params <- gev.data$results$par %>% as.numeric()
  params <- gev.fit(trans.data, show=FALSE)$mle
  
  # test gof
  gof.test <- goftest::cvm.test(scale_min.data,'pgev', params[1], params[2], params[3])
  return(gof.test$statistic)
}

opt_a0_resam <- function(target_func=NULL, scale_min.data=NULL, a0=NULL, n0=NULL, per_n0=0.8, rep_n=50){
  opt.a <- NULL
  # rep_n <- 50
  rep_max <- 100
  r.vec <- NULL
  j <- 0
  count <- 0
  while ((j < rep_n) && (count < rep_max)) {
    red_n <- n0*per_n0
    resam.min.data <- sample(scale_min.data, red_n) # sampling generating for a distribution
    
    op <- NULL
    try({op <- optim(par=1, fn=get(paste0(target_func)), sam.min.data=resam.min.data, scale_min.data=scale_min.data, lower = 0, upper = 100, method="L-BFGS-B")}, silent = TRUE)
    
    if(!is.null(op) && (op$par>0) && (round(op$par,4)!=a0) &&
       (round(op$par,4)!=0) && (round(op$par,4)!=100)){
      j <- j+1
      r.vec <- c(r.vec, op$par)
      cat("rep::",j,"est.a::",op$par,"\n")
    }
    count <- count + 1
  }
  # attention here!!!!
  if(!is.null(r.vec)){
    if(length(r.vec)>2){
      opt.a = median(r.vec)
    }
  }
  return(list(opt.a=opt.a, r.opt.a=r.vec))
}

opt_a0_addnoise <- function(target_func=NULL, scale_min.data=NULL, a0=NULL, n0=NULL, rep_n=50){
  opt.a <- NULL
  # rep_n <- 50
  rep_max <- 100
  r.vec <- NULL
  j <- 0
  count <- 0
  while ((j < rep_n) && (count < rep_max)) {
    # red_n <- n0*per_n0
    resam.min.data <- scale_min.data + rnorm(n0, mean=0, sd=0.001)
    op <- NULL
    
    try({op <- optim(par=1, fn=get(paste0(target_func)), sam.min.data=resam.min.data, scale_min.data=scale_min.data, lower = 0, upper = 100, method="L-BFGS-B")}, silent = TRUE)
    
    if(!is.null(op) && (op$par>0) && (round(op$par,4)!=a0) &&
       (round(op$par,4)!=0) && (round(op$par,4)!=100)){
      j <- j+1
      r.vec <- c(r.vec, op$par)
      cat("rep::",j,"est.a::",op$par,"\n")
    }
    count <- count + 1
  }
  # attention here!!!!
  if(!is.null(r.vec)){
    if(length(r.vec)>2){
      opt.a = median(r.vec)
    }
  }
  return(list(opt.a=opt.a, r.opt.a=r.vec))
}

opt_a0_resam_cv <- function(target_func=NULL, sam.min.data=NULL, scale_min.data=NULL, a0=NULL, n0=NULL){
  k_folds <- 5
  opt.a <- NULL
  rep_max <- 10
  r.vec <- NULL
  count <- 0
  fold_results <- c()
  
  # Create k-fold indices for cross-validation
  folds <- cut(seq(1, length(sam.min.data)), breaks=k_folds, labels=FALSE)
  
  # Loop through each fold for cross-validation
  for (i in 1:k_folds) { # k_folds=1
    # Split the data into training and validation sets
    train_indices <- which(folds != i)
    train_data <- sam.min.data[train_indices]
    
    j <- 0
    fold_r.vec <- NULL
    invalid <- TRUE
    while (((j < 1) && invalid) || ((j < 1) && (count < rep_max))) {
      op <- NULL
      try({op <- optim(par=a0, sam.min.data=train_data, scale_min.data=scale_min.data, fn=get(paste0(target_func)), method="BFGS")}, silent = TRUE)
      
      if(!is.null(op) && (op$par > 0)){
        j <- j + 1
        invalid <- FALSE
        fold_r.vec <- op$par
        cat("Fold::", i, " est.a::", op$par, "\n")
      }else{
        train_data <- sample(train_data, size=length(train_data)*0.9, replace=F) #TRUE or FALSE
      }
      count <- count + 1
    }
    fold_results <- c(fold_results, fold_r.vec)
  }
  # Return the average optimization result across folds
  if (length(fold_results) > 0) {
    opt.a <- mean(fold_results)
  }
  return(list(opt.a=opt.a, r.opt.a=fold_results))
}
