# Title: R for Distribution-Based Estimation of Hyperparameters for the Negative Power Transformation of the GEV Distribution
# Modified version: 2025-02-19

# Set path ---
setwd("/Users/thanawanp/workspace/5.InvGEV extension/Program")

# Packages ---
library(SpatialExtremes)
library(ismev)
library(extRemes)
library(evd)
library(lmom)
library(goftest)
library(dplyr)
library(stringr)
library(openxlsx)

# Source required functions ---
source("RequiredFunction.R")

# Main estimation function ---
est_funct <- function(db_method, min.data, rp0, conf0){
  # Rescale data | mu->1, sd->0.1
  trans0 <- TRUE
  sd_scale <- 0.1
  scale_min.data <- (min.data - mean(min.data)) * (sd_scale / sd(min.data)) + 1
  
  # Hyperparameter estimation--
  if(db_method == "resam"){ # 1) Resampling method
    # MLE---
    cal_mle_a <- NULL
    cal_mle_a <- opt_a0_resam(target_func="cal.f.cvm.mle", scale_min.data=scale_min.data, a0=1, n0=length(scale_min.data))
      
    if(is.null(cal_mle_a$opt.a)){
      mle_a_using_lm <- TRUE
      recal_mle_a <- NULL
      recal_mle_a <- opt_a0_resam(target_func="cal.f.cvm.lm", scale_min.data=scale_min.data, a0=1, n0=length(scale_min.data))
      
      if(is.null(recal_mle_a$opt.a)){
        recal_mle_a2 <- NULL
        recal_mle_a2 <- opt_a0_resam(target_func="cal.f.cvm.lm", 
                                     scale_min.data=sample(scale_min.data, size=(length(scale_min.data)*0.8), replace = TRUE), 
                                     a0=1, n0=length(scale_min.data)) 
        if(is.null(recal_mle_a2$opt.a)){
          mle_a <- 1
          r_mle_a <- NA
        }else{
          mle_a <- recal_mle_a2$opt.a
          r_mle_a <- recal_mle_a2$r.opt.a
        }
      }else{
        mle_a <- recal_mle_a$opt.a
        r_mle_a <- recal_mle_a$r.opt.a
      }
    }else{
      mle_a <- cal_mle_a$opt.a
      r_mle_a <- cal_mle_a$r.opt.a
      mle_a_using_lm <- FALSE
    }
    # LME---
    if(mle_a_using_lm){
      lm_a <- mle_a
      r_lm_a <- r_mle_a
    }else{
      recal_lm_a <- NULL
      recal_lm_a <- opt_a0_resam(target_func="cal.f.cvm.lm", scale_min.data=scale_min.data, a0=1, n0=length(scale_min.data))
      
      if(is.null(recal_lm_a$opt.a)){
        recal_lm_a2 <- NULL
        recal_lm_a2 <- opt_a0_resam(target_func="cal.f.cvm.lm",
                                    scale_min.data=sample(scale_min.data, size=(length(scale_min.data)*0.8), replace = TRUE), 
                                    a0=1, n0=length(scale_min.data)) 
        if(is.null(recal_lm_a2$opt.a)){
          lm_a <- 1
          r_lm_a <- NA
        }else{
          lm_a <- recal_lm_a2$opt.a
          r_lm_a <- recal_lm_a2$r.opt.a
        }
      }else{
        lm_a <- recal_lm_a$opt.a
        r_lm_a <- recal_lm_a$r.opt.a
      }
    }
  }else if(db_method == "adnoise"){ # 2) Adding small noise (AddNoise) method
    # MLE---
    cal_mle_a <- NULL
    cal_mle_a <- opt_a0_addnoise(target_func="cal.f.cvm.mle", scale_min.data=scale_min.data, a0=1, n0=length(scale_min.data))
    
    if(is.null(cal_mle_a$opt.a)){
      mle_a_using_lm <- TRUE
      recal_mle_a <- NULL
      recal_mle_a <- opt_a0_addnoise(target_func="cal.f.cvm.lm", scale_min.data=scale_min.data, a0=1, n0=length(scale_min.data))
      
      if(is.null(recal_mle_a$opt.a)){
        recal_mle_a2 <- NULL
        recal_mle_a2 <- opt_a0_addnoise(target_func="cal.f.cvm.lm", 
                                        scale_min.data=sample(scale_min.data, size=(length(scale_min.data)*0.8), replace = TRUE), 
                                        a0=1, n0=length(scale_min.data)) 
        if(is.null(recal_mle_a2$opt.a)){
          mle_a <- 1
          r_mle_a <- NA
        }else{
          mle_a <- recal_mle_a2$opt.a
          r_mle_a <- recal_mle_a2$r.opt.a
        }
      }else{
        mle_a <- recal_mle_a$opt.a
        r_mle_a <- recal_mle_a$r.opt.a
      }
    }else{
      mle_a <- cal_mle_a$opt.a
      r_mle_a <- cal_mle_a$r.opt.a
      mle_a_using_lm <- FALSE
    }
    # LME---
    if(mle_a_using_lm){
      lm_a <- mle_a
      r_lm_a <- r_mle_a
    }else{
      recal_lm_a <- NULL
      recal_lm_a <- opt_a0_addnoise(target_func="cal.f.cvm.lm", scale_min.data=scale_min.data, a0=1, n0=length(scale_min.data))
      
      if(is.null(recal_lm_a$opt.a)){
        recal_lm_a2 <- NULL
        recal_lm_a2 <- opt_a0_addnoise(target_func="cal.f.cvm.lm",
                                       scale_min.data=sample(scale_min.data, size=(length(scale_min.data)*0.8), replace = TRUE), 
                                       a0=1, n0=length(scale_min.data)) 
        if(is.null(recal_lm_a2$opt.a)){
          lm_a <- 1
          r_lm_a <- NA
        }else{
          lm_a <- recal_lm_a2$opt.a
          r_lm_a <- recal_lm_a2$r.opt.a
        }
      }else{
        lm_a <- recal_lm_a$opt.a
        r_lm_a <- recal_lm_a$r.opt.a
      }
    }
  }else if(db_method == "kfcv"){ # 3) k-folds cross-validation (kFCV) method
    # MLE---
    cal_mle_a <- NULL
    cal_mle_a <- opt_a0_resam_cv(target_func="cal.f.cvm.mle", sam.min.data=min.data, scale_min.data=scale_min.data, a0=1, n0=length(scale_min.data))
    
    if(is.null(cal_mle_a$opt.a)){
      mle_a_using_lm <- TRUE
      recal_mle_a <- NULL
      recal_mle_a <- opt_a0_resam_cv(target_func="cal.f.cvm.lm", sam.min.data=min.data, scale_min.data=scale_min.data, a0=1, n0=length(scale_min.data))
      
      if(is.null(recal_mle_a$opt.a)){
        recal_mle_a2 <- NULL
        recal_mle_a2 <- opt_a0_resam_cv(target_func="cal.f.cvm.lm", sam.min.data=min.data, a0=1, n0=length(scale_min.data), 
                                        scale_min.data=sample(scale_min.data, size=(length(scale_min.data)*0.8), replace = TRUE)) 
        if(is.null(recal_mle_a2$opt.a)){
          mle_a <- 1
          r_mle_a <- NA
        }else{
          mle_a <- recal_mle_a2$opt.a
          r_mle_a <- recal_mle_a2$r.opt.a
        }
      }else{
        mle_a <- recal_mle_a$opt.a
        r_mle_a <- recal_mle_a$r.opt.a
      }
    }else{
      mle_a <- cal_mle_a$opt.a
      r_mle_a <- cal_mle_a$r.opt.a
      mle_a_using_lm <- FALSE
    }
    # LME---
    if(mle_a_using_lm){
      lm_a <- mle_a
      r_lm_a <- r_mle_a
    }else{
      recal_lm_a <- NULL
      recal_lm_a <- opt_a0_resam_cv(target_func="cal.f.cvm.lm", sam.min.data=min.data, scale_min.data=scale_min.data, a0=1, n0=length(scale_min.data))
      
      if(is.null(recal_lm_a$opt.a)){
        recal_lm_a2 <- NULL
        recal_lm_a2 <- opt_a0_resam_cv(target_func="cal.f.cvm.lm", sam.min.data=min.data, a0=1, n0=length(scale_min.data),
                                       scale_min.data=sample(scale_min.data, size=(length(scale_min.data)*0.8), replace = TRUE)) 
        if(is.null(recal_lm_a2$opt.a)){
          lm_a <- 1
          r_lm_a <- NA
        }else{
          lm_a <- recal_lm_a2$opt.a
          r_lm_a <- recal_lm_a2$r.opt.a
        }
      }else{
        lm_a <- recal_lm_a$opt.a
        r_lm_a <- recal_lm_a$r.opt.a
      }
    }
  }else if(db_method == "boots"){ # 4) Bootstrap sampling of kFCV method
    len_B <- 50
    Boots_mle_a <- Boots_lm_a <- numeric(len_B)
    cal_mle_a <- cal_lm_a <- NULL
    for (len_b in 1:len_B) {
      # Resampling
      resam_min.data <- sample(min.data, size=length(min.data), replace = TRUE)
      # MLE
      cal_mle_a <- opt_a0_resam_cv(target_func="cal.f.cvm.mle", sam.min.data=resam_min.data, a0=1, 
                                   n0=length(scale_min.data), scale_min.data=scale_min.data)
      Boots_mle_a[len_b] <- cal_mle_a$opt.a
      # LME
      cal_lm_a <- opt_a0_resam_cv(target_func="cal.f.cvm.lm", sam.min.data=resam_min.data, a0=1, 
                                  n0=length(scale_min.data), scale_min.data=scale_min.data)
      Boots_lm_a[len_b] <- cal_lm_a$opt.a
    }
    mle_a <- mean(Boots_mle_a); r_mle_a <- Boots_mle_a
    lm_a <- mean(Boots_lm_a); r_lm_a <- Boots_lm_a
  }
  
  # 1) Maximum likelihood estimation--
  # Transformed data corresponding to est. hyperparameter
  data_m1 <- scale_min.data^(-mle_a)
  
  # Parameter estimation
  m1_fit <- m1_fit_se <- m1_fit_par <- m1_fit_nllh <- m1_fit_rl <- m1_fit_cvm <- m1_fit_ad <- NULL
  m1_fit <- gev.fit(data_m1, show=FALSE)
  m1_fit_cov <- m1_fit$cov
  if(any(is.na(m1_fit_cov)) || any(abs(m1_fit_cov) < 1e-04)){
    m1_fit_se <- boots.lme_par(data_m1, method="mle")
  }else{
    m1_fit_se <- m1_fit$se
  }
  m1_fit_par <- data.frame(est=m1_fit$mle, se=m1_fit_se)
  m1_fit_nllh<- m1_fit$nllh
  m1_fit$min.data <- min.data
  # Return level estimation
  m1_fit_rl <- getReturnLevels_npt(f=m1_fit, a=mle_a, rp=rp0, conf=conf0, sd_scale=sd_scale, trans=trans0)
  if(any(is.na(m1_fit_rl[,2]))){
    m1_fit_rl <- NULL
    m1_fit_rl <- boots.lme_rl_npt(data=data_m1, a=mle_a, min.data=min.data, rp=rp0, conf=conf0, sd_scale=sd_scale, method0="MLE", trans=trans0)
  }
  # GOF test
  m1_fit_cvm <- cvm.test(m1_fit$data, "pgev", m1_fit$mle[1], m1_fit$mle[2], m1_fit$mle[3])
  m1_fit_ad  <- ad.test(m1_fit$data, "pgev", m1_fit$mle[1], m1_fit$mle[2], m1_fit$mle[3])
  
  # 2) L-moments estimation--
  # Transformed data correspoding to est. hyperparameter
  data_m2 <- scale_min.data^(-lm_a)
  # Parameter estimation
  m2_fit <- m2_fit_se <- m2_fit_par <- m2_fit_rl <- m2_fit_cvm <- m2_fit_ad <- NULL
  m2_fit <- fevd(data_m2, method="Lmoments", type="GEV")
  m2_fit_se <- boots.lme_par(data_m2, method="l-moment")
  m2_fit_par<- data.frame(est=m2_fit$results, se=m2_fit_se)
  # Return level estimation
  m2_fit_rl <- boots.lme_rl_npt(data=data_m2, a=lm_a, min.data=min.data, rp=rp0, conf=conf0, sd_scale=sd_scale, trans=trans0)
  # GOF test
  m2_fit_cvm<- cvm.test(m2_fit$x, "pgev", m2_fit$results[1], m2_fit$results[2], m2_fit$results[3])
  m2_fit_ad <- ad.test(m2_fit$x, "pgev", m2_fit$results[1], m2_fit$results[2], m2_fit$results[3])
  
  # Combine all results----
  result <- list()
  result$mle_par <- m1_fit_par
  result$lm_par <- m2_fit_par
  result$nllh <- m1_fit_nllh
  result$mle_rl <- m1_fit_rl
  result$lm_rl <- m2_fit_rl
  result$mle_a <- mle_a
  result$lm_a <- lm_a
  result$r_mle_a <- r_mle_a
  result$r_lm_a <- r_lm_a
  
  result$gof_test <- data.frame(cvm.stat = c(m1_fit_cvm$statistic, m2_fit_cvm$statistic), 
                                cvm.pval = c(m1_fit_cvm$p.value, m2_fit_cvm$p.value),
                                ad.stat = c(m1_fit_ad$statistic, m2_fit_ad$statistic),
                                ad.pval = c(m1_fit_ad$p.value, m2_fit_ad$p.value))
  rownames(result$gof_test) <- c("MLE", "LME")
  
  return(result)
}

# Import minimum IAT data ----
id <- 90 # choose data station (90, 100, 115, 130, 146, )
data <- read.csv(paste0("/Users/thanawanp/workspace/5.InvGEV extension/Program/IAT data - 100mm/",id,".csv"))

# Extract min data
min.data <- data$min
summary(min.data) 

# Return periods and confidence interval setting
rp0 <- c(25,50,100)
conf0 <- 0.95

# Return level estimation using distribution-based (db) method ---
# Run resampling method
result_method1 <- est_funct(db_method = "resam", min.data=min.data, rp0=rp0, conf0=conf0)
result_method1$mle_par
result_method1$lm_par
result_method1$gof_test
result_method1$nllh
result_method1$mle_rl

result_method1$lm_rl
result_method1$mle_a
result_method1$lm_a
hist(result_method1$r_mle_a)
hist(result_method1$r_lm_a)

# other methods
result_method2 <- est_funct(db_method = "adnoise", min.data=min.data, rp0=rp0, conf0=conf0)  
result_method3 <- est_funct(db_method = "kfcv", min.data=min.data, rp0=rp0, conf0=conf0)
result_method4 <- est_funct(db_method = "boots", min.data=min.data, rp0=rp0, conf0=conf0)

