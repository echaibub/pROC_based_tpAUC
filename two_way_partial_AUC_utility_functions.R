
###########################
## load dependencies
###########################

library(pROC)
library(tpAUC)
library(ranger)

###########################
## utility functions
###########################

## proposed estimator function
##
pROCBasedTwoWayPartialAUC <- function(roc_object, 
                                      sensitivity_bound, 
                                      specificity_bound) {
  ## Inputs:
  ## roc_object: output of the "roc" function from the pROC package
  ## sensitivity_bound: the sensitivity threshold
  ## specificity_bound: the specificity threshold
  ##
  ## Outputs:
  ## tpAUC: estimated two-way partial AUC
  ## pAucSe: estimated partial AUC focusing on sensitivity
  ## pAucSp: estimated partial AUC focusing on specificity
  ## Auc: estimated (full) AUC
  
  ## Compute the full AUC using the "auc" functin from the pROC package.
  Auc <- pROC::auc(roc_object)[1]
  
  ## Compute the area that will be dropped.
  Arec <- sensitivity_bound * specificity_bound
  
  ## Compute the partial auc focusing on sensitivity.
  pAucSe <- pROC::auc(roc_object, 
                      partial.auc = c(1, sensitivity_bound), 
                      partial.auc.focus = "sensitivity", 
                      partial.auc.correct = FALSE)[1]
  
  ## Compute the partial auc focusing on specificity.
  pAucSp <- pROC::auc(roc_object, 
                      partial.auc = c(1, specificity_bound), 
                      partial.auc.focus = "specificity", 
                      partial.auc.correct = FALSE)[1]
  
  ## In situations where the ROC curve does not cross the
  ## are of interest (the rectangle with sides 1-sensitivity_bound 
  ## and 1-specificity_bound) we set the two-way partial AUC to 0.
  tpAUC <- 0
  
  ## Because the ROC is monotomic function a simple way 
  ## to check whether it crosses the area of interest is to simply
  ## check if the specificity value corresponding to sensitivity
  ## bound is larger than the specificity bound.
  ## The next line uses the "coords" function from the pROC package
  ## to compute the specificity value that corresponds to the 
  ## sensitivity bound.
  spec_at_sens_bound <- pROC::coords(roc_object, 
                                     x = sensitivity_bound, 
                                     input = "sensitivity", 
                                     ret = "specificity")
  
  ## Check the condition and, if it is true, compute the two-way
  ## partial AUC score.
  if (spec_at_sens_bound >= specificity_bound) {
    tpAUC <- pAucSe + pAucSp - (Auc - Arec)
  }
  
  ## Return the outputs.
  list(tpAUC = tpAUC,
       pAucSe = pAucSe, 
       pAucSp = pAucSp)
}


## alternative estimator function
##
pROCBasedTwoWayPartialAUC2 <- function(roc_object, 
                                       sensitivity_bound, 
                                       specificity_bound) {
  ## Inputs:
  ## roc_object: output of the "roc" function from the pROC package
  ## sensitivity_bound: the sensitivity threshold
  ## specificity_bound: the specificity threshold
  ##
  ## Outputs:
  ## tpAUC: estimated two-way partial AUC
  
  ## find the specificity value that corresponds to the
  ## sensitivity bound
  aux <- pROC::coords(roc_object, 
                      x = sensitivity_bound, 
                      input = "sensitivity", 
                      ret = "specificity")
  aux <- as.numeric(aux)
  
  tpauc <- 0
  if (aux >= specificity_bound) {
    
    ## compute the partial AUC
    pAucSp <- pROC::auc(roc_object, 
                        partial.auc = c(aux, specificity_bound), 
                        partial.auc.focus = "specificity", 
                        partial.auc.correct = FALSE)[1]
    
    ## compute are of the rectangle to be discarded
    A <- sensitivity_bound * (aux - specificity_bound)
    
    ## compute the two way partial AUC
    tpauc <- pAucSp - A
    
  }
  
  list(tpAUC = tpauc)
}


## function for simulating data from the binormal ROC model
##
SimulateDataBinormalROC <- function(n_N, n_P, mu_N, sig_N, mu_P, sig_P) {
  ## Inputs:
  ## n_N: number of negative cases
  ## n_P: number of positive cases
  ## mu_N, sig_N: mean and standard deviation of the negative cases
  ## mu_P, sig_P: mean and standard deviation of the positive cases
  ##
  ## Output:
  ## dat: data-frame containing labels and confidence values
  
  ## X: non-disease/negative
  ## Y: disease/positive
  X <- rnorm(n_N, mu_N, sig_N)
  Y <- rnorm(n_P, mu_P, sig_P)
  label <- c(rep(0, n_N), rep(1, n_P))
  conf <- c(X, Y)
  dat <- data.frame(label, conf)
  
  dat
}


## function for computing the true tpAUC value from
## the binormal ROC model
##
TrueBinormalTwoWayPartialAUC <- function(mu_N, 
                                         sig_N, 
                                         mu_P, 
                                         sig_P,
                                         sensitivity_bound,
                                         specificity_bound) {
  ## Inputs:
  ## n_N: number of negative cases
  ## n_P: number of positive cases
  ## mu_N, sig_N: mean and standard deviation of the negative cases
  ## mu_P, sig_P: mean and standard deviation of the positive cases
  ## sensitivity_bound: the sensitivity threshold
  ## specificity_bound: the specificity threshold
  ##
  ## Output:
  ## true_tpAUC: true two-way partial AUC value ("area A")
  ## true_AUC: true (full) AUC value
  ## area_2_drop: "area B" on the paper
  
  tpr_bound <- sensitivity_bound
  fpr_bound <- 1 - specificity_bound
  
  ## compute full AUC
  a <- (mu_P - mu_N)/sig_P
  b <- sig_N/sig_P
  true_AUC <- pnorm(a/sqrt(1 + b^2))
  
  ## find the bounds for the two-way partial AUC
  integrand <- function(z) {
    pnorm((mu_P - z)/sig_P) * dnorm((mu_N - z)/sig_N)/sig_N
  }
  lb <- mu_N - sig_N * qnorm(fpr_bound)
  ub <- mu_P - sig_P * qnorm(tpr_bound)
  
  ## fpr at the tpr_bound
  u_0 <- pnorm((mu_N - ub)/sig_N)
  
  ## compute the two way partial AUC
  true_tpAUC <- 0 
  area_2_drop <- 0
  if (fpr_bound > u_0) { ## check if the ROC crosses the area of interest
    true_pAUC <- integrate(integrand, lower = lb, upper = ub)$value
    area_2_drop <- tpr_bound * (fpr_bound - u_0)
    true_tpAUC <- true_pAUC - area_2_drop
  }
  
  list(true_tpAUC = true_tpAUC,
       true_AUC = true_AUC,
       area_2_drop = area_2_drop)
}


## function for running the first experiments comparing
## the proposed and original estimators
##
RunFirstComparisons <- function(n_P,
                                n_N,
                                n_runs) {
  ## Inputs:
  ## n_P: number of negative cases
  ## n_N: number of positive cases
  ## n_runs: number of replications
  ##
  ## Output:
  ## tpaucs: matrix of estimated tpAUC values for the proposed, alternative,
  ##         and original estimators
  ## computation_time_roc: matrix of computation times for the ROC curve 
  ## computation_time_pROC: matrix of computation times for the proposed 
  ##                        estimator
  ## computation_time_tpAUC: matrix of computation times for the original 
  ##                         estimator
  
  ## create objects for storing the results
  tpaucs <- matrix(NA, n_runs, 2)
  colnames(tpaucs) <- c("proposed", "original")
  computation_time_roc <- matrix(NA, n_runs, 3)
  colnames(computation_time_roc) <- c("user", "system", "elapsed")
  computation_time_pROC <- computation_time_roc
  computation_time_tpAUC <- computation_time_roc
  
  for (i in seq(n_runs)) {
    cat(i, "\n")
    
    ## randomly sample parameter values from the binormal ROC model
    mu_N <- runif(1, 0, 2)
    mu_P <- mu_N + runif(1, 0.1, 1.1)
    sig_N <- runif(1, 0.5, 1.5)
    sig_P <- runif(1, 0.5, 1.5)
    
    ## simulate the data
    dat <- SimulateDataBinormalROC(n_N, n_P, mu_N, sig_N, mu_P, sig_P)
    
    ## randomly sample the sensitivity and specificity thresholds
    sensitivity_bound <- runif(1, 0.2, 0.8)
    specificity_bound <- runif(1, 0.2, 0.8)
    
    ## fit the ROC curve and record the time it takes 
    t0 <- system.time(roc_fit <- pROC::roc(dat$label, 
                                           dat$conf, 
                                           auc = FALSE, 
                                           levels = c(0, 1), 
                                           direction = "<"))
    computation_time_roc[i,] <- as.numeric(t0)[1:3]
    
    ## compute the proposed estimator and record the time it takes 
    t1 <- system.time(o1 <- pROCBasedTwoWayPartialAUC(roc_fit, 
                                                      sensitivity_bound, 
                                                      specificity_bound))
    tpaucs[i, 1] <- o1$tpAUC
    ## add the time taken to compute the ROC and estimate the tpAUC
    computation_time_pROC[i,] <- as.numeric(t0)[1:3] + as.numeric(t1)[1:3]
    
    ## compute the original estimator and record the time it takes 
    t2 <- system.time(o2 <- tpAUC::tproc.est(response = as.factor(dat$label), 
                                             predictor = dat$conf, 
                                             threshold = c(1 - specificity_bound, sensitivity_bound)))
    tpaucs[i, 2] <- o2
    computation_time_tpAUC[i,] <- as.numeric(t2)[1:3]
  }
  
  list(tpaucs = tpaucs,
       computation_time_roc = computation_time_roc,
       computation_time_pROC = computation_time_pROC,
       computation_time_tpAUC = computation_time_tpAUC)
}



## function for computing the true tpAUC value from
## the binormal ROC model (it also returns a grid
## of FPR and TPR corresponding to the cutoff grid)
##
TrueBinormalTwoWayPartialAUC2 <- function(mu_N, 
                                         sig_N, 
                                         mu_P, 
                                         sig_P,
                                         sensitivity_bound,
                                         specificity_bound,
                                         cutoff_grid) {
  ## Inputs:
  ## n_N: number of negative cases
  ## n_P: number of positive cases
  ## mu_N, sig_N: mean and standard deviation of the negative cases
  ## mu_P, sig_P: mean and standard deviation of the positive cases
  ## sensitivity_bound: the sensitivity threshold
  ## specificity_bound: the specificity threshold
  ## cutoff_grid: grid of cutoff values (c)
  ##
  ## Output:
  ## true_tpAUC: true two-way partial AUC value ("area A")
  ## true_AUC: true (full) AUC value
  ## area_2_drop: "area B" on the paper
  ## TPR: TPR values corresponding to the cutoff_grid
  ## FPR: FPR values corresponding to the cutoff_grid
  
  ## compute specificities and sensitivities
  fpr <- pnorm((mu_N - cutoff_grid)/sig_N)
  tpr <- pnorm((mu_P - cutoff_grid)/sig_P)
  
  tpr_bound <- sensitivity_bound
  fpr_bound <- 1 - specificity_bound
  
  ## compute full AUC
  a <- (mu_P - mu_N)/sig_P
  b <- sig_N/sig_P
  true_AUC <- pnorm(a/sqrt(1 + b^2))
  
  ## find the bounds for the two-way partial AUC
  integrand <- function(z) {
    pnorm((mu_P - z)/sig_P) * dnorm((mu_N - z)/sig_N)/sig_N
  }
  lb <- mu_N - sig_N * qnorm(fpr_bound)
  ub <- mu_P - sig_P * qnorm(tpr_bound)
  
  ## fpr at the tpr_bound
  u_0 <- pnorm((mu_N - ub)/sig_N)
  
  ## compute the two way partial AUC
  true_tpAUC <- 0 
  area_2_drop <- 0
  if (fpr_bound > u_0) { ## check if the ROC crosses the area of interest
    true_pAUC <- integrate(integrand, lower = lb, upper = ub)$value
    area_2_drop <- tpr_bound * (fpr_bound - u_0)
    true_tpAUC <- true_pAUC - area_2_drop
  }
  
  list(true_tpAUC = true_tpAUC,
       true_AUC = true_AUC,
       area_2_drop = area_2_drop,
       TPR = tpr,
       FPR = fpr)
}


## function for estimating the bias of the tpAUC of the 
## original and proposed tpAUC estimators
##
tpAUCBias <- function(n_runs,
                      n_N,
                      n_P,
                      mu_N, 
                      sig_N, 
                      mu_P, 
                      sig_P,
                      sensitivity_bound,
                      specificity_bound) {
  ## Inputs:
  ## n_runs: number of simulated datasets used to estimate the bias
  ## n_N: number of negative cases
  ## n_P: number of positive cases
  ## mu_N, sig_N: mean and standard deviation of the negative cases
  ## mu_P, sig_P: mean and standard deviation of the positive cases
  ## sensitivity_bound: the sensitivity threshold
  ## specificity_bound: the specificity threshold
  ##
  ## Output:
  ## true_tpAUC: true two-way partial AUC value
  ## bias: vector with the estimated bias of the proposed and
  ##       original estimators
  ## tpAUC_hat: matrix storing the estimated tpAUC values for
  ##            the proposed and original estimators across all
  ##            n_runs simulations
  
  ## compute the true tpAUC value
  true_tpAUC <- TrueBinormalTwoWayPartialAUC(mu_N, 
                                             sig_N, 
                                             mu_P, 
                                             sig_P,
                                             sensitivity_bound,
                                             specificity_bound)$true_tpAUC
  tpAUC_hat <- matrix(NA, n_runs, 2)
  colnames(tpAUC_hat) <- c("pROC", "tpAUC")
  for (i in seq(n_runs)) {
    
    ## simulate data
    dat <- SimulateDataBinormalROC(n_N, n_P, mu_N, sig_N, mu_P, sig_P)
    
    ## estimate the ROC curve
    roc_fit <- roc(dat$label, dat$conf, auc = FALSE, 
                   levels = c(0, 1), direction = "<")
    
    ## compute tpAUC using the proposed estimator
    tpAUC_hat[i, 1] <- pROCBasedTwoWayPartialAUC(roc_fit, 
                                        sensitivity_bound, 
                                        specificity_bound)$tpAUC
    
    ## compute tpAUC using the original estimator
    tpAUC_hat[i, 2] <- tpAUC::tproc.est(response = as.factor(dat$label), 
                                 predictor = dat$conf, 
                                 threshold = c(1 - specificity_bound, sensitivity_bound))
    
  }
  
  ## compute the estimator's biases
  bias <- apply(tpAUC_hat, 2, mean) - true_tpAUC
  
  list(true_tpAUC = true_tpAUC,
       bias = bias,
       tpAUC_hat = tpAUC_hat)
}


## function for running the bias comparisons experiments
##
RunBiasComparisons <- function(n_sim,
                               n_runs,
                               n_N,
                               n_P,
                               mu_N_range, 
                               sig_N_range, 
                               mu_P_range, ## needs to be above mu_N 
                               sig_P_range,
                               sensitivity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               specificity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               my_seeds) {
  ## Inputs:
  ## n_sim: number of replications of the experiment
  ## n_runs: number of simulated datasets used to estimate the bias
  ## n_N: number of negative cases
  ## n_P: number of positive cases
  ## mu_N, sig_N: mean and standard deviation of the negative cases
  ## mu_P, sig_P: mean and standard deviation of the positive cases
  ## sensitivity_bound_grid: the grid of sensitivity threshold values
  ## specificity_bound_grid: the grid of specificity threshold values
  ## my_seeds: vector of seeds for random number generator
  ##
  ## Output:
  ## true_tpAUC: true two-way partial AUC value
  ## pROC_bias: matrix with the estimated bias of the proposed estimator
  ##            across all n_sim replications and combinations of values
  ##            of the sensitivity_bound_grid and specificity_bound_grid vectors
  ## tpAUC_bias: same as above for the original estimator
  ## estimates_pROC: array storing the estimated tpAUC values for
  ##                 the proposed estimator across all n_sim replications,
  ##                 all n_runs simulations, and all sensitivity/specificity
  ##                 combinations
  ## estimates_tpAUC: same as above for the original estimator
  
  ## get all possible combinations of the values in the 
  ## sensitivity_bound_grid and specificity_bound_grid vectors
  n_sensitivity_bound <- length(sensitivity_bound_grid)
  n_specificity_bound <- length(specificity_bound_grid)
  n_bound <- n_sensitivity_bound*n_specificity_bound
  sens_spec_combinations <- expand.grid(sensitivity_bound_grid, 
                                        specificity_bound_grid)
  sens_spec_combinations <- paste0("(", 
                                   sens_spec_combinations[, 2], 
                                   " , ", 
                                   sens_spec_combinations[, 1], 
                                   ")")
  
  ## create objects to store the results
  true_tpAUC <- matrix(NA, n_sim, n_bound)
  colnames(true_tpAUC) <- sens_spec_combinations
  pROC_bias <- true_tpAUC
  tpAUC_bias <- true_tpAUC
  estimates_pROC <- array(dim = c(n_sim, n_bound, n_runs), 
                          dimnames = list(NULL, sens_spec_combinations, NULL))
  estimates_tpAUC <- array(dim = c(n_sim, n_bound, n_runs), 
                           dimnames = list(NULL, sens_spec_combinations, NULL))
  
  ## run the replications
  for (k in seq(n_sim)) {
    cat(k, "\n")
    
    ## for each replication, randomly draw a distinct set of parameter values 
    ## from the binormal ROC model 
    set.seed(my_seeds[k])
    mu_N <- runif(1, mu_N_range[1], mu_N_range[2])
    sig_N <- runif(1, sig_N_range[1], sig_N_range[2])
    mu_P <- runif(1, mu_P_range[1], mu_P_range[2])
    sig_P <- runif(1, sig_P_range[1], sig_P_range[2])   
    idx <- 1
    
    ## estimate the bias for each sensitivity/specificity combination
    for (i in seq(n_sensitivity_bound)) {
      for (j in seq(n_specificity_bound)) {
        set.seed(my_seeds[k] + idx)
        aux <- tpAUCBias(n_runs,
                         n_N,
                         n_P,
                         mu_N, 
                         sig_N, 
                         mu_P, 
                         sig_P,
                         sensitivity_bound_grid[i],
                         specificity_bound_grid[j])
        true_tpAUC[k, idx] <- aux$true_tpAUC
        pROC_bias[k, idx] <- as.numeric(aux$bias["pROC"])
        tpAUC_bias[k, idx] <- as.numeric(aux$bias["tpAUC"])
        estimates_pROC[k, idx, ] <- aux$tpAUC_hat[, "pROC"]
        estimates_tpAUC[k, idx, ] <- aux$tpAUC_hat[, "tpAUC"]
        idx <- idx + 1
      }
    }
  }
  
  list(true_tpAUC = true_tpAUC,
       pROC_bias = pROC_bias,
       tpAUC_bias = tpAUC_bias,
       estimates_pROC = estimates_pROC,
       estimates_tpAUC = estimates_tpAUC)
}


## function for running the additional computation time comparison 
## experiments
##
ComputationTimeSimulations <- function(n_runs,
                                       n_N,
                                       n_P,
                                       mu_N,
                                       sig_N,
                                       mu_P,
                                       sig_P,
                                       sensitivity_bound,
                                       specificity_bound) {
  ## Inputs:
  ## n_runs: number of simulated datasets used to estimate the bias
  ## n_N: number of negative cases
  ## n_P: number of positive cases
  ## mu_N, sig_N: mean and standard deviation of the negative cases
  ## mu_P, sig_P: mean and standard deviation of the positive cases
  ## sensitivity_bound: the sensitivity threshold
  ## specificity_bound: the specificity threshold
  ##
  ## Output:
  ## tpaucs: matrix of estimated tpAUC values for the proposed, alternative,
  ##         and original estimators
  ## computation_time_roc: matrix of computation times for the ROC curve 
  ## computation_time_pROC_1: matrix of computation times for the proposed 
  ##                          estimator
  ## computation_time_pROC_2: matrix of computation times for the alternative 
  ##                          estimator
  ## computation_time_tpAUC: matrix of computation times for the original 
  ##                         estimator
   
  ## create objects to store the results
  tpaucs <- matrix(NA, n_runs, 3)
  ## pROC_1: proposed estimator
  ## pROC_2: alternative estimator
  ## tpAUC: original estimator
  colnames(tpaucs) <- c("pROC_1", "pROC_2", "tpAUC")
  computation_time_roc <- matrix(NA, n_runs, 3)
  colnames(computation_time_roc) <- c("user", "system", "elapsed")
  computation_time_pROC_1 <- computation_time_roc
  computation_time_pROC_2 <- computation_time_pROC_1
  computation_time_tpAUC <- computation_time_pROC_1
  
  for (i in seq(n_runs)) {
    cat(i, "\n")
    
    ## simulate the data
    dat <- SimulateDataBinormalROC(n_N, n_P, mu_N, sig_N, mu_P, sig_P)
    
    ## fit the ROC curve and record the time it takes 
    t0 <- system.time(roc_fit <- pROC::roc(dat$label, 
                                           dat$conf, 
                                           auc = FALSE, 
                                           levels = c(0, 1), 
                                           direction = "<"))
    computation_time_roc[i,] <- as.numeric(t0)[1:3]
    
    ## compute the proposed estimator and record the time it takes 
    t1 <- system.time(o1 <- pROCBasedTwoWayPartialAUC(roc_fit, 
                                                      sensitivity_bound, 
                                                      specificity_bound))
    tpaucs[i, 1] <- o1$tpAUC
    ## add the time taken to compute the ROC and estimate the tpAUC
    computation_time_pROC_1[i,] <- as.numeric(t0)[1:3] + as.numeric(t1)[1:3]
    
    ## compute the alternative estimator and record the time it takes
    t2 <- system.time(o2 <- pROCBasedTwoWayPartialAUC2(roc_fit, 
                                                       sensitivity_bound, 
                                                       specificity_bound))
    tpaucs[i, 2] <- o2$tpAUC
    ## add the time taken to compute the ROC and estimate the tpAUC
    computation_time_pROC_2[i,] <- as.numeric(t0)[1:3] + as.numeric(t2)[1:3]
    
    ## compute the original estimator and record the time it takes
    t3 <- system.time(o3 <- tpAUC::tproc.est(response = as.factor(dat$label), 
                                      predictor = dat$conf, 
                                      threshold = c(1 - specificity_bound, sensitivity_bound)))
    tpaucs[i, 3] <- o3
    computation_time_tpAUC[i,] <- as.numeric(t3)[1:3]
  }
  
  list(tpaucs = tpaucs,
       computation_time_roc = computation_time_roc,
       computation_time_pROC_1 = computation_time_pROC_1,
       computation_time_pROC_2 = computation_time_pROC_2,
       computation_time_tpAUC = computation_time_tpAUC)
}


## function for organizing the output of the first additional
## computation time experiment from the proposed estimator
##
ShapeUserTimepROC_1 <- function(x, n_grid) {
  ## Inputs
  ## x: list storing the outputs of all experiments
  ## n_grid: grid of sample size values
  ##
  ## Output
  ## out: matrix containing the computing times of the 
  ##      proposed estimator across all replications of
  ##      each experiment
  
  k <- length(n_grid)
  out <- matrix(NA, nrow(x[[1]]$computation_time_roc), k)
  colnames(out) <- n_grid
  for (i in seq(k)) {
    out[, i] <- x[[i]]$computation_time_roc[, "user"] + x[[i]]$computation_time_pROC_1[, "user"]
  }
  
  out
}


## function for organizing the output of the first additional
## computation time experiment from the alternative estimator
##
ShapeUserTimepROC_2 <- function(x, n_grid) {
  ## Inputs
  ## x: list storing the outputs of all experiments
  ## n_grid: grid of sample size values
  ##
  ## Output
  ## out: matrix containing the computing times of the 
  ##      alternative estimator across all replications of
  ##      each experiment
  k <- length(n_grid)
  out <- matrix(NA, nrow(x[[1]]$computation_time_roc), k)
  colnames(out) <- n_grid
  for (i in seq(k)) {
    out[, i] <- x[[i]]$computation_time_roc[, "user"] + x[[i]]$computation_time_pROC_2[, "user"]
  }
  
  out
}


## function for organizing the output of the first additional
## computation time experiment from the original estimator
##
ShapeUserTimetpAUC <- function(x, n_grid) {
  ## Inputs
  ## x: list storing the outputs of all experiments
  ## n_grid: grid of sample size values
  ##
  ## Output
  ## out: matrix containing the computing times of the 
  ##      original estimator across all replications of
  ##      each experiment
  
  k <- length(n_grid)
  out <- matrix(NA, nrow(x[[1]]$computation_time_roc), k)
  colnames(out) <- n_grid
  for (i in seq(k)) {
    out[, i] <- x[[i]]$computation_time_tpAUC[, "user"]
  }
  
  out
}


##############################################
## bootstrap confidence interval functions
##############################################

## function for computing the bootstrap assisted 
## confidence interval based on the proposed estimator
##
pROCBasedDeltaCITime <- function(n_boot,
                                 labels,
                                 confidence_1,
                                 confidence_2,
                                 sensitivity_bound, 
                                 specificity_bound,
                                 alpha = 0.05) {
  ## Inputs:
  ## n_boot: number of bootstraps
  ## labels: the label data
  ## confidence_1: confidence scores from the first classifier
  ## confidence_2: confidence scores from the second classifier
  ## sensitivity_bound: the sensitivity threshold
  ## specificity_bound: the specificity threshold
  ## alpha: significance level for the (1 - alpha) confidence interval
  ##
  ## Output:
  ## total_computing_time: the total computing time for calculating the
  ##                       observed statistics, all bootstraps, and the
  ##                       confidence interval
  ## computing_time: matrix with n_boot + 2 rows storing the computing 
  ##                 time for calculating the observed statistics, all 
  ##                 bootstraps, and the confidence interval
  ## ci_stats: confidence interval statistics (confidence interval, 
  ##           observed Delta tpAUC statistic, bootstrap variance, 
  ##           bootstrap Delta tpAUC statistics)
  ## obs_stats: observed statistics (difference of the tpAUCs of 
  ##            classifiers 1 and 2, estimated tpAUC for classifier 1,
  ##            estimated tpAUC for classifier 2)
  
  ## function for computing the difference of tpAUCs statistic
  DeltaStatistic <- function(labels,
                             confidence_1,
                             confidence_2,
                             sensitivity_bound, 
                             specificity_bound) {
    ## Inputs:
    ## labels: the label data
    ## confidence_1: confidence scores from the first classifier
    ## confidence_2: confidence scores from the second classifier
    ## sensitivity_bound: the sensitivity threshold
    ## specificity_bound: the specificity threshold
    ##
    ## Output:
    ## Delta: difference of the tpAUCs of classifiers 1 and 2
    ## tpAUC_1: estimated tpAUC for classifier 1
    ## tpAUC_2: estimated tpAUC for classifier 2
    
    ## compute ROC curves for classifiers 1 and 2
    roc_1 <- pROC::roc(labels, confidence_1, auc = FALSE, 
                       levels = c(0, 1), direction = "<")
    roc_2 <- pROC::roc(labels, confidence_2, auc = FALSE, 
                       levels = c(0, 1), direction = "<")
    
    ## compute the tpAUC values for classifiers 1 and 2 and the Delta tpAUC
    tpAUC_1 <- pROCBasedTwoWayPartialAUC(roc_1, 
                                         sensitivity_bound, 
                                         specificity_bound)$tpAUC
    tpAUC_2 <- pROCBasedTwoWayPartialAUC(roc_2, 
                                         sensitivity_bound, 
                                         specificity_bound)$tpAUC
    Delta <- tpAUC_1 - tpAUC_2
    
    list(Delta = Delta,
         tpAUC_1 = tpAUC_1,
         tpAUC_2 = tpAUC_2)
  }
  
  ## function for computing the Delta tpAUC statistic on 
  ## one bootstrap sample
  BootDeltaStatistic <- function(n,
                                 labels,
                                 confidence_1,
                                 confidence_2,
                                 sensitivity_bound, 
                                 specificity_bound) {
    ## Inputs:
    ## n: sample size
    ## labels: the label data
    ## confidence_1: confidence scores from the first classifier
    ## confidence_2: confidence scores from the second classifier
    ## sensitivity_bound: the sensitivity threshold
    ## specificity_bound: the specificity threshold
    ##
    ## Output:
    ## difference of the tpAUCs of classifiers 1 and 2 on the
    ## bootstrap sample
    
    ## resample the data with replacement
    idx <- sample(seq(n), n, replace = TRUE)
    boot_labels <- labels[idx]
    boot_conf_1 <- confidence_1[idx]
    boot_conf_2 <- confidence_2[idx]
    
    ## compute delta statistic on the bootstrap sample
    boot_roc_1 <- pROC::roc(boot_labels, boot_conf_1, auc = FALSE, 
                            levels = c(0, 1), direction = "<")
    boot_roc_2 <- pROC::roc(boot_labels, boot_conf_2, auc = FALSE, 
                            levels = c(0, 1), direction = "<")
    boot_tpAUC_1 <- pROCBasedTwoWayPartialAUC(boot_roc_1, 
                                              sensitivity_bound, 
                                              specificity_bound)$tpAUC
    boot_tpAUC_2 <- pROCBasedTwoWayPartialAUC(boot_roc_2, 
                                              sensitivity_bound, 
                                              specificity_bound)$tpAUC
    
    boot_tpAUC_1 - boot_tpAUC_2
  }
  
  ## function for computing the asymptotic confidence interval
  AsymptoticCI <- function(n,
                           obs_stat,
                           boot_stat,
                           alpha) {
    ## Inputs:
    ## n: sample size
    ## obs_stat: observed Delta tpAUC statistic
    ## boot_stat: vector with the bootstrap Delta tpAUC statistics
    ## alpha: significance level for the (1 - alpha) confidence interval
    ##
    ## Output:
    ## CI: confidence interval for the Delta tpAUC statistic
    ## obs_stat: observed Delta tpAUC statistic
    ## boot_var: bootstrap variance, 
    ##           bootstrap Delta tpAUC statistics)
    ## boot_stat: vector with the bootstrap Delta tpAUC statistics
    
    ## compute the bootstrap variance estimate of the Delta statistic
    var_boot <- mean((boot_stat - mean(boot_stat))^2)
    
    ## compute the bootstrap-assisted asymptotic confidence interval
    Z <- qnorm(1 - alpha/2, lower.tail = TRUE)
    CI <- c(obs_stat - Z * sqrt(var_boot/n),
            obs_stat + Z * sqrt(var_boot/n))
    names(CI) <- c("lower_bound", "upper_bound")
    
    list(CI = CI,
         obs_stat = obs_stat,
         var_boot = var_boot,
         boot_stat = boot_stat)
  }
  
  computing_time <- matrix(NA, n_boot + 2, 3)
  colnames(computing_time) <- c("user", "system", "elapsed")
  
  cat("observed data ", "\n")
  obs_time <- system.time(obs_stats <- DeltaStatistic(labels,
                                                     confidence_1,
                                                     confidence_2,
                                                     sensitivity_bound, 
                                                     specificity_bound))
  computing_time[1,] <- obs_time[1:3]
  
  ## run bootstraps
  n <- length(labels)
  boot_stat <- rep(NA, n_boot)
  for (i in seq(n_boot)) {
    cat("bootstrap ", i, "\n")
    
    ## compute delta statistic on the bootstrap sample
    boot_time <- system.time(boot_stat[i] <- BootDeltaStatistic(n,
                                                                labels,
                                                                confidence_1,
                                                                confidence_2,
                                                                sensitivity_bound, 
                                                                specificity_bound))
    computing_time[i + 1,] <- boot_time[1:3]
  }
  
  ## compute the bootstrap-assisted asymptotic confidence interval
  ci_time <- system.time(ci_stats <- AsymptoticCI(n,
                                                  obs_stats$Delta,
                                                  boot_stat,
                                                  alpha))
  computing_time[i + 2,] <- ci_time[1:3]
  total_computing_time <- apply(computing_time, 2, sum)
  
  list(total_computing_time = total_computing_time,
       computing_time = computing_time,
       ci_stats = ci_stats,
       obs_stats = obs_stats)
}


## function for computing the bootstrap assisted 
## confidence interval based on the proposed estimator
## (contrary to the pROCBasedDeltaCITime() function, this 
## function does not record the time taken for calculating the
## observed statistics, all bootstraps, and the confidence interval)
##
pROCBasedDeltaCI <- function(n_boot,
                             labels,
                             confidence_1,
                             confidence_2,
                             sensitivity_bound, 
                             specificity_bound,
                             alpha = 0.05) {
  ## Inputs:
  ## n_boot: number of bootstraps
  ## labels: the label data
  ## confidence_1: confidence scores from the first classifier
  ## confidence_2: confidence scores from the second classifier
  ## sensitivity_bound: the sensitivity threshold
  ## specificity_bound: the specificity threshold
  ## alpha: significance level for the (1 - alpha) confidence interval
  ##
  ## Output:
  ## CI: confidence interval
  ## obs_stat: observed Delta tpAUC statistic
  ## var_boot: bootstrap variance, 
  ## boot_stat: vector with the bootstrap Delta tpAUC statistics
  
  ## compute ROC curves for classifiers 1 and 2
  obs_roc_1 <- pROC::roc(labels, confidence_1, auc = FALSE, 
                         levels = c(0, 1), direction = "<")
  obs_roc_2 <- pROC::roc(labels, confidence_2, auc = FALSE, 
                         levels = c(0, 1), direction = "<")
  
  ## compute the tpAUC values for classifiers 1 and 2 and the Delta tpAUC
  obs_tpAUC_1 <- pROCBasedTwoWayPartialAUC(obs_roc_1, 
                                           sensitivity_bound, 
                                           specificity_bound)$tpAUC
  obs_tpAUC_2 <- pROCBasedTwoWayPartialAUC(obs_roc_2, 
                                           sensitivity_bound, 
                                           specificity_bound)$tpAUC
  obs_stat <- obs_tpAUC_1 - obs_tpAUC_2
  
  ## run bootstraps
  n <- length(labels)
  boot_stat <- rep(NA, n_boot)
  for (i in seq(n_boot)) {
    
    ## resample the data with replacement
    idx <- sample(seq(n), n, replace = TRUE)
    boot_labels <- labels[idx]
    boot_conf_1 <- confidence_1[idx]
    boot_conf_2 <- confidence_2[idx]
    
    ## compute delta statistic on the bootstrap sample
    boot_roc_1 <- pROC::roc(boot_labels, boot_conf_1, auc = FALSE, 
                            levels = c(0, 1), direction = "<")
    boot_roc_2 <- pROC::roc(boot_labels, boot_conf_2, auc = FALSE, 
                            levels = c(0, 1), direction = "<")
    boot_tpAUC_1 <- pROCBasedTwoWayPartialAUC(boot_roc_1, 
                                              sensitivity_bound, 
                                              specificity_bound)$tpAUC
    boot_tpAUC_2 <- pROCBasedTwoWayPartialAUC(boot_roc_2, 
                                              sensitivity_bound, 
                                              specificity_bound)$tpAUC
    boot_stat[i] <- boot_tpAUC_1 - boot_tpAUC_2
  }
  
  ## compute the bootstrap variance estimate of the Delta statistic
  var_boot <- mean((boot_stat - mean(boot_stat))^2)
  
  ## compute the bootstrap-assisted asymptotic confidence interval
  Z <- qnorm(alpha/2, lower.tail = FALSE)
  CI <- c(obs_stat - Z * sqrt(var_boot/n),
          obs_stat + Z * sqrt(var_boot/n))
  names(CI) <- c("lower_bound", "upper_bound")
  
  list(CI = CI,
       obs_stat = obs_stat,
       var_boot = var_boot,
       boot_stat = boot_stat)
}



## function for computing the bootstrap assisted 
## confidence interval based on the original estimator
##
tpAUCBasedDeltaCITime <- function(n_boot,
                              labels,
                              confidence_1,
                              confidence_2,
                              sensitivity_bound, 
                              specificity_bound,
                              alpha = 0.05) {
  ## Inputs:
  ## n_boot: number of bootstraps
  ## labels: the label data
  ## confidence_1: confidence scores from the first classifier
  ## confidence_2: confidence scores from the second classifier
  ## sensitivity_bound: the sensitivity threshold
  ## specificity_bound: the specificity threshold
  ## alpha: significance level for the (1 - alpha) confidence interval
  ##
  ## Output:
  ## total_computing_time: the total computing time for calculating the
  ##                       observed statistics, all bootstraps, and the
  ##                       confidence interval
  ## computing_time: matrix with n_boot + 2 rows storing the computing 
  ##                 time for calculating the observed statistics, all 
  ##                 bootstraps, and the confidence interval
  ## ci_stats: confidence interval statistics (confidence interval, 
  ##           observed Delta tpAUC statistic, bootstrap variance, 
  ##           bootstrap Delta tpAUC statistics)
  ## obs_stats: observed statistics (difference of the tpAUCs of 
  ##            classifiers 1 and 2, estimated tpAUC for classifier 1,
  ##            estimated tpAUC for classifier 2)
  
  ## function for computing the difference of tpAUCs statistic
  DeltaStatistic <- function(labels,
                             confidence_1,
                             confidence_2,
                             thr) {
    ## Inputs:
    ## labels: the label data
    ## confidence_1: confidence scores from the first classifier
    ## confidence_2: confidence scores from the second classifier
    ## thr: vector with FPR and TPR thresholds
    ##
    ## Output:
    ## Delta: difference of the tpAUCs of classifiers 1 and 2
    ## tpAUC_1: estimated tpAUC for classifier 1
    ## tpAUC_2: estimated tpAUC for classifier 2
    
    ## compute the tpAUC values for classifiers 1 and 2 and the Delta tpAUC
    tpAUC_1 <- tpAUC::tproc.est(response = as.factor(labels), 
                                    predictor = confidence_1, 
                                    threshold = thr)
    tpAUC_2 <- tpAUC::tproc.est(response = as.factor(labels), 
                                    predictor = confidence_2, 
                                    threshold = thr)
    Delta <- tpAUC_1 - tpAUC_2
    list(Delta = Delta,
         tpAUC_1 = tpAUC_1,
         tpAUC_2 = tpAUC_2)
  }
  
  ## function for computing the Delta tpAUC statistic on 
  ## one bootstrap sample
  BootDeltaStatistic <- function(n,
                                 labels,
                                 confidence_1,
                                 confidence_2,
                                 thr) {
    ## Inputs:
    ## n: sample size
    ## labels: the label data
    ## confidence_1: confidence scores from the first classifier
    ## confidence_2: confidence scores from the second classifier
    ## thr: vector with FPR and TPR thresholds
    ##
    ## Output:
    ## difference of the tpAUCs of classifiers 1 and 2 on the
    ## bootstrap sample
    
    ## resample the data with replacement
    idx <- sample(seq(n), n, replace = TRUE)
    boot_labels <- labels[idx]
    boot_conf_1 <- confidence_1[idx]
    boot_conf_2 <- confidence_2[idx]
    
    ## compute delta statistic on the bootstrap sample
    boot_tpAUC_1 <- tpAUC::tproc.est(response = as.factor(boot_labels), 
                                     predictor = boot_conf_1, 
                                     threshold = thr)
    boot_tpAUC_2 <- tpAUC::tproc.est(response = as.factor(boot_labels), 
                                     predictor = boot_conf_2, 
                                     threshold = thr)
    
    boot_tpAUC_1 - boot_tpAUC_2
  }
  
  ## function for computing the asymptotic confidence interval
  AsymptoticCI <- function(n,
                           obs_stat,
                           boot_stat,
                           alpha) {
    ## Inputs:
    ## n: sample size
    ## obs_stat: observed Delta tpAUC statistic
    ## boot_stat: vector with the bootstrap Delta tpAUC statistics
    ## alpha: significance level for the (1 - alpha) confidence interval
    ##
    ## Output:
    ## CI: confidence interval for the Delta tpAUC statistic
    ## obs_stat: observed Delta tpAUC statistic
    ## boot_var: bootstrap variance, 
    ##           bootstrap Delta tpAUC statistics)
    ## boot_stat: vector with the bootstrap Delta tpAUC statistics
    
    ## compute the bootstrap variance estimate of the Delta statistic
    var_boot <- mean((boot_stat - mean(boot_stat))^2)
    
    ## compute the bootstrap-assisted asymptotic confidence interval
    Z <- qnorm(1 - alpha/2, lower.tail = TRUE)
    CI <- c(obs_stat - Z * sqrt(var_boot/n),
            obs_stat + Z * sqrt(var_boot/n))
    names(CI) <- c("lower_bound", "upper_bound")
    
    list(CI = CI,
         obs_stat = obs_stat,
         var_boot = var_boot,
         boot_stat = boot_stat)
  }
  
  thr <- c(1 - specificity_bound, sensitivity_bound)
  computing_time <- matrix(NA, n_boot + 2, 3)
  colnames(computing_time) <- c("user", "system", "elapsed")
  
  cat("observed data ", "\n")
  obs_time <- system.time(obs_stats <- DeltaStatistic(labels,
                                                     confidence_1,
                                                     confidence_2,
                                                     thr))
  computing_time[1,] <- obs_time[1:3]
  
  ## run bootstraps
  n <- length(labels)
  boot_stat <- rep(NA, n_boot)
  for (i in seq(n_boot)) {
    cat("bootstrap ", i, "\n")
    
    ## compute delta statistic on the bootstrap sample
    boot_time <- system.time(boot_stat[i] <- BootDeltaStatistic(n,
                                                                labels,
                                                                confidence_1,
                                                                confidence_2,
                                                                thr))
    computing_time[i + 1,] <- boot_time[1:3]
  }
  
  ## compute the bootstrap-assisted asymptotic confidence interval
  ci_time <- system.time(ci_stats <- AsymptoticCI(n,
                                                  obs_stats$Delta,
                                                  boot_stat,
                                                  alpha))
  computing_time[i + 2,] <- ci_time[1:3]
  total_computing_time <- apply(computing_time, 2, sum)
  
  list(total_computing_time = total_computing_time,
       computing_time = computing_time,
       ci_stats = ci_stats,
       obs_stats = obs_stats)
}


## same as the tpAUCBasedDeltaCITime function above except
## that it saves the output after each (bootstrap) computation
##
tpAUCBasedDeltaCITimeSave <- function(n_boot,
                                  labels,
                                  confidence_1,
                                  confidence_2,
                                  sensitivity_bound, 
                                  specificity_bound,
                                  alpha = 0.05,
                                  file_name) {
  DeltaStatistic <- function(labels,
                             confidence_1,
                             confidence_2,
                             thr) {
    ## compute the tpAUC values for classifiers 1 and 2 and the Delta tpAUC
    tpAUC_1 <- tpAUC::tproc.est(response = as.factor(labels), 
                                predictor = confidence_1, 
                                threshold = thr)
    tpAUC_2 <- tpAUC::tproc.est(response = as.factor(labels), 
                                predictor = confidence_2, 
                                threshold = thr)
    Delta <- tpAUC_1 - tpAUC_2
    list(Delta = Delta,
         tpAUC_1 = tpAUC_1,
         tpAUC_2 = tpAUC_2)
  }
  BootDeltaStatistic <- function(n,
                                 labels,
                                 confidence_1,
                                 confidence_2,
                                 thr) {
    ## resample the data with replacement
    idx <- sample(seq(n), n, replace = TRUE)
    boot_labels <- labels[idx]
    boot_conf_1 <- confidence_1[idx]
    boot_conf_2 <- confidence_2[idx]
    
    ## compute delta statistic on the bootstrap sample
    boot_tpAUC_1 <- tpAUC::tproc.est(response = as.factor(boot_labels), 
                                     predictor = boot_conf_1, 
                                     threshold = thr)
    boot_tpAUC_2 <- tpAUC::tproc.est(response = as.factor(boot_labels), 
                                     predictor = boot_conf_2, 
                                     threshold = thr)
    boot_tpAUC_1 - boot_tpAUC_2
  }
  AsymptoticCI <- function(n,
                           obs_stat,
                           boot_stat,
                           alpha) {
    ## compute the bootstrap variance estimate of the Delta statistic
    var_boot <- mean((boot_stat - mean(boot_stat))^2)
    
    ## compute the bootstrap-assisted asymptotic confidence interval
    Z <- qnorm(1 - alpha/2, lower.tail = TRUE)
    CI <- c(obs_stat - Z * sqrt(var_boot/n),
            obs_stat + Z * sqrt(var_boot/n))
    names(CI) <- c("lower_bound", "upper_bound")
    
    list(CI = CI,
         obs_stat = obs_stat,
         var_boot = var_boot,
         boot_stat = boot_stat)
  }
  
  thr <- c(1 - specificity_bound, sensitivity_bound)
  computing_time <- matrix(NA, n_boot + 2, 3)
  colnames(computing_time) <- c("user", "system", "elapsed")
  
  cat("observed data ", "\n")
  obs_time <- system.time(obs_stats <- DeltaStatistic(labels,
                                                      confidence_1,
                                                      confidence_2,
                                                      thr))
  computing_time[1,] <- obs_time[1:3]
  
  save(obs_stats, file = file_name, compress = TRUE)
  
  ## run bootstraps
  n <- length(labels)
  boot_stat <- rep(NA, n_boot)
  for (i in seq(n_boot)) {
    cat("bootstrap ", i, "\n")
    
    ## compute delta statistic on the bootstrap sample
    boot_time <- system.time(boot_stat[i] <- BootDeltaStatistic(n,
                                                                labels,
                                                                confidence_1,
                                                                confidence_2,
                                                                thr))
    computing_time[i + 1,] <- boot_time[1:3]
    
    save(obs_stats, boot_stat, computing_time, file = file_name, compress = TRUE)
  }
  
  ## compute the bootstrap-assisted asymptotic confidence interval
  ci_time <- system.time(ci_stats <- AsymptoticCI(n,
                                                  obs_stats$Delta,
                                                  boot_stat,
                                                  alpha))
  computing_time[i + 2,] <- ci_time[1:3]
  total_computing_time <- apply(computing_time, 2, sum)
  
  save(obs_stats, 
       boot_stat, 
       ci_stats, 
       computing_time, 
       total_computing_time, 
       file = file_name, compress = TRUE)
  
  list(total_computing_time = total_computing_time,
       computing_time = computing_time,
       ci_stats = ci_stats,
       obs_stats = obs_stats)
}


## function for computing the bootstrap assisted 
## confidence interval based on the original estimator
## (contrary to the tpAUCBasedDeltaCITime() function, this 
## function does not record the time taken for calculating the
## observed statistics, all bootstraps, and the confidence interval)
##
tpAUCBasedDeltaCI <- function(n_boot,
                              labels,
                              confidence_1,
                              confidence_2,
                              sensitivity_bound, 
                              specificity_bound,
                              alpha = 0.05) {
  ## Inputs:
  ## n_boot: number of bootstraps
  ## labels: the label data
  ## confidence_1: confidence scores from the first classifier
  ## confidence_2: confidence scores from the second classifier
  ## sensitivity_bound: the sensitivity threshold
  ## specificity_bound: the specificity threshold
  ## alpha: significance level for the (1 - alpha) confidence interval
  ##
  ## Output:
  ## CI: confidence interval
  ## obs_stat: observed Delta tpAUC statistic
  ## var_boot: bootstrap variance, 
  ## boot_stat: vector with the bootstrap Delta tpAUC statistics
  
  thr <- c(1 - specificity_bound, sensitivity_bound)
  ## compute the tpAUC values for classifiers 1 and 2 and the Delta tpAUC
  obs_tpAUC_1 <- tpAUC::tproc.est(response = as.factor(labels), 
                                  predictor = confidence_1, 
                                  threshold = thr)
  obs_tpAUC_2 <- tpAUC::tproc.est(response = as.factor(labels), 
                                  predictor = confidence_2, 
                                  threshold = thr)
  
  obs_stat <- obs_tpAUC_1 - obs_tpAUC_2
  
  ## run bootstraps
  n <- length(labels)
  boot_stat <- rep(NA, n_boot)
  for (i in seq(n_boot)) {
    
    ## resample the data with replacement
    idx <- sample(seq(n), n, replace = TRUE)
    boot_labels <- labels[idx]
    boot_conf_1 <- confidence_1[idx]
    boot_conf_2 <- confidence_2[idx]
    
    ## compute delta statistic on the bootstrap sample
    boot_tpAUC_1 <- tpAUC::tproc.est(response = as.factor(boot_labels), 
                                     predictor = boot_conf_1, 
                                     threshold = thr)
    boot_tpAUC_2 <- tpAUC::tproc.est(response = as.factor(boot_labels), 
                                     predictor = boot_conf_2, 
                                     threshold = thr)
    boot_stat[i] <- boot_tpAUC_1 - boot_tpAUC_2
  }
  
  ## compute the bootstrap variance estimate of the Delta statistic
  var_boot <- mean((boot_stat - mean(boot_stat))^2)
  
  ## compute the bootstrap-assisted asymptotic confidence interval
  Z <- qnorm(alpha/2, lower.tail = FALSE)
  CI <- c(obs_stat - Z * sqrt(var_boot/n),
          obs_stat + Z * sqrt(var_boot/n))
  names(CI) <- c("lower_bound", "upper_bound")
  
  list(CI = CI,
       obs_stat = obs_stat,
       var_boot = var_boot,
       boot_stat = boot_stat)
}


####################################
## classification functions 
####################################

## function to fit the random forest classifier
##
FitRangerClass <- function(dat,
                           idx_train, 
                           idx_test, 
                           label_name, 
                           feature_names,
                           neg_class_name, 
                           pos_class_name) {
  ## Inputs:
  ## dat: data.frame containing the features and label data
  ## idx_train: index of the training samples
  ## idx_test: index of the test samples
  ## label_name: name of the outcome variable
  ## feature_names: names of the input variables
  ## neg_class_name: label level for the negative examples
  ## pos_class_name: label level for the positive examples
  ##
  ## Output:
  ## auc_obs: observed AUC score 
  ## pred_probs: vector with the predicted probabilities of the test set
  ##             examples being a positive case 
  ## roc_obj: roc object fit
  
  dat <- dat[, c(label_name, feature_names)]
  dat[, label_name] <- factor(as.character(dat[, label_name]), 
                              levels = c(neg_class_name, pos_class_name)) 
  my_formula <- as.formula(paste(label_name, " ~ ", 
                                 paste(feature_names, collapse = " + ")))
  fit <- ranger(my_formula, 
                data = dat[idx_train,], 
                probability = TRUE, 
                verbose = FALSE)
  pred_probs <- predict(fit, 
                        dat[idx_test, -1, drop = FALSE], 
                        type = "response")$predictions
  y_test <- dat[idx_test, 1]
  roc_obj <- roc(y_test, 
                pred_probs[, pos_class_name], 
                direction = "<", 
                levels = c(neg_class_name, pos_class_name))    
  auc_obs <- pROC::auc(roc_obj)[1]
  
  list(auc_obs = auc_obs, 
       pred_probs = pred_probs[, pos_class_name], 
       roc_obj = roc_obj)
}


## function to fit the logistic regression classifier
##
FitGlmClass <- function(dat,
                        idx_train, 
                        idx_test, 
                        label_name, 
                        feature_names,
                        neg_class_name, 
                        pos_class_name) {
  ## Inputs:
  ## dat: data.frame containing the features and label data
  ## idx_train: index of the training samples
  ## idx_test: index of the test samples
  ## label_name: name of the outcome variable
  ## feature_names: names of the input variables
  ## neg_class_name: label level for the negative examples
  ## pos_class_name: label level for the positive examples
  ##
  ## Output:
  ## auc_obs: observed AUC score 
  ## pred_probs: vector with the predicted probabilities of the test set
  ##             examples being a positive case 
  ## roc_obj: roc object fit
  
  dat <- dat[, c(label_name, feature_names)]
  dat[, label_name] <- factor(as.character(dat[, label_name]), 
                              levels = c(neg_class_name, pos_class_name)) 
  my_formula <- as.formula(paste(label_name, " ~ ", 
                                 paste(feature_names, collapse = " + ")))
  fit <- glm(my_formula, 
             data = dat[idx_train,], 
             family = "binomial")
  pred_probs <- predict(fit, 
                        dat[idx_test, -1, drop = FALSE], 
                        type = "response")
  y_test <- dat[idx_test, 1]
  roc_obj <- roc(y_test, 
                 pred_probs, 
                 direction = "<", 
                 levels = c(neg_class_name, pos_class_name))    
  auc_obs <- pROC::auc(roc_obj)[1]
  
  list(auc_obs = auc_obs, 
       pred_probs = pred_probs, 
       roc_obj = roc_obj)
}

