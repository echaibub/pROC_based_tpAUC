
## source the utility functions
source("two_way_partial_AUC_utility_functions.R")

#######################################
## run first additional experiment
#######################################

## set sample size grid
n_grid <- c(seq(1e+2, 9e+2, by = 1e+2),
            seq(1e+3, 9e+3, by = 1e+3),
            seq(1e+4, 2e+4, by = 2e+3))

## number of replications
n_runs <- 100

## set parameter values for binormal ROC model
mu_N <- 0
sig_N <- 1
mu_P <- 2
sig_P <- 1.5

## set the sensitivity and specificity thresholds
sensitivity_bound <- 0.6
specificity_bound <- 0.6

## set object for storing the results
outputs <- vector(mode = "list", length = length(n_grid))

## run the experiments
set.seed(123456789)
for (i in seq(n_grid)) {
  cat("grid", i, "\n")
  aux <- ComputationTimeSimulations(n_runs,
                                    n_N = n_grid[i]/2,
                                    n_P = n_grid[i]/2,
                                    mu_N,
                                    sig_N,
                                    mu_P,
                                    sig_P,
                                    sensitivity_bound,
                                    specificity_bound)
  
  outputs[[i]] <- aux
}



#######################################
## run second additional experiment
#######################################

## set sample size grid
n_grid_2 <- seq(1e+4, 1e+5, by = 1e+4)

## create objects that will store the results
user_time <- matrix(NA, length(n_grid_2), 3)
rownames(user_time) <- n_grid_2
colnames(user_time) <- c("pROC_1", "pROC_2", "tpAUC")
system_time <- user_time
elapsed_time <- user_time

## set parameter values for binormal ROC model
mu_N <- 0
sig_N <- 1
mu_P <- 2
sig_P <- 1.5

## set the sensitivity and specificity thresholds
sensitivity_bound <- 0.6
specificity_bound <- 0.6

## run the experiments
set.seed(123456789)
for (i in seq(n_grid_2)) {
  cat("grid", i, "\n")
  aux <- ComputationTimeSimulations(n_runs = 1,
                                    n_N = n_grid_2[i]/2,
                                    n_P = n_grid_2[i]/2,
                                    mu_N,
                                    sig_N,
                                    mu_P,
                                    sig_P,
                                    sensitivity_bound,
                                    specificity_bound)
  ## user time
  user_time[i, "pROC_1"] <- aux$computation_time_pROC_1[1, 1]
  user_time[i, "pROC_2"] <- aux$computation_time_pROC_2[1, 1]
  user_time[i, "tpAUC"] <- aux$computation_time_tpAUC[1, 1]

  ## system time
  system_time[i, "pROC_1"] <- aux$computation_time_pROC_1[1, 2]
  system_time[i, "pROC_2"] <- aux$computation_time_pROC_2[1, 2]
  system_time[i, "tpAUC"] <- aux$computation_time_tpAUC[1, 2]
  
  ## elapsed time
  elapsed_time[i, "pROC_1"] <- aux$computation_time_pROC_1[1, 3]
  elapsed_time[i, "pROC_2"] <- aux$computation_time_pROC_2[1, 3]
  elapsed_time[i, "tpAUC"] <- aux$computation_time_tpAUC[1, 3]
}


###############################################
## generate Figures 7 and 8
###############################################

## organize the outputs of the first experiment
user_time_pROC_1 <- ShapeUserTimepROC_1(x = outputs, n_grid)
user_time_pROC_2 <- ShapeUserTimepROC_2(x = outputs, n_grid)
user_time_tpAUC <- ShapeUserTimetpAUC(x = outputs, n_grid)

cp <- 1.5
cl <- 1.3

fig_path <- ""

#pdf(paste0(fig_path, "figure7.pdf"), width = 10, height = 5.25)
par(mfrow = c(1, 2), mar = c(4, 4, 1.5, 0.5), mgp = c(2.5, 0.5, 0))
boxplot(user_time_tpAUC,
        las = 2, 
        at = seq(24) - 0.1,
        col = "white", 
        border = rgb(1, 0, 0, 0.5),
        ylab = "computation time (seconds)",
        xlab = "sample size",
        main = "Computation time comparison",
        cex = 0.5,
        cex.lab = cl)
boxplot(user_time_pROC_1, 
        las = 2, 
        xaxt = "n",
        at = seq(24) + 0.1,
        col = "white", 
        border = rgb(0, 0, 1, 0.5),
        cex = 0.5,
        add = TRUE)
legend("topleft", legend = c("original estimator", "proposed estimator"),
       text.col = c("red", "darkblue"), bty = "n", cex = cl)
mtext("(a)", side = 3, at = -1, cex = cp)
####
boxplot(user_time_pROC_1, 
        las = 2, 
        at = seq(24) - 0.1, 
        col = "white", 
        border = rgb(0, 0, 1, 0.5),
        ylab = "computation time (seconds)",
        xlab = "sample size",
        main = "Computation time comparison",
        cex = 0.5,
        cex.lab = cl)
boxplot(user_time_pROC_2, 
        xaxt = "n",
        at = seq(24) + 0.1, 
        las = 2, 
        col = "white", 
        border = rgb(0, 1, 0, 0.5),
        cex = 0.5,
        add = TRUE)
lines(apply(user_time_pROC_1, 2, mean), col = "darkblue", type = "l", lwd = 2)
lines(apply(user_time_pROC_2, 2, mean), col = "darkgreen", type = "l", lwd = 2)
legend("topleft", legend = c("proposed estimator", "alternative estimator"),
       text.col = c("darkblue", "darkgreen"), bty = "n", cex = cl)
mtext("(b)", side = 3, at = -1, cex = cp)
#dev.off()


#pdf(paste0(fig_path, "figure8.pdf"), width = 10, height = 5.25)
par(mfrow = c(1, 2), mar = c(4, 4, 1.5, 0.5), mgp = c(2.75, 0.5, 0))
plot(n_grid_2, user_time[, "tpAUC"],
     xaxt = "n",
     col = "red",
     ylab = "computation time (seconds)",
     xlab = "sample size",
     main = "Computation time comparison",
     cex.lab = cl,
     type = "b",
     lwd = 2)
axis(side = 1, at = n_grid_2, labels = n_grid_2, las = 2)
lines(n_grid_2, user_time[, "pROC_1"], col = "darkblue", lwd = 2, type = "b")
legend("topleft", legend = c("original estimator", "proposed estimator"),
       text.col = c("red", "darkblue"), bty = "n", cex = cl)
mtext("(a)", side = 3, at = -1, cex = cp)
####
plot(n_grid_2, user_time[, "pROC_1"],
     xaxt = "n",
     col = "darkblue",
     ylab = "computation time (seconds)",
     xlab = "sample size",
     main = "Computation time comparison",
     cex.lab = cl,
     type = "b",
     lwd = 2)
axis(side = 1, at = n_grid_2, labels = n_grid_2, las = 2)
lines(n_grid_2, user_time[, "pROC_2"], col = "darkgreen", lwd = 2, type = "b")
legend("topleft", legend = c("proposed estimator", "alternative estimator"),
       text.col = c("darkblue", "darkgreen"), bty = "n", cex = cl)
mtext("(b)", side = 3, at = -1, cex = cp)
#dev.off()

