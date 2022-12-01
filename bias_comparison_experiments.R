
## source the utility functions
source("two_way_partial_AUC_utility_functions.R")

## set number of replications
n_sim <- 100

## set random seeds
set.seed(12345)
my_seeds <- sample(1e+4, 1e+5, n_sim)

##############################
## n = 100
##############################

bc.100.1 <- RunBiasComparisons(n_sim = n_sim,
                               n_runs = 1000,
                               n_N = 90,
                               n_P = 10,
                               mu_N_range = c(0, 1), 
                               sig_N_range = c(1, 2), 
                               mu_P_range = c(2, 3), ## needs to be above mu_N 
                               sig_P_range = c(1, 2),
                               sensitivity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               specificity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               my_seeds = my_seeds)


bc.100.2 <- RunBiasComparisons(n_sim = n_sim,
                               n_runs = 1000,
                               n_N = 50,
                               n_P = 50,
                               mu_N_range = c(0, 1), 
                               sig_N_range = c(1, 2), 
                               mu_P_range = c(2, 3), ## needs to be above mu_N 
                               sig_P_range = c(1, 2),
                               sensitivity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               specificity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               my_seeds = my_seeds)


bc.100.3 <- RunBiasComparisons(n_sim = n_sim,
                               n_runs = 1000,
                               n_N = 10,
                               n_P = 90,
                               mu_N_range = c(0, 1), 
                               sig_N_range = c(1, 2), 
                               mu_P_range = c(2, 3), ## needs to be above mu_N 
                               sig_P_range = c(1, 2),
                               sensitivity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               specificity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               my_seeds = my_seeds)



##############################
## n = 1000
##############################

bc.1000.1 <- RunBiasComparisons(n_sim = n_sim,
                               n_runs = 1000,
                               n_N = 900,
                               n_P = 100,
                               mu_N_range = c(0, 1), 
                               sig_N_range = c(1, 2), 
                               mu_P_range = c(2, 3), ## needs to be above mu_N 
                               sig_P_range = c(1, 2),
                               sensitivity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               specificity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               my_seeds = my_seeds)


bc.1000.2 <- RunBiasComparisons(n_sim = n_sim,
                               n_runs = 1000,
                               n_N = 500,
                               n_P = 500,
                               mu_N_range = c(0, 1), 
                               sig_N_range = c(1, 2), 
                               mu_P_range = c(2, 3), ## needs to be above mu_N 
                               sig_P_range = c(1, 2),
                               sensitivity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               specificity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               my_seeds = my_seeds)


bc.1000.3 <- RunBiasComparisons(n_sim = n_sim,
                               n_runs = 1000,
                               n_N = 100,
                               n_P = 900,
                               mu_N_range = c(0, 1), 
                               sig_N_range = c(1, 2), 
                               mu_P_range = c(2, 3), ## needs to be above mu_N 
                               sig_P_range = c(1, 2),
                               sensitivity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               specificity_bound_grid = c(0.2, 0.4, 0.6, 0.8),
                               my_seeds = my_seeds)



#######################################
## generate Figure 9
#######################################

bc1 <- bc.100.1
bc2 <- bc.100.2
bc3 <- bc.100.3
bc4 <- bc.1000.1
bc5 <- bc.1000.2
bc6 <- bc.1000.3

cl <- 1.5
cp <- 1.3

fig_path <- ""

#pdf(paste0(fig_path, "figure9.pdf"), width = 9.5, height = 7)
par(mfrow = c(2, 3), mar = c(7, 5.5, 1.5, 0.5), mgp = c(3, 0.5, 0))
boxplot(abs(bc1$tpAUC_bias), col = rgb(1, 0, 0, 0.5), las = 2,
        ylab = "absolute value of estimated bias", 
        main = "neg. cases = 90, pos. cases = 10")
boxplot(abs(bc1$pROC_bias), col = rgb(0, 0, 1, 0.5), las = 2, add = TRUE)
mtext("(sensitivity, specificity)", side = 1, line = 5, cex = 1)
mtext("(a)", side = 3, at = -1, cex = cp)
legend("topleft", legend = c("original estimator", "proposed estimator"),
       text.col = c("red", "blue"), bty = "n", cex = cl)
####
boxplot(abs(bc2$tpAUC_bias), col = rgb(1, 0, 0, 0.5), las = 2,
        ylab = "absolute value of estimated bias", 
        main = "neg. cases = 50, pos. cases = 50")
boxplot(abs(bc2$pROC_bias), col = rgb(0, 0, 1, 0.5), las = 2, add = TRUE)
mtext("(sensitivity, specificity)", side = 1, line = 5, cex = 1)
mtext("(b)", side = 3, at = -1, cex = cp)
legend("topright", legend = c("original estimator", "proposed estimator"),
       text.col = c("red", "blue"), bty = "n", cex = cl)
####
boxplot(abs(bc3$tpAUC_bias), col = rgb(1, 0, 0, 0.5), las = 2,
        ylab = "absolute value of estimated bias", 
        main = "neg. cases = 10, pos. cases = 90")
boxplot(abs(bc3$pROC_bias), col = rgb(0, 0, 1, 0.5), las = 2, add = TRUE)
mtext("(sensitivity, specificity)", side = 1, line = 5, cex = 1)
mtext("(c)", side = 3, at = -1, cex = cp)
legend("topright", legend = c("original estimator", "proposed estimator"),
       text.col = c("red", "blue"), bty = "n", cex = cl)
####
####
boxplot(abs(bc4$tpAUC_bias), col = rgb(1, 0, 0, 0.5), las = 2,
        ylab = "absolute value of estimated bias", 
        main = "neg. cases = 900, pos. cases = 100")
boxplot(abs(bc4$pROC_bias), col = rgb(0, 0, 1, 0.5), las = 2, add = TRUE)
mtext("(sensitivity, specificity)", side = 1, line = 5, cex = 1)
mtext("(d)", side = 3, at = -1, cex = cp)
legend("topleft", legend = c("original estimator", "proposed estimator"),
       text.col = c("red", "blue"), bty = "n", cex = cl)
####
boxplot(abs(bc5$tpAUC_bias), col = rgb(1, 0, 0, 0.5), las = 2,
        ylab = "absolute value of estimated bias", 
        main = "neg. cases = 500, pos. cases = 500")
boxplot(abs(bc5$pROC_bias), col = rgb(0, 0, 1, 0.5), las = 2, add = TRUE)
mtext("(sensitivity, specificity)", side = 1, line = 5, cex = 1)
mtext("(e)", side = 3, at = -1, cex = cp)
legend("topright", legend = c("original estimator", "proposed estimator"),
       text.col = c("red", "blue"), bty = "n", cex = cl)
####
boxplot(abs(bc6$tpAUC_bias), col = rgb(1, 0, 0, 0.5), las = 2,
        ylab = "absolute value of estimated bias", 
        main = "neg. cases = 100, pos. cases = 900")
boxplot(abs(bc6$pROC_bias), col = rgb(0, 0, 1, 0.5), las = 2, add = TRUE)
mtext("(sensitivity, specificity)", side = 1, line = 5, cex = 1)
mtext("(f)", side = 3, at = -1, cex = cp)
legend("topright", legend = c("original estimator", "proposed estimator"),
       text.col = c("red", "blue"), bty = "n", cex = cl)
par(mfrow = c(1, 1))
#dev.off()
