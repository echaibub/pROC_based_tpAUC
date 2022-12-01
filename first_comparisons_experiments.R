
## source the utility functions
source("two_way_partial_AUC_utility_functions.R")


##############
## n = 100
##############

set.seed(12345)
e1 <- RunFirstComparisons(n_N = 90, n_P = 10, n_runs = 1000)

set.seed(12345)
e2 <- RunFirstComparisons(n_N = 50, n_P = 50, n_runs = 1000)

set.seed(12345)
e3 <- RunFirstComparisons(n_N = 10, n_P = 90, n_runs = 1000)


##############
## n = 1000
##############

set.seed(12345)
e4 <- RunFirstComparisons(n_N = 900, n_P = 100, n_runs = 1000)

set.seed(12345)
e5 <- RunFirstComparisons(n_N = 500, n_P = 500, n_runs = 1000)

set.seed(12345)
e6 <- RunFirstComparisons(n_N = 100, n_P = 900, n_runs = 1000)


##############
## n = 10000
##############

set.seed(12345)
e7 <- RunFirstComparisons(n_N = 9000, n_P = 1000, n_runs = 1000)

set.seed(12345)
e8 <- RunFirstComparisons(n_N = 5000, n_P = 5000, n_runs = 1000)

set.seed(12345)
e9 <- RunFirstComparisons(n_N = 1000, n_P = 9000, n_runs = 1000)


#################################
## generate Figures 5 and 6
#################################

cp <- 1.75
cl <- 1.3
ca <- 1.3

fig_path <- ""

#pdf(paste0(fig_path, "figure5.pdf"), width = 8.5, height = 9.5)
par(mfrow = c(3, 3), mar = c(4, 3, 1.5, 0.5), mgp = c(1.75, 0.5, 0))
plot(e1$tpaucs[, 1], e1$tpaucs[, 2], 
     xlab = "original estimator", 
     ylab = "proposed estimator",
     main = "neg.cases = 90, pos. cases = 10",
     cex.lab = cl, cex.axis = ca)
abline(a = 0, b = 1, col = "red")
legend("topleft", legend = "(a)", cex = cp, bty = "n")
####
plot(e2$tpaucs[, 1], e2$tpaucs[, 2], 
     xlab = "original estimator", 
     ylab = "proposed estimator",
     main = "neg.cases = 50, pos. cases = 50",
     cex.lab = cl, cex.axis = ca)
abline(a = 0, b = 1, col = "red")
legend("topleft", legend = "(b)", cex = cp, bty = "n")
####
plot(e3$tpaucs[, 1], e3$tpaucs[, 2], 
     xlab = "original estimator", 
     ylab = "proposed estimator",
     main = "neg.cases = 10, pos. cases = 90",
     cex.lab = cl, cex.axis = ca)
abline(a = 0, b = 1, col = "red")
legend("topleft", legend = "(c)", cex = cp, bty = "n")
####
####
plot(e4$tpaucs[, 1], e4$tpaucs[, 2], 
     xlab = "original estimator", 
     ylab = "proposed estimator",
     main = "neg.cases = 900, pos. cases = 100",
     cex.lab = cl, cex.axis = ca)
abline(a = 0, b = 1, col = "red")
legend("topleft", legend = "(d)", cex = cp, bty = "n")
####
plot(e5$tpaucs[, 1], e5$tpaucs[, 2], 
     xlab = "original estimator", 
     ylab = "proposed estimator",
     main = "neg.cases = 500, pos. cases = 500",
     cex.lab = cl, cex.axis = ca)
abline(a = 0, b = 1, col = "red")
legend("topleft", legend = "(e)", cex = cp, bty = "n")
####
plot(e6$tpaucs[, 1], e6$tpaucs[, 2], 
     xlab = "original estimator", 
     ylab = "proposed estimator",
     main = "neg.cases = 100, pos. cases = 900",
     cex.lab = cl, cex.axis = ca)
abline(a = 0, b = 1, col = "red")
legend("topleft", legend = "(f)", cex = cp, bty = "n")
####
####
plot(e7$tpaucs[, 1], e7$tpaucs[, 2], 
     xlab = "original estimator", 
     ylab = "proposed estimator",
     main = "neg.cases = 9000, pos. cases = 1000",
     cex.lab = cl, cex.axis = ca)
abline(a = 0, b = 1, col = "red")
legend("topleft", legend = "(g)", cex = cp, bty = "n")
####
plot(e8$tpaucs[, 1], e8$tpaucs[, 2], 
     xlab = "original estimator", 
     ylab = "proposed estimator",
     main = "neg.cases = 5000, pos. cases = 5000",
     cex.lab = cl, cex.axis = ca)
abline(a = 0, b = 1, col = "red")
legend("topleft", legend = "(h)", cex = cp, bty = "n")
####
plot(e9$tpaucs[, 1], e9$tpaucs[, 2], 
     xlab = "original estimator", 
     ylab = "proposed estimator",
     main = "neg.cases = 1000, pos. cases = 9000",
     cex.lab = cl, cex.axis = ca)
abline(a = 0, b = 1, col = "red")
legend("topleft", legend = "(i)", cex = cp, bty = "n")
par(mfrow = c(1, 1))
#dev.off()


e.o <- data.frame(e1$computation_time_tpAUC[, 1], 
                  e2$computation_time_tpAUC[, 1],
                  e3$computation_time_tpAUC[, 1],
                  e4$computation_time_tpAUC[, 1],
                  e5$computation_time_tpAUC[, 1],
                  e6$computation_time_tpAUC[, 1],
                  e7$computation_time_tpAUC[, 1],
                  e8$computation_time_tpAUC[, 1],
                  e9$computation_time_tpAUC[, 1])

e.p <- data.frame(e1$computation_time_pROC[, 1], 
                  e2$computation_time_pROC[, 1],
                  e3$computation_time_pROC[, 1],
                  e4$computation_time_pROC[, 1],
                  e5$computation_time_pROC[, 1],
                  e6$computation_time_pROC[, 1],
                  e7$computation_time_pROC[, 1],
                  e8$computation_time_pROC[, 1],
                  e9$computation_time_pROC[, 1])


#pdf(paste0(fig_path, "figure6.pdf"), width = 8, height = 4.5)
par(mfrow = c(1, 2), mar = c(5, 4, 1.5, 0.5), mgp = c(2, 0.5, 0))
boxplot(e.o, 
        names = paste0("exp. ", seq(9)), 
        ylab = "computation time (seconds)", 
        main = "original", 
        border = rep(c("darkorange", "purple", "black"), 3),
        las = 2, cex = 0.5, col = "white")
mtext(side = 1, at = 2, line = 3, text = "n = 100")
mtext(side = 1, at = 5, line = 3, text = "n = 1,000")
mtext(side = 1, at = 8, line = 3, text = "n = 10,000")
mtext(side = 1, at = 5, line = 4, text = "sample size")
abline(v = 3.5, col = "grey")
abline(v = 6.5, col = "grey")
mtext("(a)", side = 3, at = 0.5, line = 0.25, cex = 1.5)
legend("topleft", 
       legend = c("90% neg, 10% pos", "50% neg, 50% pos", "10% neg, 90% pos"),
       text.col = c("darkorange", "purple", "black"))
boxplot(e.p, 
        names = paste0("exp. ", seq(9)), 
        ylab = "computation time (seconds)", 
        main = "proposed", 
        border = rep(c("darkorange", "purple", "black"), 3),
        las = 2, cex = 0.5, col = "white")
mtext(side = 1, at = 2, line = 3, text = "n = 100")
mtext(side = 1, at = 5, line = 3, text = "n = 1,000")
mtext(side = 1, at = 8, line = 3, text = "n = 10,000")
mtext(side = 1, at = 5, line = 4, text = "sample size")
abline(v = 3.5, col = "grey")
abline(v = 6.5, col = "grey")
mtext("(b)", side = 3, at = 0.5, line = 0.25, cex = 1.5)
legend("topleft", 
       legend = c("90% neg, 10% pos", "50% neg, 50% pos", "10% neg, 90% pos"),
       text.col = c("darkorange", "purple", "black"))
#dev.off()

