## source the utility functions
source("two_way_partial_AUC_utility_functions.R")

data_path <- ""

##########################################
## Diabetic retinopathy data
##########################################

## get the data from the UCI Machine Learning database
fname_1 <- url("http://archive.ics.uci.edu/ml/machine-learning-databases/00329/messidor_features.arff")
dat_1 <- read.csv(fname_1, header=FALSE, comment.char = "@")

dim(dat_1)

## the last column contains the label we are interested in predicting
## and we rename it as "Class"
names(dat_1)[20] <- "Class" 


## randomly split the data indexes into training and test set indexes 
set.seed(1234567)
idx_train_1 <- sample(seq(nrow(dat_1)), round(nrow(dat_1)/2), replace = FALSE)
idx_test_1 <- setdiff(seq(nrow(dat_1)), idx_train_1)


## fit the random forest classifier
rf_1 <- FitRangerClass(dat_1,
                     idx_train_1, 
                     idx_test_1, 
                     label_name = "Class", 
                     feature_names = names(dat_1)[-c(20)],
                     neg_class_name = "0", 
                     pos_class_name = "1")

## fit the logistic regression classifier
lr_1 <- FitGlmClass(dat_1,
                  idx_train_1, 
                  idx_test_1, 
                  label_name = "Class", 
                  feature_names = names(dat_1)[-c(20)],
                  neg_class_name = "0", 
                  pos_class_name = "1")

rf_1$auc_obs
lr_1$auc_obs

## compute the bootstrap assisted asymptotic confidence interval
## using the proposed estimator
set.seed(1234567)
out_1_p <- pROCBasedDeltaCITime(n_boot = 1000,
                              labels = dat_1$Class[idx_test_1],
                              confidence_1 = lr_1$pred_probs,
                              confidence_2 = rf_1$pred_probs,
                              sensitivity_bound = 0.4, 
                              specificity_bound = 0.4,
                              alpha = 0.05)

out_1_p$total_computing_time
lapply(out_1_p$ci_stats[1:3], round, 6)
lapply(out_1_p$obs_stats, round, 6)

## alternative way to compute the bootstrap assisted asymptotic
## confidence interval using the proposed estimator
## (instead of recording the time of each computation and
## adding them all the get the total time, here we record the
## total time taken to run all the calculations directly)
set.seed(1234567)
system.time(out_1_p2 <- pROCBasedDeltaCI(n_boot = 1000,
                                        labels = dat_1$Class[idx_test_1],
                                        confidence_1 = lr_1$pred_probs,
                                        confidence_2 = rf_1$pred_probs,
                                        sensitivity_bound = 0.4, 
                                        specificity_bound = 0.4,
                                        alpha = 0.05))

round(out_1_p2$CI, 6)
round(out_1_p2$obs_stat, 6)

## compute the bootstrap assisted asymptotic confidence interval
## using the original estimator
set.seed(1234567)
out_1_o <- tpAUCBasedDeltaCITime(n_boot = 1000,
                               labels = dat_1$Class[idx_test_1],
                               confidence_1 = lr_1$pred_probs,
                               confidence_2 = rf_1$pred_probs,
                               sensitivity_bound = 0.4, 
                               specificity_bound = 0.4,
                               alpha = 0.05)

out_1_o$total_computing_time
lapply(out_1_o$ci_stats[1:3], round, 6)
lapply(out_1_o$obs_stats, round, 6)

## alternative way to compute the bootstrap assisted asymptotic
## confidence interval using the original estimator
## (instead of recording the time of each computation and
## adding them all the get the total time, here we record the
## total time taken to run all the calculations directly)
set.seed(1234567)
system.time(out_1_o2 <- tpAUCBasedDeltaCI(n_boot = 1000,
                                         labels = dat_1$Class[idx_test_1],
                                         confidence_1 = lr_1$pred_probs,
                                         confidence_2 = rf_1$pred_probs,
                                         sensitivity_bound = 0.4, 
                                         specificity_bound = 0.4,
                                         alpha = 0.05))

round(out_1_o2$CI, 6)
round(out_1_o2$obs_stat, 6)


##########################################
## Sepsis survival data
##########################################

## To get the sepsis survival data, download the zip file from the UCI database:
## "http://archive.ics.uci.edu/ml/machine-learning-databases/00628/s41598-020-73558-3_sepsis_survival_dataset.zip
## and then save the csv file s41598-020-73558-3_sepsis_survival_primary_cohort.csv

## read the saved csv file
data_path <- ""
dat_2 <- read.csv(paste0(data_path, "s41598-020-73558-3_sepsis_survival_primary_cohort.csv"), header = TRUE)

dim(dat_2)

## keep only the first record of sepsis of each individual
dat_2 <- dat_2[dat_2$episode_number == 1,] 
dim(dat_2)

## remove the episode_number column
dat_2 <- dat_2[, -3]

## randomly split the data indexes into training and test set indexes 
set.seed(1234567)
idx_train_2 <- sample(seq(nrow(dat_2)), round(nrow(dat_2)/2), replace = FALSE)
idx_test_2 <- setdiff(seq(nrow(dat_2)), idx_train_2)

## fit the random forest classifier
rf_2 <- FitRangerClass(dat_2,
                     idx_train_2, 
                     idx_test_2, 
                     label_name = "hospital_outcome_1alive_0dead", 
                     feature_names = names(dat_2)[-c(3)],
                     neg_class_name = "0", 
                     pos_class_name = "1")

## fit the logistic regression classifier
lr_2 <- FitGlmClass(dat_2,
                  idx_train_2, 
                  idx_test_2, 
                  label_name = "hospital_outcome_1alive_0dead", 
                  feature_names = names(dat_2)[-c(3)],
                  neg_class_name = "0", 
                  pos_class_name = "1")

rf_2$auc_obs
lr_2$auc_obs

## compute the bootstrap assisted asymptotic confidence interval
## using the proposed estimator
set.seed(1234567)
out_2_p <- pROCBasedDeltaCITime(n_boot = 1000,
                              labels = dat_2$hospital_outcome_1alive_0dead[idx_test_2],
                              confidence_1 = lr_2$pred_probs,
                              confidence_2 = rf_2$pred_probs,
                              sensitivity_bound = 0.4, 
                              specificity_bound = 0.4,
                              alpha = 0.05)
out_2_p$total_computing_time
lapply(out_2_p$ci_stats[1:3], round, 6)
lapply(out_2_p$obs_stats, round, 6)


## compute the bootstrap assisted asymptotic confidence interval
## using the original estimator (because this take a very long time
## to compute this function saves the results at each bootstrap iteration)
f_name <- "output_sepsis_survival_original_estimator.RData"
set.seed(1234567)
out_2_o <- tpAUCBasedDeltaCITimeSave(n_boot = 1000,
                               labels = dat_2$hospital_outcome_1alive_0dead[idx_test_2],
                               confidence_1 = lr_2$pred_probs,
                               confidence_2 = rf_2$pred_probs,
                               sensitivity_bound = 0.4, 
                               specificity_bound = 0.4,
                               alpha = 0.05,
                               file_name = f_name)

out_2_o$total_computing_time
out_2_o$total_computing_time[1]/(60*60*24) ## get computing time in days
lapply(out_2_o$ci_stats[1:3], round, 6)
lapply(out_2_o$obs_stats, round, 6)


##########################################
## generate Figure 10
##########################################

spB <- 0.4
seB <- 0.4
lwd_roc <- 1
cp <- 1.3

fig_path <- ""

#pdf(paste0(fig_path, "figure10.pdf"), width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 1.5, 0.5), mgp = c(2.5, 0.5, 0))
plot(rf_1$roc_obj, col = "purple", lwd = lwd_roc,
     main = "Diabetic Retinopathy")
plot(lr_1$roc_obj, col = "black", lwd = lwd_roc, add = TRUE)
legend("bottomright", legend = c("logistic regression", "random forest"),
       text.col = c("black", "purple"), bty = "n")
segments(x0 = 1, x1 = spB, y0 = seB, y1 = seB, col = "red", lwd = 1)
segments(x0 = 1, x1 = spB, y0 = 1, y1 = 1, col = "red", lwd = 1)
segments(x0 = spB, x1 = spB, y0 = seB, y1 = 1, col = "red", lwd = 1)
segments(x0 = 1, x1 = 1, y0 = seB, y1 = 1, col = "red", lwd = 1)
mtext("(a)", side = 3, at = 1.1, cex = cp)
####
plot(rf_2$roc_obj, col = "purple", lwd = lwd_roc,
     main = "Sepsis Survival")
plot(lr_2$roc_obj, col = "black", lwd = lwd_roc, add = TRUE)
legend("bottomright", legend = c("logistic regression", "random forest"),
       text.col = c("black", "purple"), bty = "n")
segments(x0 = 1, x1 = spB, y0 = seB, y1 = seB, col = "red", lwd = 1)
segments(x0 = 1, x1 = spB, y0 = 1, y1 = 1, col = "red", lwd = 1)
segments(x0 = spB, x1 = spB, y0 = seB, y1 = 1, col = "red", lwd = 1)
segments(x0 = 1, x1 = 1, y0 = seB, y1 = 1, col = "red", lwd = 1)
mtext("(b)", side = 3, at = 1.1, cex = cp)
#dev.off()

