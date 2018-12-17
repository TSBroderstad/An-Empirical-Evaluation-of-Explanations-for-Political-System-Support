########################################################################################
#################################### SVF - 8009 ########################################
########################################################################################

library(haven)
library(naniar)
library(mlr)
library(mlbench)
library(kernlab)
library(ranger)
library(cmaes)
library(mmpf)
library(haven)

ESS6 <- read.csv("ESS6_cleaned")

#### LOOP FOR VAR. IMPORTANCE ####

#PREAMBLE LOOP  
cntry <- as.character(unique(ESS6$cntry))
regr.task <- list()
mod <- list()
lrn <- list()
importance <- list()
set.seed(1234)
lrn <- makeLearner("regr.ranger", importance = c("permutation"))

## START LOOP
for (i in 1: length(cntry)){
  data <- ESS6
  data <- data[data$cntry == cntry[i], ]
  
  
  
  # DEFINE TASK 
  regr.task[[i]] <-  makeRegrTask(data = data, target = "stfdem")
  
  # FEATURE IMPORTANCE 
  mod[[i]] <- train(lrn, regr.task[[i]])
  importance[[i]] <- getFeatureImportance(mod[[i]], nmc = 50L, interaction = TRUE, local = TRUE)
  importance[[i]] <- gather(importance[[i]]$res, "var.name", "importance")
  # importance[[i]] <- rename(importance[[i]][,1:2], c("importance" = cntry[i]))
} ## END LOOP

# MAKE IMPORTANCE DF
importance2 <- data.frame(importance[1])
for (i in 2:length(cntry)){
  importance2[,i+1] <- importance[[i]]$importance
} ## END LOOP

#POST
colnames(importance2) <- c("var.name", cntry)

importance2 <- importance2[-1,]

##########################################################################
############################### RMSE FOR COUNTRY ########################
##########################################################################

#### TUNING ####
# RUN ONCE

set.seed(1234)
cntry <- as.character(unique(ESS6$cntry))
regr.task <- list()
bmr <- list()
perf <- list()

# Define a search space for each learner's parameter
ps_ksvm = makeParamSet(
  makeNumericParam("sigma", lower = -12, upper = 12, trafo = function(x) 2^x)
)
ps_rf = makeParamSet(
  makeIntegerParam("num.trees", lower = 1L, upper = 200L)
)
# Choose a resampling strategy
rdesc = makeResampleDesc("CV", iters = 10)
# Choose a performance measure
meas = rmse
# Choose a tuning method
ctrl = makeTuneControlCMAES(budget = 100L)
# Make tuning wrappers
tuned.ksvm = makeTuneWrapper(learner = "regr.ksvm", resampling = rdesc, measures = meas,
                             par.set = ps_ksvm, control = ctrl, show.info = FALSE)
tuned.rf = makeTuneWrapper(learner = "regr.ranger", resampling = rdesc, measures = meas,
                           par.set = ps_rf, control = ctrl, show.info = FALSE)
# Four learners to be compared
lrns = list(makeLearner("regr.lm"), tuned.ksvm,
            tuned.rf)


## START LOOP
for (i in 1: length(cntry)){
  data <- ESS6
  data <- data[data$cntry == cntry[i], ]
  data$cntry <- NULL
  
  # Define task
  regr.task[[i]] <-  makeRegrTask(data = data, target = "stfdem")
  
  # Conduct the benchmark experiment
  bmr[[i]] <-  benchmark(learners = lrns, tasks = regr.task[[i]], 
                         resamplings = rdesc, measures = rmse, show.info = FALSE,
                         keep.pred = FALSE, models = FALSE)
}


#### EXTRACT AGG. RMSE FOR COUNTRY LEVEL ####
regr.lm <- vector()
regr.ranger <- vector()
regr.ksvm <- vector()
for (i in 1: length(bmr)){
  
  regr.lm <- c(regr.lm, bmr[[i]][["results"]][["data"]][["regr.lm"]][["aggr"]][["rmse.test.rmse"]])
  regr.ranger <- c(regr.ranger, bmr[[i]][["results"]][["data"]][["regr.ksvm.tuned"]][["aggr"]][["rmse.test.rmse"]])
  regr.ksvm <- c(regr.ksvm, bmr[[i]][["results"]][["data"]][["regr.ranger.tuned"]][["aggr"]][["rmse.test.rmse"]])
}

# Build DF for RMSE
rmse_df <- as.data.frame(cbind(regr.lm, regr.ranger, regr.ksvm))


rmse_ci <- as.data.frame(rmse_ci)
rmse_ci <- t(rmse_ci)

rmse_ci$variable <-  rownames(rmse_ci)

rmse_ci <- as.data.frame(rmse_ci)



colnames(rmse_ci) <- rownames(rmse_ci)
rownames(rmse_ci) <- colnames(rmse_ci)
rmse_ci$variable <- rownames(rmse_ci)

rmse_ci$variable <- NULL

#### PLOT RMSE FOR COUNTRY ####


rmse_plot <- rmse_ci %>%
  ggplot(aes (x = Mean, y = reorder(variable, Mean))) +
  geom_point() +
  geom_errorbarh(aes(xmax = UpperCI, xmin = LowerCI), height = 0.1) +
  ylab("") +
  xlab("RMSE") +
  geom_vline(xintercept = 1.722, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0,2, 0.05),
                     limits = c(1.5,2))  +
  theme_m()
rmse_plot 




#### EXTRACT RMSE FOR IND. CNRTY ####

lm <- vector()
ranger <- vector()
ksvm <- vector()
for (i in 1: length(bmr)){
  
  lm <- c(lm, bmr[[i]][["results"]][["data"]][["regr.lm"]][["measures.test"]][["rmse"]])
  ranger <- c(ranger, bmr[[i]][["results"]][["data"]][["regr.ranger.tuned"]][["measures.test"]][["rmse"]])
  ksvm <- c(ksvm, bmr[[i]][["results"]][["data"]][["regr.ksvm.tuned"]][["measures.test"]][["rmse"]])
}
# Build DF for RMSE
rmsedf <- as.data.frame(cbind(lm, ranger, ksvm))

rmsedf$id <- 1:nrow(rmsedf)

rmsedf <- rmsedf %>% 
  mutate(
    cntry = case_when(
      id == 1:10 ~ "AL",
      id == 11:20 ~ "BE",
      id == 21:30 ~ "BG",
      id == 31:40 ~ "CY",
      id == 41:50 ~ "CZ",
      id == 51:60 ~ "DK",
      id == 61:70 ~ "EE",
      id == 71:80 ~ "FI",
      id == 81:90 ~ "FR",
      id == 91:100 ~ "DE",
      id == 101:110 ~ "HU",
      id == 111:120 ~ "IS",
      id == 121:130 ~ "IE",
      id == 131:140 ~ "IL",
      id == 141:150 ~ "IT",
      id == 151:160 ~ "XK",
      id == 161:170 ~ "LT",
      id == 171:180 ~ "NL",
      id == 181:190 ~ "NO",
      id == 191:200 ~ "PL",
      id == 201:210 ~ "PT",
      id == 211:220 ~ "RU",
      id == 221:230 ~ "SK",
      id == 231:240 ~ "SI",
      id == 241:250 ~ "ES",
      id == 251:260 ~ "SE",
      id == 261:270 ~ "CH",
      id == 271:280 ~ "UA",
      id == 281:290 ~ "GB"
    )
  )

# COMPUTE CI 

rmsedf_cntry <- rmsedf %>% 
  group_by(cntry) %>% 
  summarise_at(vars(lm, ranger, ksvm), funs(mean, lower, upper))

#### PLOT RMSE BY COUNTRY ####

# RANGER
ranger_plot_cntry <- rmsedf_cntry %>%
  ggplot(aes (x = ranger_mean, y = reorder(cntry, ranger_mean))) +
  geom_point() +
  geom_errorbarh(aes(xmax = ranger_upper, xmin = ranger_lower), height = 0.2) +
  ylab("") +
  xlab("") +
  geom_vline(xintercept = 1.723, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0,3, 0.5),
                     limits = c(1,2.5))  +
  ggtitle("RF") +
  theme_m()
ranger_plot_cntry

# LM PLOT
lm_plot_cntry <- rmsedf_cntry %>%
  ggplot(aes (x = lm_mean, y = reorder(cntry, lm_mean))) +
  geom_point() +
  geom_errorbarh(aes(xmax = lm_upper, xmin = lm_lower), height = 0.2) +
  ylab("") +
  xlab("") +
  geom_vline(xintercept = 1.617, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0,3, 0.5),
                     limits = c(1,2.5))  +
  ggtitle("LM") +
  theme_m(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
# lm_plot_cntry

# KSVM PLOT
ksvm_plot_cntry <- rmsedf_cntry %>%
  ggplot(aes (x = ksvm_mean, y = reorder(cntry, ksvm_mean))) +
  geom_point() +
  geom_errorbarh(aes(xmax = ksvm_upper, xmin = ksvm_lower), height = 0.2) +
  ylab("") +
  xlab("RMSE") +
  geom_vline(xintercept = 1.714, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0,3, 0.5),
                     limits = c(1,2.5))  +
  ggtitle("KSVM") +
  theme_m(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
# ksvm_plot_cntry

library(cowplot)

plot_grid(ranger_plot_cntry, ksvm_plot_cntry, lm_plot_cntry, ncol = 3)


#### FULL MODEL ####

# # Define a search space for each learner's parameter
# ps_ksvm = makeParamSet(
#   makeNumericParam("sigma", lower = -12, upper = 12, trafo = function(x) 2^x)
# )
# ps_rf = makeParamSet(
#   makeIntegerParam("num.trees", lower = 1L, upper = 200L)
# )
# # Choose a resampling strategy
# rdesc = makeResampleDesc("CV", iters = 10)
# # Choose a performance measure
# meas = rmse
# # Choose a tuning method
# ctrl = makeTuneControlCMAES(budget = 100L)
# # Make tuning wrappers
# tuned.ksvm = makeTuneWrapper(learner = "regr.ksvm", resampling = rdesc, measures = meas,
#                              par.set = ps_ksvm, control = ctrl, show.info = FALSE)
# tuned.rf = makeTuneWrapper(learner = "regr.ranger", resampling = rdesc, measures = meas,
#                            par.set = ps_rf, control = ctrl, show.info = FALSE)
# # Four learners to be compared
# lrns = list(makeLearner("regr.lm"), tuned.ksvm,
#             tuned.rf)



# Define task
regr.task_c <-  makeRegrTask(data = ESS6, target = "stfdem")

# Conduct the benchmark experiment
bmr_c <-  benchmark(learners = lrns, tasks = regr.task, 
                    resamplings = rdesc, measures = rmse, show.info = FALSE,
                    keep.pred = FALSE, models = FALSE)

# FEATURE IMPORTANCE 
lrn <- makeLearner("regr.ranger", importance = c("permutation"))
mod_c <- train(lrn, regr.task)
importance <- getFeatureImportance(mod[[i]], nmc = 50L, interaction = TRUE, local = TRUE)
# importance <- gather(importance[[i]]$res, "var.name", "importance")

imp_c <- gather(importance_c$res, "var.name", "importance")
imp_c

imp_c$var.name <- recode(imp_c$var.name,"fairelcc" = "FairElection")
imp_c$var.name <- recode(imp_c$var.name,"cttresa" = "FairCourt")
imp_c$var.name <- recode(imp_c$var.name,"fairelc" = "ElectionImp")
imp_c$var.name <- recode(imp_c$var.name,"incpart" = "IncumbentClose")
imp_c$var.name <- recode(imp_c$var.name,"GovCommunication" = "GovCom")

imp_c %>% 
  ggplot(aes(x = importance, y = reorder(var.name, importance))) +
  geom_point() +
  xlab("Importance") +
  ylab("") +
  theme_classic()


######## DESC. STATS ###########
### SUMMARY FUNCTION CI #######

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
 
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

ESS_agg <- summarySE(ESS6, measurevar = "stfdem", groupvars = "cntry")

#### Desc stat. table ####

swd_cntry <- ESS_agg %>%
  ggplot(aes (x = stfdem, y = reorder(cntry, stfdem))) +
  geom_point() +
  geom_errorbarh(aes(xmax = stfdem + (1.96 * se),
                     xmin = stfdem - (1.96 * se)), height = 0.1) +
  ylab("") +
  xlab("") +
  ggtitle("Satisfaction with democracy") +
  scale_x_continuous(
    limits = c(3, 8),
    breaks = round(seq(3, 8, .5), 2)) +
  geom_vline(xintercept = 5.419, linetype = "dashed")

swd_cntry

ESS6 %>% group_by(cntry) %>% 
  stargazer(ESS6, summary.stat = c("n", "sd", "mean", "min", "max"))
