######################################################################
# Patterns of variance  in Saponin abundance among & within leaf age #
######################################################################
# N = 84 saponins in these analyes (out of N = 86) - some removed that were very rare, not present in enough leaves to calculate variability. All N = 86 will be in the other turnover / diversity analyses

# Analyses in this script match structure of non-chem trait % change analyses (e.g. "~/Google Drive/Projects/Alfalfa project/Code/Non-chemistry traits/02.1 - sla_analyses_final.csv" and others in series)

# How does the abundance of particular compounds differ among leaves of wild vs domestic plants?
# This is the analysis of N = 86 confirmed saponins only

# For abundance analyses, difficulty of incorporating zeros. We tried a few approaches. To avoid removing these samples, Will ran left-censored models in brms to avoid removing so much data. However, it was tough because although you can estimate 0 response variables (e.g. 0 values for SD), you can't incorporate covariates as easily (e.g. the mean), which is troublesome. 
# So, it seems like we need to work with the zeros. 

# Approach: 
# a) Remove zeros
# b) For measurement of whole-plant variability, restrict analysis to compounds found at least once in each leaf age class. This is also biologically informed; if wild plants have more frequent 'absences', it is likely that we will be calculating SD for a given compound from a subset of leaves after excluding those zeros. This subset might be nonrandom in regard to leaf age class (e.g. certain leaf ages more likely to lack the compounds than others). If we want to interpret leaf age class is an important driver of the whole-plant variation (as we do for the other traits), it would be better to only consider compounds that are found across all those leaf stages. I dealt with a similar issue with the C and N data, where we sometimes lacked values for certain leaf stages (because of sample loss/accident during shipping); in that case, I excluded any plants where we didn't have a representative leaf for each leaf stage (for measures of SD at the whole-plant stage). Within leaf stages, I think this will show up as zeros (e.g. less than 2 leaves ==> can't calculate SD). The presence/absence data will be incorporated into the betadisp analyses.
# c) Calculate SD and mean from the raw data (additive variation). Inspect residuals. Best models had log(sd) or log(mean) as response, and the scaled log(mean) as a predictor.

############ Packages ################
library(ggplot2)
library(Rmisc)
library(lme4)
library(emmeans)
library(ggeffects)
library(lmerTest)
library(bbmle)
library(glmmTMB)
library(tidyr)
library(dplyr)
library(MuMIn)
library(boot)

############ Data ################
data_saps <- read.csv("~/saponins_ql_pg_final.csv")

# Other plant info (e.g. size, etc)
covariates <- read.csv("~/Plant_level_covariates.csv")

# Reminder of data processing prior to this stage
# 0.0. wrangle into same format; remove remove DDA, pooled samples, some blanks; separate "sample_id" into it's different bits
# 1.0. QL only: the reported retention time in the RT column should be within +/- 0.03 (i.e., 3.41-3.47 for at 3.44 min peak). 
# 2.0. QL only. Subtract Peak Start Time from Peak End Time and then set a minimum peak width of 0.1 min (this matches exactly what we set for the peak picking in Progenesis)
# 3.0 PG + QL. Minimum peak area: 10,000. Below this threshold, high risk of false positives. Apply to data from both methods, ensure same range. 
# 4.0 PG + QL. Normalize all peak areas to Digitoxin, in samples and in blanks
#     *Remember to exclude the EA adduct from QL data (already excluded from PG data)
#     *Average digitoxin peak areas between the two methods (QL and PG)
#     *Normalize all peak intensities, across integration methods, by this single average value of digitoxin
# 5.0. PG + QL. Remove mean amount of each compound in blanks from its sample. Do this separately for PG and QL. Save these areas separately from the original / raw peak areas. Re-Rbind the data together.

# Chemical data check
length(unique(data_saps$compound)) # 86 saponins
length(unique(data_saps$plant_pop)) # 60 plant individuals
length(unique(data_saps$leaf_stage)) # 3 leaf stages

######## Merge plant covariates with chemical data ##############
head(data_saps)
unique(data_saps$plant_pop)
data_saps2 <- merge(data_saps, covariates, by="plant_pop", all=TRUE); nrow(data_saps2); nrow(data_saps) # 46440 (if don't match, likely typo with 41_D_UIN)

######## Remove zeros ##############
data_saps_nozeros <- subset(data_saps2, peak_intensity_norm_noNoise != 0); nrow(data_saps_nozeros) # 26776 rows

# Re-level domestic & wild for downstream comparisons
data_saps_nozeros$dom_wild <- factor(data_saps_nozeros$dom_wild, levels=c("wild", "domestic"))

# Also make a column for the sqrt-transformed peak intensity data
data_saps_nozeros$peak_intensity_norm_noNoise_sqrt <- sqrt(data_saps_nozeros$peak_intensity_norm_noNoise)

##################################################################
# Next, calculate mean & SD 
# Do so at two scales: across all leaves & within leaf age class
# also: do so from raw peak intensity data and from sqrt-transformed peak intensity data
##################################################################
head(data_saps_nozeros)
#-------------#
#  Per-plant  #
#-------------#
# First, remove plants X compound combos where the compound wasn't found at least once in all leaf stages (see logic @ top of script)
head(data_saps_nozeros); data_saps_nozeros$tally <- 1
data_saps_nozeros$compound.plant_pop <- paste(data_saps_nozeros$compound, data_saps_nozeros$plant_pop, sep=".")

tally1 <- aggregate(tally ~ compound + plant_pop + leaf_stage, data=data_saps_nozeros, length); head(tally1)
tally1$tally <- 1
tally2 <- aggregate(tally ~ compound + plant_pop, data=tally1, length); head(tally2)
tally2$compound.plant_pop <- paste(tally2$compound, tally2$plant_pop, sep="."); head(tally2)
keep <- unique(subset(tally2, tally == 3))$compound.plant_pop

data_subset1 <- subset(data_saps_nozeros, compound.plant_pop %in% keep); head(data_subset1); nrow(data_subset1) # down to 24759 - not bad!

# Calculate among-all-leaf mean and sd from this subset of data (only in cases where a compound was found across all leaf ages)
# From raw peak intensities
data.pp <- do.call(data.frame, (aggregate(peak_intensity_norm_noNoise ~ compound + integration_method + plant_pop + dom_wild + plant.volume_centered + repro.stage + repro.stage.simp + percent_FAL + percent_HEM + percent_CAER, data=data_subset1, function(x) c(mean = mean(x), sd = sd(x), N = length(x))))); head(data.pp); str(data.pp); nrow(data.pp)
names(data.pp) <- c("compound", "integration_method", "plant_pop", "dom_wild", "plant.volume_centered", "repro.stage", "repro.stage.simp", "percent_FAL","percent_HEM","percent_CAER", "mean", "sd", "N"); head(data.pp)
hist(data.pp$sd); hist(data.pp$mean)
hist(log(data.pp$sd)); hist(log(data.pp$mean))

# From sqrt-transormed peak intensities
#data.pp <- do.call(data.frame, (aggregate(peak_intensity_norm_noNoise_sqrt ~ compound + integration_method + plant_pop + dom_wild + plant.volume_centered + repro.stage + repro.stage.simp, data=data_subset1, function(x) c(mean = mean(x), sd = sd(x), N = length(x))))); head(data.pp); str(data.pp); nrow(data.pp)
#names(data.pp) <- c("compound", "integration_method", "plant_pop", "dom_wild", "plant.volume_centered", "repro.stage", "repro.stage.simp", "mean", "sd", "N"); head(data.pp)
#hist(data.pp$sd); hist(data.pp$mean)

# for each compound, this dataframe has its mean and sd per plant (= among 9 leaves), only for plants in which the compound was found among all leaf ages (e.g. in at least 1 young, middle, and old leaf).

length(unique(data.pp$compound))  # 84 - this removed a couple compounds! Interesting!! Makes sense - for some very rare compounds, SD of peak intensity won't make as much sense to calculate. Presence/absence makes more sense.
head(data.pp); str(data.pp)

#-------------------------------#
#  Per leaf relative age class  #
#-------------------------------#
# As above, but this time remove any plant x leaf stage x compound combos with less than 2 reps
data_saps_nozeros$compound.plant_pop.leaf_stage <- paste(data_saps_nozeros$compound.plant_pop, data_saps_nozeros$leaf_stage, sep="."); 
tally1 <- aggregate(tally ~ compound.plant_pop.leaf_stage, data=data_saps_nozeros, length); head(tally1)
keep2 <- unique(subset(tally1, tally != 1))$compound.plant_pop.leaf_stage

data_subset2 <- subset(data_saps_nozeros, compound.plant_pop.leaf_stage %in% keep2); head(data_subset2); nrow(data_subset2) # down to 25939 - not bad at all!

# across all leaves of all relative age classes); keep relevant plant-level covariates
# With raw peak intensities
data.ps <- do.call(data.frame, (aggregate(peak_intensity_norm_noNoise ~ compound + integration_method + plant_pop + dom_wild + leaf_stage + plant.volume_centered + repro.stage + repro.stage.simp + percent_FAL + percent_HEM + percent_CAER, data=data_subset2, function(x) c(mean = mean(x), sd = sd(x), N = length(x))))); head(data.ps); str(data.ps); nrow(data.ps)
names(data.ps) <- c("compound", "integration_method", "plant_pop", "dom_wild", "leaf_stage", "plant.volume_centered", "repro.stage", "repro.stage.simp", "percent_FAL","percent_HEM","percent_CAER", "mean", "sd", "N"); head(data.ps)
data.ps$leaf_stage <- factor(data.ps$leaf_stage, levels=c("young", "middle", "old"))

# With sqrt-transformed peak intensities
#data.ps <- do.call(data.frame, (aggregate(peak_intensity_norm_noNoise_sqrt ~ compound + integration_method + plant_pop + dom_wild + leaf_stage + plant.volume_centered + repro.stage + repro.stage.simp, data=data_subset2, function(x) c(mean = mean(x), sd = sd(x), N = length(x))))); head(data.ps); str(data.ps); nrow(data.ps)
#names(data.ps) <- c("compound", "integration_method", "plant_pop", "dom_wild", "leaf_stage", "plant.volume_centered", "repro.stage", "repro.stage.simp", "mean", "sd", "N"); head(data.ps)
#unique(data.ps$N) # good, only 2 or 3 leaves (no singletons)

length(unique(data.ps$compound)) # also 84 saponins
head(data.ps); str(data.ps)

#####################################################
# 1) Model fitting: whole-plant / across all leaves #
#####################################################
head(data.pp); str(data.pp)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# First, fit models with both mean and domestication status #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# What is the change in variability with domestication, accounting for effects of the mean?
data.pp$sd.log <- log(data.pp$sd); head(data.pp) # bootMer needs pre-transformed response var
# log the mean (improved model residuals a lot vs runnning the raw mean)
data.pp$mean.log <- log(data.pp$mean)
# scale (mean-center) the mean
data.pp$mean.log_centered <- data.pp$mean.log - mean(data.pp$mean.log)

##### A) Tweak ranFX structure; allow the mean~var relationship to vary within domestication status, within compound, within plant individual
m1 <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + (1 + mean.log_centered|dom_wild) + (1 + mean|compound) + (1 + mean.log_centered|plant_pop), data=data.pp) # failed to converge
m2 <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp) # converged OK
m3 <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + (1|repro.stage) + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp) # converged OK
m4 <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + (1|repro.stage) + (1 + mean.log_centered|dom_wild) + (1 + mean.log_centered|compound), data=data.pp) # failed to converge
m5 <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + (1|repro.stage) + (1 + mean.log_centered|dom_wild) + (1 + mean.log_centered|plant_pop), data=data.pp) #  failed to converge
m6 <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + (1|compound) + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp) # failed to converge
AICctab(m2, m3) # m2 best; mean~variance relationships differ within compound and within plant individual. 

# Quick-check if adding reproductive stage is important
data.pp$repro.stage_dummy <- as.factor(sample(unique(data.pp$repro.stage), replace=T, nrow(data.pp))); head(data.pp)
m2_repro <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + (1|repro.stage) + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp) # converged OK
m2_dummy <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + (1|repro.stage_dummy) + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp) # converged OK
AICctab(m2_repro, m2_dummy) # the same - so, having reproductive stage doesn't add anything to the model

##### B) Tweak fixed FX structure of m2
# First, quick-check if adding integration method is important
m2_noInt <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp) 
AICctab(m2, m2_noInt) # still good to have integration method in there; makes sense, as QuanLynx was focused on lower-abundance compounds. I think I could justify models with or without that fixed effect
# Dredge m2
m2.1 <- glmmTMB(log(sd) ~ dom_wild*mean.log_centered*integration_method + plant.volume_centered + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp) 
m2.2 <- glmmTMB(log(sd) ~ dom_wild*mean.log_centered + integration_method + plant.volume_centered + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp) 
m2.3 <- glmmTMB(log(sd) ~ dom_wild*mean.log_centered + integration_method + plant.volume_centered + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp) 
m2.4 <- glmmTMB(log(sd) ~ dom_wild*integration_method + mean.log_centered + plant.volume_centered + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp) 
AICctab(m2, m2.1, m2.2, m2.3, m2.4) # no interactions - go with m2

# See if adding percent wild ssp improves model
m2_both <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + percent_FAL + percent_CAER + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp)
m2_caer <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + percent_CAER + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp)
m2_fal <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + percent_FAL + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp)
AICctab(m2, m2_both, m2_caer, m2_fal) # m2 still best

# Observation-level random effect
data.pp$observation <- paste(data.pp$compound, data.pp$plant_pop, sep="-")
m2.obs <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop) + + (1|observation), data=data.pp)
AICctab(m2, m2.obs) # m2 best

# Name the top models
topmodel.pp <- glmmTMB(log(sd) ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp)
# Version for bootMer
topmodel.pp2 <- glmmTMB(sd.log ~ dom_wild + mean.log_centered + plant.volume_centered + integration_method + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.pp)

DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.pp)) # hmm -  significant deviations.. maybe just because I have SO much power here? Significance comes from many tiny deviations... so, yes model isn't perfect, but doesn't look like it's performing terribly. Modelling a lot of data, and a lot of power in the random effects. How are residuals weighted (with this much data, could be spurious p-value). Probably just taking random effects into account, but not using them to calculate deviations between predicted and observed - that's a lot of data points and probably many small deviations. Shapiro-Wilkes test often comes down to sample size -- so, hard to interpret. How does the model do with prediction?
hist(resid(topmodel.pp))
plot(residuals(topmodel.pp) ~ predict(topmodel.pp))
plot(residuals(topmodel.pp) ~ data.pp$mean.log_centered); lines(lowess(resid(topmodel.pp) ~ data.pp$mean.log_centered), col="red") # this looks ok, pretty close to zero and flat

performance::icc(topmodel.pp)

# Best to predict the data, and see how the model does.
obs.sd.log <- log(data.pp$sd)
pred.sd_log <- predict(topmodel.pp)
data <- data.frame(obs.sd_log, pred.sd_log)
ggplot(data, aes(x=obs.sd_log, y=pred.sd_log)) +
  geom_point(alpha=.5)+
  geom_smooth()+
  geom_abline(slope=1, intercept=0, col="red")  # overpredicting variability (SD) for low-var compounds, but generally OK

##### C) Pull out predicted means for wild & domestic + model nobs + df
model.ests <- data.frame(ggpredict(topmodel.pp, terms="dom_wild"))
ggplot(model.ests, aes(x=x, y=predicted)) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), size=.5, width=.2, color="black")

wild.pred <- model.ests$predicted[which(model.ests$x == "wild")]; wild.pred
wild.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "wild")]; wild.pred_conf.low
wild.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "wild")]; wild.pred_conf.high

dom.pred <- model.ests$predicted[which(model.ests$x == "domestic")]; dom.pred
dom.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "domestic")]; dom.pred_conf.low
dom.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "domestic")]; dom.pred_conf.high

model_nobs <- nobs(topmodel.pp)
model_df.residual <- df.residual(topmodel.pp)

ests <- data.frame(model_nobs, model_df.residual, wild.pred, wild.pred_conf.low, wild.pred_conf.high, dom.pred, dom.pred_conf.low, dom.pred_conf.high); ests

##### D) Calculate / append significance  
sig1 <- data.frame(car::Anova(topmodel.pp, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- cbind(ests, sig)

##### E) Bootstrap CI around beta
# This is from the March 15, 2020 glmmTMB vignette. See "isLMM.glmmTMB" section. Instead of "$zi", in their example, I am interested in the cond estimate for domestic - so I name that specifically within the function
fixef(topmodel.pp)
fixef(topmodel.pp2)
fixef(topmodel.pp2)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object
b1 <- lme4::bootMer(topmodel.pp2, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")
boot <- boot.ci(b1,type="perc")

# Pull out the effect (beta) from the model
beta <- as.numeric(fixef(topmodel.pp)$cond["dom_wilddomestic"])
# Pull out the bootstrapped 95% CI around that beta. I need to get this from the "percent" part of the bootstrapping output. The items within aren't named, which is annoying. From bootMer, below "These latter four components will be matrices with 5 columns, the first column containing the level, the next two containing the indices of the order statistics used in the calculations and the final two the calculated endpoints themselves." So I need to pull the latter two.
boot; str(boot); class(boot)
out <- data.frame(boot$percent); out # so the last two components are the upper and lower CIs, respecively
beta.CI_upper.boot <- as.numeric(out[length(out)]); beta.CI_upper.boot 
beta.CI_lower.boot <- as.numeric(out[length(out)-1]); beta.CI_lower.boot
# Put them together
ci <- data.frame(beta, beta.CI_lower.boot, beta.CI_upper.boot)

##### F) Calculate % changes
# This is a little different, because of the log-transformation: See here
# https://stats.stackexchange.com/questions/362556/linear-mixed-effect-model-interpretation-with-log-transformed-dependent-variable; also https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
# "In addition, you should interpret your transformed beta coefficients (and associated confidence intervals) as multiplicative changes once you transform to the original units. E.g., beta from the model = 0.0085; exp(.008585) ~= 1.0086 = 0.86% change"... "Only the dependent/response variable is log-transformed. Exponentiate the coefficient, subtract one from this number, and multiply by 100. This gives the percent increase (or decrease) in the response for every one-unit increase in the independent variable. Example: the coefficient is 0.198. (exp(0.198) – 1) * 100 = 21.9. For every one-unit increase in the independent variable, our dependent variable increases by about 22%."
change.pp <- cbind(pred_sig, ci); change.pp
change.pp$exp.beta <- exp(change.pp$beta)
change.pp$exp.beta.lower <- exp(change.pp$beta.CI_lower.boot)
change.pp$exp.beta.upper <- exp(change.pp$beta.CI_upper.boot)
change.pp
change.pp$percent.change <- (change.pp$exp.beta-1)*100
change.pp$percent.change.lower <-(change.pp$exp.beta.lower-1)*100
change.pp$percent.change.upper <- (change.pp$exp.beta.upper-1)*100

change.pp$leaf_stage <- "all"
change.pp$response <- "sd"
change.pp$scale <- "per-plant"
change.pp$model <- "status+mean"

##### G) Put it all together & re-order for merging
change.pp <- change.pp[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "wild.pred_conf.low", "wild.pred_conf.high", "dom.pred", "dom.pred_conf.low", "dom.pred_conf.high", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.pp

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Next, fit models with just domestication status # 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# What is the change in variability with domestication, without accounting for effects of the mean?
##### A) Tweak ranFX structure NA
##### B) Tweak fixed FX structure NA
# So this is the log of data, with the original peak intensities sqrt-transformed first
topmodel.pp_nomean <- glmmTMB(log(sd) ~ dom_wild + plant.volume_centered + integration_method + (1|compound) + (1|plant_pop), data=data.pp)
topmodel.pp_nomean2 <- glmmTMB(sd.log ~ dom_wild + plant.volume_centered + integration_method + (1|compound) + (1|plant_pop), data=data.pp) # for bootMer
hist(resid(topmodel.pp_nomean))
plot(resid(topmodel.pp_nomean))

##### C) Pull out predicted means for wild & domestic + model nobs + df
model.ests <- data.frame(ggpredict(topmodel.pp_nomean, terms="dom_wild"))
ggplot(model.ests, aes(x=x, y=predicted)) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), size=.5, width=.2, color="black")

wild.pred <- model.ests$predicted[which(model.ests$x == "wild")]; wild.pred
wild.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "wild")]; wild.pred_conf.low
wild.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "wild")]; wild.pred_conf.high

dom.pred <- model.ests$predicted[which(model.ests$x == "domestic")]; dom.pred
dom.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "domestic")]; dom.pred_conf.low
dom.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "domestic")]; dom.pred_conf.high

model_nobs <- nobs(topmodel.pp_nomean)
model_df.residual <- df.residual(topmodel.pp_nomean)

ests <- data.frame(model_nobs, model_df.residual, wild.pred, wild.pred_conf.low, wild.pred_conf.high, dom.pred, dom.pred_conf.low, dom.pred_conf.high); ests

##### D) Calculate / append significance  
sig1 <- data.frame(car::Anova(topmodel.pp_nomean, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- cbind(ests, sig)
pred_sig

##### E) Bootstrap CI 
# This is from the March 15, 2020 glmmTMB vignette. See "isLMM.glmmTMB" section. Instead of "$zi", in their example, I am interested in the cond estimate for domestic - so I name that specifically within the function
fixef(topmodel.pp_nomean); fixef(topmodel.pp_nomean2)
fixef(topmodel.pp_nomean2)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object
b1 <- lme4::bootMer(topmodel.pp_nomean2, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")
boot <- boot.ci(b1,type="perc")

# Pull out the effect (beta) from the model
beta <- as.numeric(fixef(topmodel.pp_nomean)$cond["dom_wilddomestic"])
# Pull out the bootstrapped 95% CI around that beta. I need to get this from the "percent" part of the bootstrapping output. The items within aren't named, which is annoying. From bootMer, below "These latter four components will be matrices with 5 columns, the first column containing the level, the next two containing the indices of the order statistics used in the calculations and the final two the calculated endpoints themselves." So I need to pull the latter two.
boot; str(boot); class(boot)
out <- data.frame(boot$percent); out # so the last two components are the upper and lower CIs, respecively
beta.CI_upper.boot <- as.numeric(out[length(out)]); beta.CI_upper.boot
beta.CI_lower.boot <- as.numeric(out[length(out)-1]); beta.CI_lower.boot
# Put them together
ci <- data.frame(beta, beta.CI_lower.boot, beta.CI_upper.boot)

##### F) Calculate % changes
# This is a little different, because of the log-transformation: See here
# https://stats.stackexchange.com/questions/362556/linear-mixed-effect-model-interpretation-with-log-transformed-dependent-variable; also https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
# "In addition, you should interpret your transformed beta coefficients (and associated confidence intervals) as multiplicative changes once you transform to the original units. E.g., beta from the model = 0.0085; exp(.008585) ~= 1.0086 = 0.86% change"... "Only the dependent/response variable is log-transformed. Exponentiate the coefficient, subtract one from this number, and multiply by 100. This gives the percent increase (or decrease) in the response for every one-unit increase in the independent variable. Example: the coefficient is 0.198. (exp(0.198) – 1) * 100 = 21.9. For every one-unit increase in the independent variable, our dependent variable increases by about 22%."
change.pp_nomean <- cbind(pred_sig, ci); change.pp_nomean
change.pp_nomean$exp.beta <- exp(change.pp_nomean$beta)
change.pp_nomean$exp.beta.lower <- exp(change.pp_nomean$beta.CI_lower.boot)
change.pp_nomean$exp.beta.upper <- exp(change.pp_nomean$beta.CI_upper.boot)
change.pp_nomean
change.pp_nomean$percent.change <- (change.pp_nomean$exp.beta-1)*100
change.pp_nomean$percent.change.lower <-(change.pp_nomean$exp.beta.lower-1)*100
change.pp_nomean$percent.change.upper <- (change.pp_nomean$exp.beta.upper-1)*100

change.pp_nomean$leaf_stage <- "all"
change.pp_nomean$response <- "sd"
change.pp_nomean$scale <- "per-plant"
change.pp_nomean$model <- "status_only"

##### G) Put it all together & re-order for merging
change.pp_nomean <- change.pp_nomean[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "wild.pred_conf.low", "wild.pred_conf.high", "dom.pred", "dom.pred_conf.low", "dom.pred_conf.high", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.pp_nomean

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#--#-#-#-#-#-#
# Next, fit models to show shift in mean #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# What is the change in mean with domestication?
head(data.pp)

##### A) Tweak ranFX structure; allow the mean~var relationship to vary within domestication status, within compound, within plant individual
m1 <- glmmTMB(log(mean) ~ dom_wild + plant.volume_centered + integration_method + (1|repro.stage) + (1|compound) + (1|plant_pop), data=data.pp) 
m2 <- glmmTMB(log(mean) ~ dom_wild + plant.volume_centered + integration_method + (1|compound) + (1|plant_pop), data=data.pp) 
m3 <- glmmTMB(log(mean) ~ dom_wild + plant.volume_centered + integration_method + (1|repro.stage) + (1|plant_pop), data=data.pp) 
m4 <- glmmTMB(log(mean) ~ dom_wild + plant.volume_centered + integration_method + (1|compound), data=data.pp) 
m5 <- glmmTMB(log(mean) ~ dom_wild + plant.volume_centered + integration_method + (1|plant_pop), data=data.pp) 
AICctab(m1, m2, m3, m4, m5) # m2 best; compound + plant individual

##### B) Tweak fixed FX structure of m2
# First, quick-check if adding integration method is important
m2_noInt <- glmmTMB(log(mean) ~ dom_wild + plant.volume_centered + (1|compound) + (1|plant_pop), data=data.pp) 
AICctab(m2, m2_noInt) # hmm not much difference, actually.
dredge(glmmTMB(mean.log ~ dom_wild*plant.volume_centered + integration_method + (1|compound) + (1|plant_pop), data=data.pp))
# Within 2 AIC, all terms, no interactions

# Check if ancestry / wild subspecies % improves model
m2_both <- glmmTMB(log(mean) ~ dom_wild + plant.volume_centered + integration_method + percent_CAER + percent_FAL + (1|compound) + (1|plant_pop), data=data.pp)
m2_caer <- glmmTMB(log(mean) ~ dom_wild + plant.volume_centered + integration_method + percent_CAER + (1|compound) + (1|plant_pop), data=data.pp)
m2_fal <- glmmTMB(log(mean) ~ dom_wild + plant.volume_centered + integration_method + percent_FAL + (1|compound) + (1|plant_pop), data=data.pp)
AICctab(m2, m2_both,m2_caer, m2_fal) # m2 best

topmodel.mean <- glmmTMB(log(mean) ~ dom_wild + plant.volume_centered + integration_method + (1|compound) + (1|plant_pop), data=data.pp)
# version for bootMer
topmodel.mean2 <- glmmTMB(mean.log ~ dom_wild + plant.volume_centered + integration_method + (1|compound) + (1|plant_pop), data=data.pp)

plot(residuals(topmodel.mean) ~ predict(topmodel.mean))
hist(residuals(topmodel.mean))

##### C) Pull out predicted means for wild & domestic + model nobs + df
model.ests <- data.frame(ggpredict(topmodel.mean, terms="dom_wild"))
ggplot(model.ests, aes(x=x, y=predicted)) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), size=.5, width=.2, color="black")

wild.pred <- model.ests$predicted[which(model.ests$x == "wild")]; wild.pred
wild.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "wild")]; wild.pred_conf.low
wild.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "wild")]; wild.pred_conf.high

dom.pred <- model.ests$predicted[which(model.ests$x == "domestic")]; dom.pred
dom.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "domestic")]; dom.pred_conf.low
dom.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "domestic")]; dom.pred_conf.high

model_nobs <- nobs(topmodel.mean)
model_df.residual <- df.residual(topmodel.mean)

ests <- data.frame(model_nobs, model_df.residual, wild.pred, wild.pred_conf.low, wild.pred_conf.high, dom.pred, dom.pred_conf.low, dom.pred_conf.high); ests

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.mean, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- cbind(ests, sig)

##### E) Bootstrap CI 
# This is from the March 15, 2020 glmmTMB vignette. See "isLMM.glmmTMB" section. Instead of "$zi", in their example, I am interested in the cond estimate for domestic - so I name that specifically within the function
fixef(topmodel.mean); fixef(topmodel.mean)
fixef(topmodel.mean2)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object
b1 <- lme4::bootMer(topmodel.mean2, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")
boot <- boot.ci(b1,type="perc")

# Pull out the effect (beta) from the model
beta <- as.numeric(fixef(topmodel.mean)$cond["dom_wilddomestic"])
# Pull out the bootstrapped 95% CI around that beta. I need to get this from the "percent" part of the bootstrapping output. The items within aren't named, which is annoying. From bootMer, below "These latter four components will be matrices with 5 columns, the first column containing the level, the next two containing the indices of the order statistics used in the calculations and the final two the calculated endpoints themselves." So I need to pull the latter two.
boot; str(boot); class(boot)
out <- data.frame(boot$percent); out # so the last two components are the upper and lower CIs, respecively
beta.CI_upper.boot <- as.numeric(out[length(out)]); beta.CI_upper.boot
beta.CI_lower.boot <- as.numeric(out[length(out)-1]); beta.CI_lower.boot
# Put them together
ci <- data.frame(beta, beta.CI_lower.boot, beta.CI_upper.boot)

##### F) Calculate % changes
# This is a little different, because of the log-transformation: See here
# https://stats.stackexchange.com/questions/362556/linear-mixed-effect-model-interpretation-with-log-transformed-dependent-variable; also https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
# "In addition, you should interpret your transformed beta coefficients (and associated confidence intervals) as multiplicative changes once you transform to the original units. E.g., beta from the model = 0.0085; exp(.008585) ~= 1.0086 = 0.86% change"... "Only the dependent/response variable is log-transformed. Exponentiate the coefficient, subtract one from this number, and multiply by 100. This gives the percent increase (or decrease) in the response for every one-unit increase in the independent variable. Example: the coefficient is 0.198. (exp(0.198) – 1) * 100 = 21.9. For every one-unit increase in the independent variable, our dependent variable increases by about 22%."
change.mean <- cbind(pred_sig, ci); change.mean
change.mean$exp.beta <- exp(change.mean$beta)
change.mean$exp.beta.lower <- exp(change.mean$beta.CI_lower.boot)
change.mean$exp.beta.upper <- exp(change.mean$beta.CI_upper.boot)
change.mean
change.mean$percent.change <- (change.mean$exp.beta-1)*100
change.mean$percent.change.lower <-(change.mean$exp.beta.lower-1)*100
change.mean$percent.change.upper <- (change.mean$exp.beta.upper-1)*100

change.mean$leaf_stage <- "all"
change.mean$response <- "mean"
change.mean$scale <- "per-plant"
change.mean$model <- "status_only"

##### G) Put it all together & re-order for merging
change.mean <- change.mean[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "wild.pred_conf.low", "wild.pred_conf.high", "dom.pred", "dom.pred_conf.low", "dom.pred_conf.high", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.mean

#--------------------------------------------------------------#
# Rbind all % changes for the whole-plant / across all leaves  #
#--------------------------------------------------------------#
estimates.pp <- rbind(change.pp, change.pp_nomean, change.mean); estimates.pp

#######################################
# 2) Model fitting: by leaf age class #
#######################################
head(data.ps); str(data.ps)

# What is the change in variability with domestication, accounting for effects of the mean?
data.ps$sd.log <- log(data.ps$sd); head(data.ps) # bootMer needs pre-transformed response var
# log the mean (improved model residuals a lot)
data.ps$mean.log <- log(data.ps$mean)
# scale (mean-center) the mean (following meeting w/ Will)
data.ps$mean.log_centered <- data.ps$mean.log - mean(data.ps$mean.log)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# First, fit models with both mean and domestication status #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# What is the change in variability within each leaf age class with domestication, accounting for effects of the mean?

##### A) Tweak ranFX structure
m1 <- glmmTMB(log(sd) ~ dom_wild + leaf_stage + mean.log_centered + integration_method + (1 + mean.log_centered|dom_wild) + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.ps)  # failed to converge
m2 <- glmmTMB(log(sd) ~ dom_wild + leaf_stage + mean.log_centered + integration_method + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.ps) # converged
m3 <- glmmTMB(log(sd) ~ dom_wild + leaf_stage + mean.log_centered + integration_method + (1 + mean.log_centered|dom_wild) + (1 + mean.log_centered|compound), data=data.ps) # failed to converge
m4 <- glmmTMB(log(sd) ~ dom_wild + leaf_stage + mean.log_centered + integration_method + (1 + mean.log_centered|dom_wild) + (1 + mean.log_centered|plant_pop), data=data.ps) # failed to converge
m5 <- glmmTMB(log(sd) ~ dom_wild + leaf_stage + mean.log_centered + integration_method + (1 + mean.log_centered|leaf_stage) + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.ps) # converged
AICctab(m2, m5) # Go with m2

# Manipulate fixed-effect interaction structure
m2.1 <- glmmTMB(sd.log ~ dom_wild*leaf_stage*mean.log_centered + integration_method + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.ps)
m2.2 <- glmmTMB(sd.log ~ dom_wild*leaf_stage + mean.log_centered + integration_method + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.ps)
m2.3 <- glmmTMB(sd.log ~ dom_wild*mean.log_centered + leaf_stage + integration_method + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.ps)
m2.4 <- glmmTMB(sd.log ~ dom_wild + leaf_stage + mean.log_centered + integration_method + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.ps)
AICctab(m2.1, m2.2, m2.3, m2.4) # m2.2 better, with the 2-way interaction
# Confirm that integration method is important
m2.5 <- glmmTMB(sd.log ~ dom_wild*leaf_stage + mean.log_centered + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.ps)
AICctab(m2.2, m2.5) # Yep, m2.2 better

#Check for role of wild ssp ancestry
m2.2_both <- glmmTMB(sd.log ~ dom_wild*leaf_stage + mean.log_centered + integration_method + percent_CAER + percent_FAL + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.ps)
m2.2_caer <- glmmTMB(sd.log ~ dom_wild*leaf_stage + mean.log_centered + integration_method + percent_CAER + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.ps)
m2.2_fal <- glmmTMB(sd.log ~ dom_wild*leaf_stage + mean.log_centered + integration_method + percent_FAL + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.ps)
AICctab(m2.2, m2.2_both, m2.2_caer, m2.2_fal) # Yep including FAL best

summary(m2.2)
DHARMa::testResiduals(DHARMa::simulateResiduals(m2.2))
plot(residuals(m2.2) ~ predict(m2.2))
lines(lowess(resid(m2.2) ~ data.ps$mean.log_centered), col="red") # this looks ok, pretty close to zero and flat
hist(resid(m2.2))

topmodel.ps <- glmmTMB(log(sd) ~ dom_wild*leaf_stage + mean.log_centered + integration_method + percent_FAL + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.ps)
# For bootMer
topmodel.ps2 <- glmmTMB(sd.log ~ dom_wild*leaf_stage + mean.log_centered + integration_method + percent_FAL + (1 + mean.log_centered|compound) + (1 + mean.log_centered|plant_pop), data=data.ps) 

##### C) Pull out predicted means for wild & domestic + model nobs + df
emmeans(topmodel.ps, ~ dom_wild|leaf_stage, type="response")
model.ests <- data.frame(ggpredict(topmodel.ps, terms=c("dom_wild", "leaf_stage")))
ggplot(model.ests, aes(x=x, y=predicted)) +
  facet_grid(~group) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), size=.5, width=.2, color="black")

##---- Young ----##
wild.pred <- model.ests$predicted[which(model.ests$x == "wild" & model.ests$group == "young")]; wild.pred
wild.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "wild" & model.ests$group == "young")]; wild.pred_conf.low
wild.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "wild" & model.ests$group == "young")]; wild.pred_conf.high
dom.pred <- model.ests$predicted[which(model.ests$x == "domestic" & model.ests$group == "young")]; dom.pred
dom.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "domestic" & model.ests$group == "young")]; dom.pred_conf.low
dom.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "domestic" & model.ests$group == "young")]; dom.pred_conf.high
young <- data.frame(wild.pred, wild.pred_conf.low, wild.pred_conf.high, dom.pred, dom.pred_conf.low, dom.pred_conf.high)
young$leaf_stage <- "young"

##---- Middle ----##
wild.pred <- model.ests$predicted[which(model.ests$x == "wild" & model.ests$group == "middle")]; wild.pred
wild.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "wild" & model.ests$group == "middle")]; wild.pred_conf.low
wild.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "wild" & model.ests$group == "middle")]; wild.pred_conf.high
dom.pred <- model.ests$predicted[which(model.ests$x == "domestic" & model.ests$group == "middle")]; dom.pred
dom.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "domestic" & model.ests$group == "middle")]; dom.pred_conf.low
dom.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "domestic" & model.ests$group == "middle")]; dom.pred_conf.high
middle <- data.frame(wild.pred, wild.pred_conf.low, wild.pred_conf.high, dom.pred, dom.pred_conf.low, dom.pred_conf.high)
middle$leaf_stage <- "middle"

##---- Old ----##
wild.pred <- model.ests$predicted[which(model.ests$x == "wild" & model.ests$group == "old")]; wild.pred
wild.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "wild" & model.ests$group == "old")]; wild.pred_conf.low
wild.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "wild" & model.ests$group == "old")]; wild.pred_conf.high
dom.pred <- model.ests$predicted[which(model.ests$x == "domestic" & model.ests$group == "old")]; dom.pred
dom.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "domestic" & model.ests$group == "old")]; dom.pred_conf.low
dom.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "domestic" & model.ests$group == "old")]; dom.pred_conf.high
old <- data.frame(wild.pred, wild.pred_conf.low, wild.pred_conf.high, dom.pred, dom.pred_conf.low, dom.pred_conf.high)
old$leaf_stage <- "old"

ests <- rbind(young,middle,old); ests
# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.ps)
ests$model_df.residual <- df.residual(topmodel.ps)

##### D) Calculate / append significance 
sig0 <- data.frame(pairs(emmeans(topmodel.ps, ~ dom_wild|leaf_stage))) # post-hoc test 
names(sig0)[names(sig0) == 't.ratio'] <- 'statistic'
names(sig0)[names(sig0) == 'df'] <- 'statistic.df'
sig0$test.statistic <- "t.ratio"
sig <- subset(sig0, select=c("leaf_stage", "test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- merge(ests, sig, by="leaf_stage")

##### E) Bootstrap CI 
# Here, I'll extract the effect of domestication within each leaf stage, due to the significant interaction 
# (Used this great post as a resource: https://www.r-bloggers.com/bootstrapping-follow-up-contrasts-for-within-subject-anovas-part-2/)
# Also used this helpful page to re-level the order of comparisons in emmeans (for some reason it always considers domestic the 'control' group, so estimates are negative rather than positive). 
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
# For example - 
normal <- pairs(emmeans(topmodel.ps2, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.ps2, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.ps2, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
summ <- data.frame(summary(b1)); summ
names(summ)[names(summ) == 'original'] <- 'beta'
summ$leaf_stage <- rownames(summ)
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
conf$leaf_stage <- rownames(conf)
ci0 <- merge(summ, conf, by="leaf_stage")
ci <- subset(ci0, select=c("leaf_stage", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot")); ci

##### F) Calculate % changes
# This is a little different, because of the log-transformation: See here
# https://stats.stackexchange.com/questions/362556/linear-mixed-effect-model-interpretation-with-log-transformed-dependent-variable; also https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
# "In addition, you should interpret your transformed beta coefficients (and associated confidence intervals) as multiplicative changes once you transform to the original units. E.g., beta from the model = 0.0085; exp(.008585) ~= 1.0086 = 0.86% change"... "Only the dependent/response variable is log-transformed. Exponentiate the coefficient, subtract one from this number, and multiply by 100. This gives the percent increase (or decrease) in the response for every one-unit increase in the independent variable. Example: the coefficient is 0.198. (exp(0.198) – 1) * 100 = 21.9. For every one-unit increase in the independent variable, our dependent variable increases by about 22%."
change.ps <- merge(pred_sig, ci, by="leaf_stage"); change.ps
change.ps$exp.beta <- exp(change.ps$beta)
change.ps$exp.beta.lower <- exp(change.ps$beta.CI_lower.boot)
change.ps$exp.beta.upper <- exp(change.ps$beta.CI_upper.boot)
change.ps
change.ps$percent.change <- (change.ps$exp.beta-1)*100
change.ps$percent.change.lower <-(change.ps$exp.beta.lower-1)*100
change.ps$percent.change.upper <- (change.ps$exp.beta.upper-1)*100

change.ps$response <- "sd"
change.ps$scale <- "per-age-class"
change.ps$model <- "status+mean"

##### G) Put it all together & re-order for merging
change.ps <- change.ps[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "wild.pred_conf.low", "wild.pred_conf.high", "dom.pred", "dom.pred_conf.low", "dom.pred_conf.high", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.ps

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Next, fit models with just domestication status #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# What is the change in variability within each leaf age class with domestication, without accounting for the mean?
##### A) Tweak ranFX structure NA
##### B) Tweak fixed FX structure
m1 <- glmmTMB(log(sd) ~ dom_wild*leaf_stage + integration_method + (1|compound) + (1|plant_pop), data=data.ps)
m2 <- glmmTMB(log(sd) ~ dom_wild + leaf_stage + integration_method + (1|compound) + (1|plant_pop), data=data.ps)
AICctab(m1, m2)

topmodel.ps_nomean <- glmmTMB(log(mean) ~ dom_wild*leaf_stage + integration_method + percent_FAL + (1|compound) + (1|plant_pop), data=data.ps)
# For bootMer
topmodel.ps2_nomean <- glmmTMB(mean.log ~ dom_wild*leaf_stage + integration_method + percent_FAL + (1|compound) + (1|plant_pop), data=data.ps)
hist(resid(topmodel.ps_nomean))
plot(residuals(topmodel.ps_nomean) ~ predict(topmodel.ps_nomean))

##### C) Pull out predicted means for wild & domestic + model nobs + df
model.ests <- data.frame(ggpredict(topmodel.ps_nomean, terms=c("dom_wild", "leaf_stage")))
ggplot(model.ests, aes(x=x, y=predicted)) +
  facet_grid(~group) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), size=.5, width=.2, color="black")

##---- Young ----##
wild.pred <- model.ests$predicted[which(model.ests$x == "wild" & model.ests$group == "young")]; wild.pred
wild.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "wild" & model.ests$group == "young")]; wild.pred_conf.low
wild.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "wild" & model.ests$group == "young")]; wild.pred_conf.high
dom.pred <- model.ests$predicted[which(model.ests$x == "domestic" & model.ests$group == "young")]; dom.pred
dom.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "domestic" & model.ests$group == "young")]; dom.pred_conf.low
dom.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "domestic" & model.ests$group == "young")]; dom.pred_conf.high
young <- data.frame(wild.pred, wild.pred_conf.low, wild.pred_conf.high, dom.pred, dom.pred_conf.low, dom.pred_conf.high)
young$leaf_stage <- "young"

##---- Middle ----##
wild.pred <- model.ests$predicted[which(model.ests$x == "wild" & model.ests$group == "middle")]; wild.pred
wild.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "wild" & model.ests$group == "middle")]; wild.pred_conf.low
wild.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "wild" & model.ests$group == "middle")]; wild.pred_conf.high
dom.pred <- model.ests$predicted[which(model.ests$x == "domestic" & model.ests$group == "middle")]; dom.pred
dom.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "domestic" & model.ests$group == "middle")]; dom.pred_conf.low
dom.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "domestic" & model.ests$group == "middle")]; dom.pred_conf.high
middle <- data.frame(wild.pred, wild.pred_conf.low, wild.pred_conf.high, dom.pred, dom.pred_conf.low, dom.pred_conf.high)
middle$leaf_stage <- "middle"

##---- Old ----##
wild.pred <- model.ests$predicted[which(model.ests$x == "wild" & model.ests$group == "old")]; wild.pred
wild.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "wild" & model.ests$group == "old")]; wild.pred_conf.low
wild.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "wild" & model.ests$group == "old")]; wild.pred_conf.high
dom.pred <- model.ests$predicted[which(model.ests$x == "domestic" & model.ests$group == "old")]; dom.pred
dom.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "domestic" & model.ests$group == "old")]; dom.pred_conf.low
dom.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "domestic" & model.ests$group == "old")]; dom.pred_conf.high
old <- data.frame(wild.pred, wild.pred_conf.low, wild.pred_conf.high, dom.pred, dom.pred_conf.low, dom.pred_conf.high)
old$leaf_stage <- "old"

ests <- rbind(young,middle,old); ests
# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.ps_nomean)
ests$model_df.residual <- df.residual(topmodel.ps_nomean)

##### D) Calculate / append significance 
sig0 <- data.frame(pairs(emmeans(topmodel.ps_nomean, ~ dom_wild|leaf_stage))) # post-hoc test 
names(sig0)[names(sig0) == 't.ratio'] <- 'statistic'
names(sig0)[names(sig0) == 'df'] <- 'statistic.df'
sig0$test.statistic <- "t.ratio"
sig <- subset(sig0, select=c("leaf_stage", "test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- merge(ests, sig, by="leaf_stage")

##### E) Bootstrap CI 
# Here, instead of extracting the single fixed effect for domestication across all leaf stages, I need to extract the leaf stage-specific effects due to the significant interaction 
# (Used this great post as a resource: https://www.r-bloggers.com/bootstrapping-follow-up-contrasts-for-within-subject-anovas-part-2/)
# Also used this helpful page to re-level the order of comparisons in emmeans (for some reason it always considers domestic the 'control' group, so estimates are negative rather than positive). 
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
# For example - 
normal <- pairs(emmeans(topmodel.ps_nomean, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.ps_nomean, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.ps2_nomean, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
# Pull out the betas for each leaf stage + their 95% CI 
summ <- data.frame(summary(b1)); summ$leaf_stage <- rownames(summ); summ
names(summ)[names(summ) == 'original'] <- 'beta'
conf <- data.frame(confint(b1, type = "perc")); conf$leaf_stage <- rownames(conf); conf
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci0 <- merge(summ, conf, by="leaf_stage")
ci <- subset(ci0, select=c("leaf_stage", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot")); ci

##### F) Calculate % changes
# This is a little different, because of the log-transformation: See here
# https://stats.stackexchange.com/questions/362556/linear-mixed-effect-model-interpretation-with-log-transformed-dependent-variable; also https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
# "In addition, you should interpret your transformed beta coefficients (and associated confidence intervals) as multiplicative changes once you transform to the original units. E.g., beta from the model = 0.0085; exp(.008585) ~= 1.0086 = 0.86% change"... "Only the dependent/response variable is log-transformed. Exponentiate the coefficient, subtract one from this number, and multiply by 100. This gives the percent increase (or decrease) in the response for every one-unit increase in the independent variable. Example: the coefficient is 0.198. (exp(0.198) – 1) * 100 = 21.9. For every one-unit increase in the independent variable, our dependent variable increases by about 22%."
change.ps_nomean <- merge(pred_sig, ci, by="leaf_stage"); change.ps_nomean
change.ps_nomean$exp.beta <- exp(change.ps_nomean$beta)
change.ps_nomean$exp.beta.lower <- exp(change.ps_nomean$beta.CI_lower.boot)
change.ps_nomean$exp.beta.upper <- exp(change.ps_nomean$beta.CI_upper.boot)
change.ps_nomean
change.ps_nomean$percent.change <- (change.ps_nomean$exp.beta-1)*100
change.ps_nomean$percent.change.lower <-(change.ps_nomean$exp.beta.lower-1)*100
change.ps_nomean$percent.change.upper <- (change.ps_nomean$exp.beta.upper-1)*100

change.ps_nomean$response <- "sd"
change.ps_nomean$scale <- "per-age-class"
change.ps_nomean$model <- "status_only"

##### G) Put it all together & re-order for merging
change.ps_nomean <- change.ps_nomean[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "wild.pred_conf.low", "wild.pred_conf.high", "dom.pred", "dom.pred_conf.low", "dom.pred_conf.high", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.ps_nomean

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#
# Finally, fit models to show shift in mean #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# What is the change in mean by leaf age class with domestication?
##### A) Tweak ranFX structure NA
##### B) Tweak fixed FX structure
m1 <- glmmTMB(log(mean) ~ dom_wild*leaf_stage + integration_method + (1|compound) + (1|plant_pop), data=data.ps)
m2 <- glmmTMB(log(mean) ~ dom_wild*leaf_stage + (1|compound) + (1|plant_pop), data=data.ps)
m3 <- glmmTMB(log(mean) ~ dom_wild + leaf_stage + integration_method + (1|compound) + (1|plant_pop), data=data.ps)
m4 <- glmmTMB(log(mean) ~ dom_wild + leaf_stage + (1|compound) + (1|plant_pop), data=data.ps)
AICctab(m1, m2, m3, m4) # m3 best

#Check if ancestry improves model
m3_both <- glmmTMB(log(mean) ~ dom_wild + leaf_stage + integration_method + percent_CAER + percent_FAL + (1|compound) + (1|plant_pop), data=data.ps)
m3_caer <- glmmTMB(log(mean) ~ dom_wild + leaf_stage + integration_method + percent_CAER + (1|compound) + (1|plant_pop), data=data.ps)
m3_fal <- glmmTMB(log(mean) ~ dom_wild + leaf_stage + integration_method + percent_FAL + (1|compound) + (1|plant_pop), data=data.ps)
AICctab(m3, m3_both, m3_caer, m3_fal) # m3 is best - no ancestry

topmodel.mean <- glmmTMB(log(mean) ~ dom_wild + leaf_stage + integration_method + (1|compound) + (1|plant_pop), data=data.ps)
# For bootMer
topmodel.mean2 <- glmmTMB(mean.log ~ dom_wild + leaf_stage + integration_method + (1|compound) + (1|plant_pop), data=data.ps)
hist(resid(topmodel.mean))
plot(residuals(topmodel.mean) ~ predict(topmodel.mean))

##### C) Pull out predicted means for wild & domestic + model nobs + df
model.ests <- data.frame(ggpredict(topmodel.mean, terms=c("dom_wild", "leaf_stage")))
ggplot(model.ests, aes(x=x, y=predicted)) +
  facet_grid(~group) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), size=.5, width=.2, color="black")

##---- Young ----##
wild.pred <- model.ests$predicted[which(model.ests$x == "wild" & model.ests$group == "young")]; wild.pred
wild.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "wild" & model.ests$group == "young")]; wild.pred_conf.low
wild.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "wild" & model.ests$group == "young")]; wild.pred_conf.high
dom.pred <- model.ests$predicted[which(model.ests$x == "domestic" & model.ests$group == "young")]; dom.pred
dom.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "domestic" & model.ests$group == "young")]; dom.pred_conf.low
dom.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "domestic" & model.ests$group == "young")]; dom.pred_conf.high
young <- data.frame(wild.pred, wild.pred_conf.low, wild.pred_conf.high, dom.pred, dom.pred_conf.low, dom.pred_conf.high)
young$leaf_stage <- "young"

##---- Middle ----##
wild.pred <- model.ests$predicted[which(model.ests$x == "wild" & model.ests$group == "middle")]; wild.pred
wild.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "wild" & model.ests$group == "middle")]; wild.pred_conf.low
wild.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "wild" & model.ests$group == "middle")]; wild.pred_conf.high
dom.pred <- model.ests$predicted[which(model.ests$x == "domestic" & model.ests$group == "middle")]; dom.pred
dom.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "domestic" & model.ests$group == "middle")]; dom.pred_conf.low
dom.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "domestic" & model.ests$group == "middle")]; dom.pred_conf.high
middle <- data.frame(wild.pred, wild.pred_conf.low, wild.pred_conf.high, dom.pred, dom.pred_conf.low, dom.pred_conf.high)
middle$leaf_stage <- "middle"

##---- Old ----##
wild.pred <- model.ests$predicted[which(model.ests$x == "wild" & model.ests$group == "old")]; wild.pred
wild.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "wild" & model.ests$group == "old")]; wild.pred_conf.low
wild.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "wild" & model.ests$group == "old")]; wild.pred_conf.high
dom.pred <- model.ests$predicted[which(model.ests$x == "domestic" & model.ests$group == "old")]; dom.pred
dom.pred_conf.low <- model.ests$conf.low[which(model.ests$x == "domestic" & model.ests$group == "old")]; dom.pred_conf.low
dom.pred_conf.high <- model.ests$conf.high[which(model.ests$x == "domestic" & model.ests$group == "old")]; dom.pred_conf.high
old <- data.frame(wild.pred, wild.pred_conf.low, wild.pred_conf.high, dom.pred, dom.pred_conf.low, dom.pred_conf.high)
old$leaf_stage <- "old"

ests <- rbind(young,middle,old); ests
# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.mean)
ests$model_df.residual <- df.residual(topmodel.mean)

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.mean, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig # same for all leaf stages

# Put together the model-predicted means + significance
pred_sig <- merge(ests, sig) # same for each leaf stage (no interaction), so can just merge across

##### E) Bootstrap CI 
# Here, instead of extracting the single fixed effect for domestication across all leaf stages, I need to extract the leaf stage-specific effects due to the significant interaction 
# (Used this great post as a resource: https://www.r-bloggers.com/bootstrapping-follow-up-contrasts-for-within-subject-anovas-part-2/)
# Also used this helpful page to re-level the order of comparisons in emmeans (for some reason it always considers domestic the 'control' group, so estimates are negative rather than positive). 
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
# For example - 
normal <- pairs(emmeans(topmodel.mean, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.mean, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.mean2, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
# Pull out the betas for each leaf stage + their 95% CI 
summ <- data.frame(summary(b1)); summ$leaf_stage <- rownames(summ); summ
names(summ)[names(summ) == 'original'] <- 'beta'
conf <- data.frame(confint(b1, type = "perc")); conf$leaf_stage <- rownames(conf); conf
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci0 <- merge(summ, conf, by="leaf_stage")
ci <- subset(ci0, select=c("leaf_stage", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot")); ci

##### F) Calculate % changes
# This is a little different, because of the log-transformation: See here
# https://stats.stackexchange.com/questions/362556/linear-mixed-effect-model-interpretation-with-log-transformed-dependent-variable; also https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
# "In addition, you should interpret your transformed beta coefficients (and associated confidence intervals) as multiplicative changes once you transform to the original units. E.g., beta from the model = 0.0085; exp(.008585) ~= 1.0086 = 0.86% change"... "Only the dependent/response variable is log-transformed. Exponentiate the coefficient, subtract one from this number, and multiply by 100. This gives the percent increase (or decrease) in the response for every one-unit increase in the independent variable. Example: the coefficient is 0.198. (exp(0.198) – 1) * 100 = 21.9. For every one-unit increase in the independent variable, our dependent variable increases by about 22%."
change.mean <- merge(pred_sig, ci, by="leaf_stage"); change.mean
change.mean$exp.beta <- exp(change.mean$beta)
change.mean$exp.beta.lower <- exp(change.mean$beta.CI_lower.boot)
change.mean$exp.beta.upper <- exp(change.mean$beta.CI_upper.boot)
change.mean
change.mean$percent.change <- (change.mean$exp.beta-1)*100
change.mean$percent.change.lower <-(change.mean$exp.beta.lower-1)*100
change.mean$percent.change.upper <- (change.mean$exp.beta.upper-1)*100

change.mean$response <- "mean"
change.mean$scale <- "per-age-class"
change.mean$model <- "status_only"

##### G) Put it all together & re-order for merging
change.mean <- change.mean[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "wild.pred_conf.low", "wild.pred_conf.high", "dom.pred", "dom.pred_conf.low", "dom.pred_conf.high", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.mean

#-------------------------------------------------------------------#
# Rbind all % changes for variability quantified within leaf stage  #
#-------------------------------------------------------------------#
estimates.ps <- rbind(change.ps, change.ps_nomean, change.mean); estimates.ps

#------------------------------------------------------------------#
#------------------------------------------------------------------#
# Rbind the two scales together : whole-plant + across all leaves  #
#------------------------------------------------------------------#
#------------------------------------------------------------------#
both <- rbind(estimates.pp, estimates.ps); both

# Rename / Relevel some variables for plotting
both$hypothesis <- "Var ~ Domestication + Mean"
both$hypothesis[which(both$model == "status_only" & both$response == "mean")] <- "Mean ~ Domestication"
both$hypothesis[which(both$model == "status_only" & both$response == "sd")] <- "Var ~ Domestication"
both$hypothesis <- factor(both$hypothesis, levels=c("Var ~ Domestication", "Var ~ Domestication + Mean", "Mean ~ Domestication"))
both$leaf_stage <- factor(both$leaf_stage, levels=c("all", "young", "middle", "old"))
# Add a character column for significance
# Add a character column for significance
# . = trend, 0.1 level; * = .05 level; ** = .01; *** = <.001
both$sig.asterisk <- ""
both$sig.asterisk[which(both$p.value <.15)] <- "."
both$sig.asterisk[which(both$p.value <.055)] <- "*"
both$sig.asterisk[which(both$p.value <.015)] <- "**"
both$sig.asterisk[which(both$p.value <.0015)] <- "***"
both

# Add the trait 
both$trait <- "saponins_84"
# Rename
saps <- both; head(saps)

# Quick check that a plot looks OK
ggplot(data=saps, aes(x=leaf_stage, y=percent.change, group=hypothesis)) +
  geom_hline(yintercept=0, lty=2, color = "gray40") +
  geom_errorbar(aes(ymin=percent.change.lower, ymax=percent.change.upper, color=leaf_stage), size=4, alpha=0.3, width=0, position=position_dodge(width=.5)) +
  geom_point(aes(color=leaf_stage, shape=hypothesis, size=hypothesis), stroke=1.5, position=position_dodge(width=.5)) +
  scale_shape_manual("Response", values=c("Var ~ Domestication"=16, "Var ~ Domestication + Mean"=21, "Mean ~ Domestication"=15)) +
  scale_size_manual(values=c("Var ~ Domestication"=3.5, "Var ~ Domestication + Mean"=2.7, "Mean ~ Domestication"=2.7)) +
  scale_color_manual("Scale of Variability", values=c("young" = "yellowgreen", "middle" = "mediumseagreen", "old" = "forestgreen", "all" = "burlywood3")) +
  theme_bw()+
  geom_text(aes(label = sig.asterisk, y=max(saps$percent.change.upper)+25), size=5, color="gray32", position=position_dodge(width=.5), angle=30) +
  theme(panel.border = element_rect(colour="black"), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), strip.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Scale of Variability") +
  ylab("Predicted % Change\nWild ---> Domestic") +
  scale_x_discrete(labels = c("All leaves","Apical\n(youngest)", "Expanded\n(middle)", "Basal\n(oldest)")) +
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=10), strip.text.x=element_text(size=15, face="bold"), axis.title.y=element_text(size=15,face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x=element_text(size=15,face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  guides(color=FALSE, size=FALSE)
