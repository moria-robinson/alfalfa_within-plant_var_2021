#######################################
# Leaf Water Content (LWC) - analyses #
#######################################
# Goals: 
# 1) Calculate variability (sd) at the whole-plant and the per-leaf-age-class scale
# 2) Fit models to estimate effect of domestication on variability at each of these scales
# 3) Pull out % change and Credible Interval for each comparison of interest
# 4) Export those values for use in figure

######## packages #############
library(tidyr)
library(bbmle)
library(MuMIn)
library(DHARMa)
library(lme4)
library(boot)
library(glmmTMB)
library(ggeffects)
library(ggplot2)
library(emmeans)
library(sjPlot)

####### DATA ########
data <- read.csv("~/lwc_data.csv")

head(data); str(data)
## Right-away: re-level the factors of interest for comparisons later
unique(data$dom_wild)
data$dom_wild <- factor(data$dom_wild, levels=c("wild", "domestic")); head(data)
unique(data$leaf_stage)
data$leaf_stage <- factor(data$leaf_stage, levels=c("young", "middle", "old")); head(data)

##################################################################
# Next, calculate mean & SD 
# Do so at two scales: across all leaves & within leaf age class
##################################################################
head(data)
#-------------#
#  Per-plant  #
#-------------#
# across all leaves of all relative age classes); keep relevant plant-level covariates
data.pp <- do.call(data.frame, aggregate(LWC ~ plant_pop + dom_wild + plant.volume_centered + repro.stage + repro.stage.simp + percent_FAL + percent_HEM + percent_CAER, data=data, function(x) c(mean = mean(x), sd = sd(x), N = length(x)))); head(data.pp); nrow(data.pp) # 60 - one measure of mean and sd per plant
names(data.pp) <- c("plant_pop", "dom_wild", "plant.volume_centered", "repro.stage", "repro.stage.simp", "percent_FAL","percent_HEM","percent_CAER","mean", "sd", "N")
head(data.pp); str(data.pp)

#-------------------------------#
#  Per leaf relative age class  #
#-------------------------------#
# across all leaves of all relative age classes); keep relevant plant-level covariates
data.ps <- do.call(data.frame, aggregate(LWC ~ plant_pop + dom_wild + plant.volume_centered + repro.stage + repro.stage.simp + leaf_stage + percent_FAL + percent_HEM + percent_CAER, data=data, function(x) c(mean = mean(x), sd = sd(x), N = length(x)))); head(data.ps); nrow(data.ps) # 180 - 3 measures of mean and sd per plant (one per relative leaf age class)
names(data.ps) <- c("plant_pop", "dom_wild", "plant.volume_centered", "repro.stage", "repro.stage.simp", "leaf_stage","percent_FAL","percent_HEM","percent_CAER", "mean", "sd", "N")
head(data.ps); str(data.ps)

#####################################################
# 1) Model fitting: whole-plant / across all leaves #
#####################################################
head(data.pp); str(data.pp)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# First, fit models with both mean and domestication status #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# What is the change in variability with domestication, accounting for effects of the mean?

##### A) Tweak ranFX structure
m1.1 <- glmmTMB(sd ~ mean + dom_wild + plant.volume_centered + (1|repro.stage), data=data.pp, na.action = "na.fail")
m1.2 <- glmmTMB(sd ~ mean + dom_wild + plant.volume_centered + (1|repro.stage.simp), data=data.pp, na.action = "na.fail")
AICctab(m1.1, m1.2)
# Exactly the same.
# Make uninformative dummy random effects to compare repro.stage with, decide to keep or not in model
data.pp$repro.stage_dummy <- as.factor(sample(unique(data.pp$repro.stage), replace=T, nrow(data.pp))); head(data.pp)
m1.3 <- glmmTMB(sd ~ mean + dom_wild + plant.volume_centered + (1|repro.stage_dummy), data=data.pp, na.action = "na.fail")
AICctab(m1.1, m1.3) # also the same. So, having reproductive stage doesn't add anything to the model

##### B) Tweak fixed FX structure
m1.4 <- glmmTMB(sd ~ mean*dom_wild*plant.volume_centered + percent_FAL + percent_CAER, data=data.pp, na.action = "na.fail")
dredge(m1.4) # Top models (within 2AIC) have mean and dom_wild, no plant size + no interaction
topmodel.pp <- glmmTMB(sd ~ mean + dom_wild, data=data.pp, na.action = "na.fail")
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.pp))

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
fixef(topmodel.pp)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object
b1 <- lme4::bootMer(topmodel.pp, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")
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
change.pp <- cbind(pred_sig, ci); change.pp
change.pp$percent.change <- change.pp$beta/change.pp$wild.pred*100
change.pp$percent.change.lower <- change.pp$beta.CI_lower.boot/change.pp$wild.pred*100
change.pp$percent.change.upper <- change.pp$beta.CI_upper.boot/change.pp$wild.pred*100

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
topmodel.pp_nomean <- glmmTMB(sd ~ dom_wild, data=data.pp, na.action = "na.fail")
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.pp_nomean))

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
fixef(topmodel.pp_nomean)
fixef(topmodel.pp_nomean)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object
b1 <- lme4::bootMer(topmodel.pp_nomean, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")
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
change.pp_nomean <- cbind(pred_sig, ci); change.pp_nomean
change.pp_nomean$percent.change <- change.pp_nomean$beta/change.pp_nomean$wild.pred*100
change.pp_nomean$percent.change.lower <- change.pp_nomean$beta.CI_lower.boot/change.pp_nomean$wild.pred*100
change.pp_nomean$percent.change.upper <- change.pp_nomean$beta.CI_upper.boot/change.pp_nomean$wild.pred*100

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

##### A) Tweak random FX structure
m3.1 <- glmmTMB(mean ~ dom_wild + plant.volume_centered + (1|repro.stage), data=data.pp, na.action = "na.fail")
m3.2 <- glmmTMB(mean ~ dom_wild + plant.volume_centered + (1|repro.stage.simp), data=data.pp, na.action = "na.fail")
AICctab(m3.1, m3.2)
# < 2AIC difference
# Make uninformative dummy random effects to compare repro.stage with, decide to keep or not in model
m3.3 <- glmmTMB(mean ~ dom_wild + plant.volume_centered + (1|repro.stage_dummy), data=data.pp, na.action = "na.fail")
AICctab(m3.1, m3.3) # also the same; adding reproductive stage doesn't improve model

##### B) Tweak fixed FX structure
m4.1 <- glmmTMB(mean ~ dom_wild*plant.volume_centered + percent_FAL + percent_CAER, data=data.pp, na.action = "na.fail")
dredge(m4.1) # Top model (by > 2AIC) has only dom_wild, no plant size
topmodel.mean <- glmmTMB(mean ~ dom_wild, data=data.pp, na.action = "na.fail")
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.mean))

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
fixef(topmodel.mean)
fixef(topmodel.mean)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object
b1 <- lme4::bootMer(topmodel.mean, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")
boot <- boot.ci(b1,type="perc")

# Pull out the effect (beta) from the model
beta <- as.numeric(fixef(topmodel.mean)$cond["dom_wilddomestic"])
# Pull out the bootstrapped 95% CI around that beta. I need to get this from the "percent" part of the bootstrapping output. The items within aren't named, which is annoying. From bootMer, below "These latter four components will be matrices with 5 columns, the first column containing the level, the next two containing the indices of the order statistics used in the calculations and the final two the calculated endpoints themselves." So I need to pull the latter two.
boot; str(boot); class(boot)
out <- data.frame(boot$percent); out # so the last two components are the upper and lower CIs, respectively
beta.CI_upper.boot <- as.numeric(out[length(out)]); beta.CI_upper.boot 
beta.CI_lower.boot <- as.numeric(out[length(out)-1]); beta.CI_lower.boot
# Put them together
ci <- data.frame(beta, beta.CI_lower.boot, beta.CI_upper.boot)

##### F) Calculate % changes
change.mean <- cbind(pred_sig, ci); change.mean
change.mean$percent.change <- change.mean$beta/change.mean$wild.pred*100
change.mean$percent.change.lower <- change.mean$beta.CI_lower.boot/change.mean$wild.pred*100
change.mean$percent.change.upper <- change.mean$beta.CI_upper.boot/change.mean$wild.pred*100

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

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# First, fit models with both mean and domestication status #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# What is the change in variability within each leaf age class with domestication, accounting for effects of the mean?

##### A) Tweak random FX structure
m5.1 <- glmmTMB(sd ~ dom_wild + leaf_stage + mean + plant.volume_centered + (1|plant_pop) + (1|repro.stage), data=data.ps) # did not converge
m5.2 <- glmmTMB(sd ~ dom_wild + leaf_stage + mean + plant.volume_centered + (1|plant_pop), data=data.ps) 
m5.3 <- glmmTMB(sd ~ dom_wild + leaf_stage + mean + plant.volume_centered + (1|repro.stage), data=data.ps) 
AICctab(m5.2, m5.3) # m5.2 and m5.3 are the same; use plant_pop

##### B) Tweak fixed FX structure
m5.4 <- glmmTMB(sd ~ dom_wild*leaf_stage*mean + plant.volume_centered + percent_FAL + percent_CAER + (1|plant_pop), data=data.ps)
dredge(m5.4) # top models
topmodel.ps <- glmmTMB(sd ~ dom_wild + leaf_stage + mean + (1|plant_pop), data=data.ps) # top models (within 2AIC) had all the fixed effects except plant size, and 3rd-best has a leaf stage*mean interaction, but not an interaction of interest.  
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.ps))

##### C) Pull out predicted means for wild & domestic + model nobs + df
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
# No interactions, so same effect of dom_wild when SD is quantified within leaf age, regardless of age class
sig1 <- data.frame(car::Anova(topmodel.ps, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- merge(ests, sig)
# pred_sig <- merge(ests, sig, by="leaf_stage") # don't merge by leaf stage; same significance for each leaf stage (no intxn)

##### E) Bootstrap CI 
# No interaction with leaf age, so once again just pull the effect of domestication quantified at this scale
fixef(topmodel.ps)
fixef(topmodel.ps)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object
b1 <- lme4::bootMer(topmodel.ps, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")
boot <- boot.ci(b1,type="perc")

# Pull out the effect (beta) from the model
beta <- as.numeric(fixef(topmodel.ps)$cond["dom_wilddomestic"])
# Pull out the bootstrapped 95% CI around that beta. I need to get this from the "percent" part of the bootstrapping output. The items within aren't named, which is annoying. From bootMer, below "These latter four components will be matrices with 5 columns, the first column containing the level, the next two containing the indices of the order statistics used in the calculations and the final two the calculated endpoints themselves." So I need to pull the latter two.
boot; str(boot); class(boot)
out <- data.frame(boot$percent); out # so the last two components are the upper and lower CIs, respecively
beta.CI_upper.boot <- as.numeric(out[length(out)]); beta.CI_upper.boot 
beta.CI_lower.boot <- as.numeric(out[length(out)-1]); beta.CI_lower.boot
# Put them together
ci <- data.frame(beta, beta.CI_lower.boot, beta.CI_upper.boot)

##### F) Calculate % changes
change.ps <- cbind(pred_sig, ci); change.ps
change.ps$percent.change <- change.ps$beta/change.ps$wild.pred*100
change.ps$percent.change.lower <- change.ps$beta.CI_lower.boot/change.ps$wild.pred*100
change.ps$percent.change.upper <- change.ps$beta.CI_upper.boot/change.ps$wild.pred*100

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
topmodel.ps_nomean <- glmmTMB(sd ~ dom_wild*leaf_stage + (1|plant_pop), data=data.ps) # Here, there is a significant interaction with leaf stage
summary(topmodel.ps_nomean)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.ps_nomean))

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
b1 <- lme4::bootMer(topmodel.ps_nomean, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
# Pull out the betas for each leaf stage + their 95% CI 
summ <- data.frame(summary(b1)); summ$leaf_stage <- rownames(summ); summ
names(summ)[names(summ) == 'original'] <- 'beta'
conf <- data.frame(confint(b1, type = "perc")); conf$leaf_stage <- rownames(conf); conf
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci0 <- merge(summ, conf, by="leaf_stage")
ci <- subset(ci0, select=c("leaf_stage", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot")); ci

##### F) Calculate % changes
change.ps_nomean <- cbind(pred_sig, ci); change.ps_nomean
change.ps_nomean$percent.change <- change.ps_nomean$beta/change.ps_nomean$wild.pred*100
change.ps_nomean$percent.change.lower <- change.ps_nomean$beta.CI_lower.boot/change.ps_nomean$wild.pred*100
change.ps_nomean$percent.change.upper <- change.ps_nomean$beta.CI_upper.boot/change.ps_nomean$wild.pred*100

change.ps_nomean$response <- "sd"
change.ps_nomean$scale <- "per-age-class"
change.ps_nomean$model <- "status_only"

##### G) Put it all together & re-order for merging
change.ps_nomean <- change.ps_nomean[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "wild.pred_conf.low", "wild.pred_conf.high", "dom.pred", "dom.pred_conf.low", "dom.pred_conf.high", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.ps_nomean

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#
# Finally, fit models to show shift in mean #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# What is the change in mean by leaf agge class with domestication?
##### A) Tweak random FX structure
m6.1 <- glmmTMB(mean ~ dom_wild + leaf_stage + plant.volume_centered + (1|plant_pop) + (1|repro.stage), data=data.ps, na.action = "na.fail")
m6.2 <- glmmTMB(mean ~ dom_wild + leaf_stage + plant.volume_centered + (1|repro.stage), data=data.ps, na.action = "na.fail")
m6.3 <- glmmTMB(mean ~ dom_wild + leaf_stage + plant.volume_centered + (1|plant_pop), data=data.ps, na.action = "na.fail")
AICctab(m6.1, m6.2, m6.3) # just the model with plant_pop best > 2AIC difference

##### B) Tweak fixed FX structure
m6.4 <- glmmTMB(mean ~ dom_wild*leaf_stage + percent_FAL + percent_CAER + (1|plant_pop), data=data.ps, na.action = "na.fail")
dredge(m6.4) # Top model has only dom_wild + leaf stage + their interaction, and no plant size
topmodel.mean <- glmmTMB(mean ~ dom_wild*leaf_stage + (1|plant_pop), data=data.ps, na.action = "na.fail")
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.mean))

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
sig0 <- data.frame(pairs(emmeans(topmodel.mean, ~ dom_wild|leaf_stage))) # post-hoc test 
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
normal <- pairs(emmeans(topmodel.mean, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.mean, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.mean, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
# Pull out the betas for each leaf stage + their 95% CI 
summ <- data.frame(summary(b1)); summ$leaf_stage <- rownames(summ); summ
names(summ)[names(summ) == 'original'] <- 'beta'
conf <- data.frame(confint(b1, type = "perc")); conf$leaf_stage <- rownames(conf); conf
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci0 <- merge(summ, conf, by="leaf_stage")
ci <- subset(ci0, select=c("leaf_stage", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot")); ci

##### F) Calculate % changes
change.mean <- merge(pred_sig, ci, by="leaf_stage"); change.mean
change.mean$percent.change <- change.mean$beta/change.mean$wild.pred*100
change.mean$percent.change.lower <- change.mean$beta.CI_lower.boot/change.mean$wild.pred*100
change.mean$percent.change.upper <- change.mean$beta.CI_upper.boot/change.mean$wild.pred*100

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
both$trait <- "lwc"
# Rename
lwc <- both; head(lwc)

# Quick check that a plot looks OK
ggplot(data=lwc, aes(x=leaf_stage, y=percent.change, group=hypothesis)) +
  geom_hline(yintercept=0, lty=2, color = "gray40") +
  geom_errorbar(aes(ymin=percent.change.lower, ymax=percent.change.upper, color=leaf_stage), size=4, alpha=0.3, width=0, position=position_dodge(width=.5)) +
  geom_point(aes(color=leaf_stage, shape=hypothesis, size=hypothesis), stroke=1.5, position=position_dodge(width=.5)) +
  scale_shape_manual("Response", values=c("Var ~ Domestication"=16, "Var ~ Domestication + Mean"=21, "Mean ~ Domestication"=15)) +
  scale_size_manual(values=c("Var ~ Domestication"=3.5, "Var ~ Domestication + Mean"=2.7, "Mean ~ Domestication"=2.7)) +
  scale_color_manual("Scale of Variability", values=c("young" = "yellowgreen", "middle" = "mediumseagreen", "old" = "forestgreen", "all" = "burlywood3")) +
  theme_bw()+
  geom_text(aes(label = sig.asterisk, y=max(lwc$percent.change.upper)+25), size=5, color="gray32", position=position_dodge(width=.5), angle=30) +
  theme(panel.border = element_rect(colour="black"), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), strip.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Scale of Variability") +
  ylab("Predicted % Change\nWild ---> Domestic") +
  scale_x_discrete(labels = c("All leaves","Apical\n(youngest)", "Expanded\n(middle)", "Basal\n(oldest)")) +
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=10), strip.text.x=element_text(size=15, face="bold"), axis.title.y=element_text(size=15,face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x=element_text(size=15,face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  guides(color=FALSE, size=FALSE)

