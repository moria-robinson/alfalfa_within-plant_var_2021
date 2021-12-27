##########################################
# Diversity of compounds - saponins only # 
##########################################
# Alpha / gamma diversity: chemical diversity of plants at two scales. 

# 1) Aggregated across all leaves (=gamma); yields a single value of diversity for a single plant. Abundance covariate = single measure of total compound abundance produced by the plant. Could use total abundance of all compounds across all leaves, or the average across all leaves (should give same answer, latter is just divvied by 9). Easier to use average, same measure used to quantify beta div across all leaves.

# 2) Within a single leaf (=alpha); yields a single value of diversity for each leaf. Thus, as an abundance covariate, most power / makes most sense to use the single value of total compound abundance per leaf. Before, was using the average among young, middle, or older-aged leaves; but this maps the same value of total compound abundance to 3 unique values of diversity, which isn't right.

# These questions get at two very related questions, though at two scales; do plant individuals that produce more compounds overall produce more compounds overall, and do individual leaves that produce more compounds overall produce more diverse cocktails of compounds. Presumably these two scales should parallel each other - but perhaps within leaf age class, there are different strengths of relationship between overall production & diversity, etc.

# As the total abundance covariate, models use average leaf areas per plant for (gamma), same as for beta diversity. But for alpha diversity, instead of the average within leaf age class, since each leaf has a unique value of alpha, we can nget more info by using their corresponding unique total concentration (do leaves with more compounds have more diversity). So, for alpha, we simply use the total raw peak intensities per leaf.

# No transformation here, as I want to pick up differences in evenness (and not exclude zeros, which log-transformation would require).

# Thoughts for text, and reason to analyze this way; if there is a difference in diversity (e.g. richness, Shannon) - is this because one group of plants makes more compounds overall (e.g. effect goes away with total peak area as a covariate), or does the group have greater diversity per unit of compound abundance? 
# Note that this is per unit of tissue, since solvent & digitoxin were added in proportion to the amount of tissue. So, from the herbivore's perspective, the overall diversity results (e.g. without total peak area in the model) describe a "real" difference between a wild and domestic plant, in a comparable bite. Even if total peak area removes the domestication effect, this simply means that herbivores are experiencing both greater diversity and greater amounts simultaneously; similarly to when mean and SD shift jointly for univariate traits. The model comparisons are getting more at mechanisms that generate chemical diversity within the plant. 

# 1) TOTAL PEAK AREA
#      1a) Per-plant. Average total per-leaf peak area, averaged across all leaves per-plant. Compare total peak area between wild and domestic plants (this is the average value used as a covariate in models of gamma div + beta across all leaves)
#      2a) Per leaf age class-1. Total per-leaf peak area. Compare in the same way we would alpha (this is the only 'mean' where we actually have one value per leaf; so it's not really a mean, but rather the appropriate scale to look at abundance ~ diversity relationships for alpha. Need to make sure to explain well in supp.)

# 2) RICHNESS. For each scale (per-plant, per-leaf age class) run models with and without total peak area to explore relationships described above. 
#      2a) - Per-plant. GAMMA RICHNESS (pooled across all leaves)
#      2b) - Per leaf age class. ALPHA RICHNESS (per-leaf)

# 3) SHANNON DIVERSITY. For each scale (per-plant, per-leaf age class), run models with and without total peak area to explore relationships described above. 
#      3a) - Per-plant. GAMMA SHANNON DIVERSITY (pooled across all leaves)
#      3b) - Per leaf age class. ALPHA SHANNON DIVERSITY (per-leaf)

############ Packages ############
library(lme4)
library(lmerTest)
library(vegan)
library(labdsv)
library(Rmisc)
library(ggplot2)
library(car)
library(tidyr)
library(glmmTMB)
library(MuMIn)
library(ggeffects)
library(emmeans)
library(bbmle)
library(sjPlot)
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

# Re-level domestic & wild for downstream comparisons
data_saps2$dom_wild <- factor(data_saps2$dom_wild, levels=c("wild", "domestic"))
# Re-level leaf stage for downstream comparisons
data_saps2$leaf_stage <- factor(data_saps2$leaf_stage, levels=c("young", "middle", "old"))

data <- data_saps2; nrow(data_saps2); nrow(data) # 46440

########################################################
# Calculate diversity - per plant & per leaf age class #
########################################################
head(data); str(data); length(unique(data$sample_id)); length(unique(data$compound)) # 540 leaves. 9 per plant.

#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Caculate gamma diversity #
#-#-#-#-#-#-#-#-#-#-#-#-#-#
# sum up raw concentrations of each compound at the level of the plant (add all the leaves together)
data.pp <- aggregate(peak_intensity_norm_noNoise ~ plant_pop + compound, data=data, sum); head(data.pp); nrow(data.pp) # 86 saponins*60 plants = 5160 rows

# Add a presence/absence column (to tally for richness)
data.pp$present.absent <- 0
data.pp$present.absent[which(data.pp$peak_intensity_norm_noNoise != 0)] <- 1
head(data.pp); data.pp$present.absent

# For each plant, calculate: shannon diversity, simpson's diversity, compound richness. Also exp(Shannon) - effective numbers - per Matt's suggestion. Calculate all with raw abundances
# First, get data into correct form to calculate Shannon, Simpson
temp1 <- subset(data.pp, select=c("plant_pop", "compound", "peak_intensity_norm_noNoise")); str(temp1)
div <- matrify(temp1); div[c(1:10), c(1:10)]

# Shannon
H <- diversity(div); head(H); str(H)

# Simpson
simp <- diversity(div, "simpson"); str(simp)

# Richness
temp2 <- subset(data.pp, select=c("plant_pop", "compound", "present.absent"))
rich <- matrify(temp2)
richness <- rowSums(rich); str(richness)

# Cbind all the diversity metrics per sample
diversity1 <- data.frame(cbind(H, simp = simp[names(H)], richness = richness[names(H)])); head(diversity1)

# Effective numbers (I think less important, because % change should be the same as if we use Shannon, just harder to interpret in terms of original units)
diversity1$expH <- exp(diversity1$H); head(diversity1)

# Save plant_pop ID
diversity1$plant_pop <- rownames(diversity1); head(diversity1); rownames(diversity1) <- NULL
head(diversity1)

# For gamma diversity, we will use the average per-leaf total compound concentrations as the mean covariate.
# First, calculate the total compound concentrations per leaf. This will be used as a covariate for the alpha div measures.
area_pL <- do.call(data.frame, (aggregate(peak_intensity_norm_noNoise ~ plant_pop + dom_wild + sample_id + leaf_stage + plant.volume_centered + repro.stage.simp + percent_FAL + percent_HEM + percent_CAER, data=data, function(x) c(sum = sum(x), N = length(x))))); head(area_pL); str(area_pL); nrow(area_pL) # 540 rows
head(area_pL) # N columnn - good, summed across 86 compounds per leaf
area_pL$peak_intensity_norm_noNoise.N <- NULL
names(area_pL) <- c("plant_pop", "dom_wild", "sample_id", "leaf_stage", "plant.volume_centered", "repro.stage.simp", "percent_FAL","percent_HEM","percent_CAER", "total_peak_area_perLeaf") 
head(area_pL); nrow(area_pL) # 540 measures; one per leaf

# Next, get the mean peak area per plant, averaged across the 9 leaves/plant. This mean is the covariate to be used with gamma diversity models. 
mean.peak.area_raw.pp <- do.call(data.frame, (aggregate(total_peak_area_perLeaf ~ plant_pop + dom_wild + plant.volume_centered + repro.stage.simp + percent_FAL + percent_HEM + percent_CAER, data=area_pL, function(x) c(mean = mean(x), N = length(x))))); head(mean.peak.area_raw.pp); str(mean.peak.area_raw.pp); nrow(mean.peak.area_raw.pp) # 60 rows; one per plant; check that the N=9. Yes!
mean.peak.area_raw.pp$total_peak_area_perLeaf.N <- NULL; names(mean.peak.area_raw.pp) <- c("plant_pop", "dom_wild", "plant.volume_centered", "repro.stage.simp", "percent_FAL","percent_HEM","percent_CAER", "total_peak_area_perLeaf_mean"); head(mean.peak.area_raw.pp); nrow(mean.peak.area_raw.pp) 

# NOTE: mean peak area per age class is calculated in the next section, and used to report effects of domestication on total compound production per leaf age class. But it's not used as a covariate for anything in this script; only in the beta diversity analyses per plant/per leaf age class.

# Work with this dataframe - all diversity measures calculated per-plant
diversity.pp <- merge(diversity1, mean.peak.area_raw.pp)
head(diversity.pp); nrow(diversity.pp) # 60 (one measure per plant)

#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Caculate alpha diversity #
#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Calculate metrics of compound diversity for each leaf
head(data); str(data)

# Add a presence/absence column (to tally for richness)
data$present.absent <- 0
data$present.absent[which(data$peak_intensity_norm_noNoise != 0)] <- 1
head(data); data$present.absent

# For each single-leaf sample, calculate: shannon diversity, simpsons's diversity, compound richness. Also exp(Shannon) - effective numbers - per Matt's suggestion. Calculate all with raw abundances

# Get data into correct form to calculate Shannon, Simpson per individual leaf
temp1 <- subset(data, select=c("sample_id", "compound", "peak_intensity_norm_noNoise")); str(temp1)
div <- matrify(temp1); div[c(1:10), c(1:10)]

# Shannon
H <- diversity(div); head(H); str(H)

# Simpson
simp <- diversity(div, "simpson"); str(simp)

# Richness
temp2 <- subset(data, select=c("sample_id", "compound", "present.absent"))
rich <- matrify(temp2)
richness <- rowSums(rich); str(richness)

# Cbind all the diversity metrics per sample
diversity2 <- data.frame(cbind(H, simp = simp[names(H)], richness = richness[names(H)])); head(diversity2)

# Add effective numbers 
diversity2$expH <- exp(diversity2$H)

# Save sample_id
diversity2$sample_id <- rownames(diversity2); head(diversity2); rownames(diversity2) <- NULL
head(diversity2); nrow(diversity2)

# Merge with total average concentrations/leaf
head(area_pL); nrow(area_pL) 
diversity.ps <- merge(diversity2, area_pL, by=c("sample_id"))
# Work with this dataframe - all diversity measures calculated per-leaf
head(diversity.ps); nrow(diversity.ps) # 540 (60 plants * 9 leaves/plant); unique alpha div + total concentrations per leaf

#######################################
#######################################
#      TOTAL COMPOUND PRODUCTION      # 
#######################################
#######################################
# In addition to having these compound amounts appended to the diversity dataframe & used as covariates, we will show how the average total compound production per leaf has increased with domestication. We will visualize this in the figure as the average total compound production/leaf averaged across all leaves, and the average among the three leaves/leaf age class. These averages correspond to the covariate used in gamma diversity + beta diversity across all leaves, and beta diversity within leaf age class, respectively. Alpha diversity won't have a perfectly analagous shift-in-means figure, but it's OK I think since this is a secondary focus/scale of the paper e.g. (within rather than among-leaf diversity)

# This will also be useful in the main figure for interpretation of results, as it can also show how making more compounds  overall leads to more variable amounts among leaves. This, in turn, can links to results showing links between total amount and alpha, for example.

head(area_pL); nrow(area_pL) # 540 rows; one per leaf; use total_peak_area_perLeaf

#-------------#
#  Per-plant  #
#-------------#
peak.pp <- do.call(data.frame, aggregate(total_peak_area_perLeaf ~ plant_pop + dom_wild + plant.volume_centered + repro.stage.simp + percent_FAL + percent_HEM + percent_CAER, data=area_pL, function(x) c(mean = mean(x), sd = sd(x), N = length(x)))); head(peak.pp); nrow(peak.pp) # 60 - one measure of mean and sd per plant
names(peak.pp) <- c("plant_pop", "dom_wild", "plant.volume_centered", "repro.stage.simp", "percent_FAL","percent_HEM","percent_CAER", "mean", "sd", "N")
head(peak.pp); str(peak.pp)

#-------------------------------#
#  Per leaf relative age class  #
#-------------------------------#
peak.ps <- do.call(data.frame, aggregate(total_peak_area_perLeaf ~ plant_pop + dom_wild + leaf_stage + plant.volume_centered + repro.stage.simp + percent_FAL + percent_HEM + percent_CAER, data=area_pL, function(x) c(mean = mean(x), sd = sd(x), N = length(x)))); head(peak.ps); nrow(peak.ps) # 180 - one measure of mean and sd per leaf age class, per plant
names(peak.ps) <- c("plant_pop", "dom_wild", "leaf_stage", "plant.volume_centered", "repro.stage.simp","percent_FAL","percent_HEM","percent_CAER", "mean", "sd", "N")
head(peak.ps); str(peak.ps)

#####################################################
# 1) Model fitting: whole-plant / across all leaves #
#####################################################
head(peak.pp); str(peak.pp)
hist(peak.pp$mean) # pretty normal, actually
hist(peak.pp$sd) # more skewed

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# First, fit models with both mean and domestication status #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# What is the change in variability with domestication, accounting for effects of the mean?

##### A) Tweak ranFX structure
m1.1 <- glmmTMB(sd ~ mean*dom_wild*plant.volume_centered + (1|repro.stage.simp), data=peak.pp, na.action = "na.fail")
# Make uninformative dummy random effects to compare repro.stage.simp with, decide to keep or not in model
peak.pp$repro.stage.simp_dummy <- as.factor(sample(unique(peak.pp$repro.stage.simp), replace=T, nrow(peak.pp))); head(peak.pp)
m1.2 <- glmmTMB(sd ~ mean*dom_wild*plant.volume_centered + (1|repro.stage.simp_dummy), data=peak.pp, na.action = "na.fail")
AICctab(m1.1, m1.2) # the same. So, having reproductive stage doesn't add anything to the model
m1.3 <- glmmTMB(sd ~ mean*dom_wild*plant.volume_centered, data=peak.pp, na.action = "na.fail")

##### B) Tweak fixed FX structure
# Check if ancestry improves the model
m1.4 <- glmmTMB(sd ~ mean + dom_wild + percent_FAL + percent_CAER, data=peak.pp, na.action = "na.fail")
m1.5 <- glmmTMB(sd ~ mean + dom_wild + percent_FAL, data=peak.pp, na.action = "na.fail")
m1.6 <- glmmTMB(sd ~ mean + dom_wild + percent_CAER, data=peak.pp, na.action = "na.fail")
m1.7 <- glmmTMB(sd ~ mean + dom_wild + percent_CAER, data=peak.pp, na.action = "na.fail")

AICctab(m1.4, m1.5, m1.6, m1.7) # m1.5 best, with Falcata

# No sig interactions; top models (within 2AIC) have mean, dom_wild
topmodel.pp <- glmmTMB(sd ~ mean + dom_wild + percent_FAL, data=peak.pp, na.action = "na.fail")
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.pp)) # pretty good, one outlier
check_collinearity(topmodel.pp)

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.pp, terms=c("dom_wild"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted")); ests1
ests <- spread(ests1, x, predicted); ests
names(ests) <- c("wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.pp)
ests$model_df.residual <- df.residual(topmodel.pp); ests

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
change.pp <- change.pp[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.pp

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Next, fit models with just domestication status #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# What is the change in variability with domestication, without accounting for effects of the mean?
##### A) Tweak ranFX structure NA
##### B) Tweak fixed FX structure NA
topmodel.pp_nomean <- glmmTMB(sd ~ dom_wild + percent_FAL, data=peak.pp, na.action = "na.fail")
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.pp_nomean))

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.pp_nomean, terms=c("dom_wild"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted")); ests1
ests <- spread(ests1, x, predicted); ests
names(ests) <- c("wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.pp_nomean)
ests$model_df.residual <- df.residual(topmodel.pp_nomean); ests

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
change.pp_nomean <- change.pp_nomean[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.pp_nomean

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#--#-#-#-#-#-#
# Next, fit models to show shift in mean # 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# What is the change in mean with domestication?
# THIS = THE COVARIATE IN GAMMA & WHOLE-PLANT BETA DIV MODELS

##### A) Tweak random FX structure
m3.1 <- glmmTMB(mean ~ dom_wild*plant.volume_centered + (1|repro.stage.simp), data=peak.pp, na.action = "na.fail")
m3.2 <- glmmTMB(mean ~ dom_wild*plant.volume_centered + (1|repro.stage.simp_dummy), data=peak.pp, na.action = "na.fail")
AICctab(m3.1, m3.2) # Within 2 AIC. So, having reproductive stage doesn't add anything to the model

##### B) Tweak fixed FX structure
# Check ancestry 
m3.3 <- glmmTMB(mean ~ dom_wild*plant.volume_centered + percent_FAL + percent_CAER, data=peak.pp, na.action = "na.fail")
dredge(m3.3) # Top models (by > 2AIC) include dom_wild + plant size +percent FAL
topmodel.mean <- glmmTMB(mean ~ dom_wild + plant.volume_centered + percent_FAL, data=peak.pp, na.action = "na.fail")
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.mean))

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.mean, terms=c("dom_wild"))); ests0
ggplot(ests0, aes(x=x, y=predicted)) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), size=.5, width=.2, color="black")
ests1 <- subset(ests0, select=c("x", "predicted")); ests1
ests <- spread(ests1, x, predicted); ests
names(ests) <- c("wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.mean)
ests$model_df.residual <- df.residual(topmodel.mean); ests

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
change.mean <- change.mean[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.mean

#--------------------------------------------------------------#
# Rbind all % changes for the whole-plant / across all leaves  #
#--------------------------------------------------------------#
estimates.pp <- rbind(change.pp, change.pp_nomean, change.mean); estimates.pp

#######################################
# 2) Model fitting: by leaf age class #
#######################################
head(peak.ps); str(peak.ps)

# transforming SD works better for some models; need a column already-transformed for bootstrapping
range(peak.ps$sd)
peak.ps$sd.log <- log(peak.ps$sd)
peak.ps$mean.log <- log(peak.ps$mean)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# First, fit models with both mean and domestication status #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
##### A) Tweak ranFX structure
peak.ps$plant_pop_dummy <- as.factor(sample(unique(peak.ps$plant_pop), replace=T, nrow(peak.ps))); head(peak.ps)
m5.1 <- glmmTMB(sd ~ dom_wild + leaf_stage + mean + plant.volume_centered + (1|plant_pop), data=peak.ps) 
m5.2 <- glmmTMB(sd ~ dom_wild + leaf_stage + mean + plant.volume_centered + (1|repro.stage.simp), data=peak.ps)
m5.1_dummy <- glmmTMB(sd ~ dom_wild + leaf_stage + mean + plant.volume_centered + (1|plant_pop_dummy), data=peak.ps)

AICctab(m5.1, m5.2, m5.1_dummy) # m5.1 slightly better; use plant_pop

##### B) Tweak fixedFX structure
m5.3 <- glmmTMB(sd ~ dom_wild*leaf_stage*mean + plant.volume_centered + percent_FAL + percent_CAER + (1|plant_pop), data=peak.ps)
dredge(m5.3) # top models (< 2 AIC) have all the main effects, no intxns

topmodel.ps <- glmmTMB(log(sd) ~ dom_wild + leaf_stage + mean + plant.volume_centered + percent_FAL + (1|plant_pop), data=peak.ps) 
topmodel.ps2 <- glmmTMB(sd.log ~ dom_wild + leaf_stage + mean + plant.volume_centered + percent_FAL+(1|plant_pop), data=peak.ps) # bootMer needs to map the response var to the original data (e.g. pre-transformed) 
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.ps)) # better with log transformation

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.ps, terms=c("dom_wild", "leaf_stage")))
ggplot(ests0, aes(x=x, y=predicted)) +
  facet_grid(~group) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), size=.5, width=.2, color="black")
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.ps)
ests$model_df.residual <- df.residual(topmodel.ps); ests

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
b1 <- lme4::bootMer(topmodel.ps2, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")
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
# This is a little different, because of the log-transformation: See here
# https://stats.stackexchange.com/questions/362556/linear-mixed-effect-model-interpretation-with-log-transformed-dependent-variable; also https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
# "In addition, you should interpret your transformed beta coefficients (and associated confidence intervals) as multiplicative changes once you transform to the original units. E.g., beta from the model = 0.0085; exp(.008585) ~= 1.0086 = 0.86% change"... "Only the dependent/response variable is log-transformed. Exponentiate the coefficient, subtract one from this number, and multiply by 100. This gives the percent increase (or decrease) in the response for every one-unit increase in the independent variable. Example: the coefficient is 0.198. (exp(0.198) – 1) * 100 = 21.9. For every one-unit increase in the independent variable, our dependent variable increases by about 22%."
change.ps <- cbind(pred_sig, ci); change.ps
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
change.ps <- change.ps[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.ps

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Next, fit models with just domestication status #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
##### A) Tweak ranFX structure NA
##### B) Tweak fixed FX structure NA
# What is the change in variability within each leaf age class with domestication, without accounting for the mean?
topmodel.ps_nomean <- glmmTMB(log(sd) ~ dom_wild*leaf_stage + plant.volume_centered + percent_FAL + (1|plant_pop), data=peak.ps) 
topmodel.ps_nomean2 <- glmmTMB(sd.log ~ dom_wild*leaf_stage + plant.volume_centered + percent_FAL + (1|plant_pop), data=peak.ps) 
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.ps_nomean))

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.ps_nomean, terms=c("dom_wild", "leaf_stage")))
ggplot(ests0, aes(x=x, y=predicted)) +
  facet_grid(~group) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), size=.5, width=.2, color="black")
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.ps_nomean)
ests$model_df.residual <- df.residual(topmodel.ps_nomean); ests

##### D) Calculate & append significance
sig0 <- data.frame(pairs(emmeans(topmodel.ps_nomean, ~ dom_wild|leaf_stage))) # post-hoc test 
names(sig0)[names(sig0) == 't.ratio'] <- 'statistic'
names(sig0)[names(sig0) == 'df'] <- 'statistic.df'
sig0$test.statistic <- "t.ratio"
sig <- subset(sig0, select=c("leaf_stage", "test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- merge(ests, sig, by="leaf_stage")

##### E) Bootstrap CI 
# Here, I'll extract the effect of domestication within each leaf stage
# Also used this helpful page to re-level the order of comparisons in emmeans (for some reason it always considers domestic the 'control' group, so estimates are negative rather than positive). 
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
# For example - 
normal <- pairs(emmeans(topmodel.ps_nomean, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.ps_nomean, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.ps_nomean2, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
summ <- data.frame(beta = b1$t0); summ
summ$leaf_stage <- rownames(summ)
conf <- data.frame(confint(b1, type = "perc")); conf
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
conf$leaf_stage <- rownames(conf)
ci <- merge(summ, conf, by="leaf_stage"); ci

##### F) Calculate % changes
# This is a little different, because of the log-transformation: See here
# https://stats.stackexchange.com/questions/362556/linear-mixed-effect-model-interpretation-with-log-transformed-dependent-variable; also https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
change.ps_nomean <- cbind(pred_sig, ci); change.ps_nomean
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
change.ps_nomean <- change.ps_nomean[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.ps_nomean

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#
# Finally, fit models to show shift in mean #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# What is the change in mean by leaf age class with domestication?
##### A) Tweak random FX structure
m6.1 <- glmmTMB(mean ~ dom_wild + leaf_stage + plant.volume_centered + (1|plant_pop) + (1|repro.stage.simp), data=peak.ps, na.action = "na.fail")
m6.2 <- glmmTMB(mean ~ dom_wild + leaf_stage + plant.volume_centered + (1|repro.stage.simp), data=peak.ps, na.action = "na.fail")
m6.3 <- glmmTMB(mean ~ dom_wild + leaf_stage + plant.volume_centered + (1|plant_pop), data=peak.ps, na.action = "na.fail")
AICctab(m6.1, m6.2, m6.3) # the model with just plant_pop best > 2AIC difference

##### B) Tweak fixed FX structure
m6.4 <- glmmTMB(mean ~ dom_wild*leaf_stage + plant.volume_centered + percent_FAL + percent_CAER + (1|plant_pop), data=peak.ps, na.action = "na.fail")
dredge(m6.4) # Top models have interaction + plant size
topmodel.mean <- glmmTMB(log(mean) ~ dom_wild*leaf_stage + plant.volume_centered + percent_FAL + (1|plant_pop), data=peak.ps, na.action = "na.fail")
topmodel.mean2 <- glmmTMB(mean.log ~ dom_wild*leaf_stage + plant.volume_centered + percent_FAL + (1|plant_pop), data=peak.ps, na.action = "na.fail")
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.mean))

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.mean, terms=c("dom_wild", "leaf_stage")))
ggplot(ests0, aes(x=x, y=predicted)) +
  facet_grid(~group) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), size=.5, width=.2, color="black")
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.mean)
ests$model_df.residual <- df.residual(topmodel.mean); ests

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
change.mean <- cbind(pred_sig, ci); change.mean
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
change.mean <- change.mean[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
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
# add trait
both$trait <- "Total saponins per leaf: avg"
# Rename
change.abu.all <- both; head(change.abu.all)

#######################################
#######################################
#         DIVERSITY ANALYSES          # 
#######################################
#######################################
head(diversity.pp) # total compounds covariate: total_peak_area_perLeaf_mean (among-leaf avg)
head(diversity.ps) # total compounds covariate: total_peak_area_perLeaf (per-leaf)

# Center the total compound amounts, untransformed
hist(diversity.pp$total_peak_area_perLeaf_mean)
diversity.pp$total_peak_area_perLeaf_mean.cent <- diversity.pp$total_peak_area_perLeaf_mean - mean(diversity.pp$total_peak_area_perLeaf_mean)
hist(diversity.pp$total_peak_area_perLeaf_mean.cent)

hist(diversity.ps$total_peak_area_perLeaf)
diversity.ps$total_peak_area_perLeaf.cent <- diversity.ps$total_peak_area_perLeaf - mean(diversity.ps$total_peak_area_perLeaf)
hist(diversity.ps$total_peak_area_perLeaf.cent)

# Also log-transform and then center; this really improves model Q-Q plots for richness 
diversity.pp$total_peak_area_perLeaf_mean.log <- log(diversity.pp$total_peak_area_perLeaf_mean)
diversity.pp$total_peak_area_perLeaf_mean.log.cent <- diversity.pp$total_peak_area_perLeaf_mean.log - mean(diversity.pp$total_peak_area_perLeaf_mean.log)
hist(diversity.pp$total_peak_area_perLeaf_mean.log.cent)

diversity.ps$total_peak_area_perLeaf.log <- log(diversity.ps$total_peak_area_perLeaf)
diversity.ps$total_peak_area_perLeaf.log.cent <- diversity.ps$total_peak_area_perLeaf.log - mean(diversity.ps$total_peak_area_perLeaf.log)
hist(diversity.ps$total_peak_area_perLeaf.log.cent)

#-#-#-#-#-#-#-#-#-#-#-#
#---------------------#
#      RICHNESS      #
#---------------------#
#-#-#-#-#-#-#-#-#-#-#-#
#------------------------------------#
# Abundance - richness relationships #
#------------------------------------#
# For fun, check out relationships between richness and  peak area
# Per plant
ggplot(data=diversity.pp, aes(x = total_peak_area_perLeaf_mean, y=richness)) +
  geom_point(size=2, alpha=0.5) +
  facet_grid(~dom_wild) + 
  geom_smooth(size=.5, se=FALSE, span=5) +
  theme(panel.background = element_rect(fill="white", color="black")) +
  theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
  #  theme_classic() +
  theme(legend.position = "none") +
  ylab("Richness of compounds in a plant\n(summed across N=9 leaves)") +
  xlab("Compound production\n(total compounds/leaf, avg across N=9 leaves/plant)") +
  theme(legend.position = "none") +
  theme(axis.title.x=element_text(size=13, face="bold", margin = ggplot2::margin(t = 10, r = 10, b = 0, l = 0)), axis.title.y=element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 10, vjust=.5), axis.text.y=element_text(size=10))

# Per leaf & age class
ggplot(data=diversity.ps, aes(x = total_peak_area_perLeaf, y=richness)) +
  geom_point(aes(color=leaf_stage), size=2, alpha=0.5) +
  scale_color_manual(values=c("young" = "olivedrab1", "middle" = "limegreen", "old" = "darkgreen")) +
  facet_grid(~dom_wild) + 
  geom_smooth(aes(group=leaf_stage, col=leaf_stage), size=.5, se=FALSE, span=5) +
  geom_smooth(aes(group=dom_wild), size=.8, col="gray10", se=FALSE, span=5) +
  theme(panel.background = element_rect(fill="white", color="black")) +
  theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
  #  theme_classic() +
  theme(legend.position = "none") +
  ylab("Richness of compounds in a leaf") +
  xlab("Total peak area in a leaf") +
  theme(legend.position = "none") +
  theme(axis.title.x=element_text(size=13, face="bold", margin = ggplot2::margin(t = 10, r = 10, b = 0, l = 0)), axis.title.y=element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 10, vjust=.5), axis.text.y=element_text(size=10))
# These suggest that the slope of the peak area - richness relationship might vary within leaf age group... should try to fit models that allow this

# # # # # # # # # # # # # # # 
#        Whole-plant        #
# # # # # # # # # # # # # # # 
#---------------------------------------#
#  WITHOUT average peak area/leaf/plant #
#---------------------------------------#
head(diversity.pp); hist(diversity.pp$richness) # actually quite normal! Bodes well!
m7 <- glmmTMB(richness ~ dom_wild + plant.volume_centered + repro.stage.simp+ percent_FAL + percent_CAER, data=diversity.pp)
dredge(m7)
topmodel.rich.pp <- glmmTMB(richness ~ dom_wild + plant.volume_centered + percent_FAL, data=diversity.pp)
summary(topmodel.rich.pp)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.rich.pp)) # all very good

##### Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.rich.pp, terms=c("dom_wild"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted")); ests1
ests <- spread(ests1, x, predicted); ests
names(ests) <- c("wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.rich.pp)
ests$model_df.residual <- df.residual(topmodel.rich.pp); ests

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.rich.pp, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- cbind(ests, sig); pred_sig

##### E) Bootstrap CI 
# This is from the March 15, 2020 glmmTMB vignette. See "isLMM.glmmTMB" section. Instead of "$zi", in their example, I am interested in the cond estimate for domestic - so I name that specifically within the function
fixef(topmodel.rich.pp)
fixef(topmodel.rich.pp)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object

b1 <- lme4::bootMer(topmodel.rich.pp, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")

summ <- data.frame(beta = b1$t0); summ
conf <- data.frame(confint(b1, type = "perc")); conf
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci <- merge(summ, conf); ci

##### F) Calculate % changes
change.rich.pp <- cbind(pred_sig, ci); change.rich.pp
change.rich.pp$percent.change <- change.rich.pp$beta/change.rich.pp$wild.pred*100
change.rich.pp$percent.change.lower <- change.rich.pp$beta.CI_lower.boot/change.rich.pp$wild.pred*100
change.rich.pp$percent.change.upper <- change.rich.pp$beta.CI_upper.boot/change.rich.pp$wild.pred*100

change.rich.pp$leaf_stage <- "all"
change.rich.pp$response <- "richness"
change.rich.pp$scale <- "per-plant"
change.rich.pp$model <- "status_only"

##### G) Put it all together & re-order for merging
change.rich.pp <- change.rich.pp[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.rich.pp

#------------------------------------#
#  WITH average peak area/leaf/plant #
#------------------------------------#
topmodel.rich.pp_mean <- glmmTMB(richness ~ dom_wild + total_peak_area_perLeaf_mean.log.cent + plant.volume_centered  + percent_FAL, data=diversity.pp)
summary(topmodel.rich.pp_mean)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.rich.pp_mean)) # Also all diagnostics good

##### Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.rich.pp_mean, terms=c("dom_wild"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted")); ests1
ests <- spread(ests1, x, predicted); ests
names(ests) <- c("wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.rich.pp_mean)
ests$model_df.residual <- df.residual(topmodel.rich.pp_mean); ests

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.rich.pp_mean, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- cbind(ests, sig); pred_sig

##### E) Bootstrap CI 
# This is from the March 15, 2020 glmmTMB vignette. See "isLMM.glmmTMB" section. Instead of "$zi", in their example, I am interested in the cond estimate for domestic - so I name that specifically within the function
fixef(topmodel.rich.pp_mean)
fixef(topmodel.rich.pp_mean)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object

b1 <- lme4::bootMer(topmodel.rich.pp_mean, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")

summ <- data.frame(beta = b1$t0); summ
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci <- merge(summ, conf); ci

##### F) Calculate % changes
change.rich.pp_mean <- cbind(pred_sig, ci); change.rich.pp_mean
change.rich.pp_mean$percent.change <- change.rich.pp_mean$beta/change.rich.pp_mean$wild.pred*100
change.rich.pp_mean$percent.change.lower <- change.rich.pp_mean$beta.CI_lower.boot/change.rich.pp_mean$wild.pred*100
change.rich.pp_mean$percent.change.upper <- change.rich.pp_mean$beta.CI_upper.boot/change.rich.pp_mean$wild.pred*100

change.rich.pp_mean$leaf_stage <- "all"
change.rich.pp_mean$response <- "richness"
change.rich.pp_mean$scale <- "per-plant"
change.rich.pp_mean$model <- "status+meanpL"

##### G) Put it all together & re-order for merging
change.rich.pp_mean <- change.rich.pp_mean[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.rich.pp_mean

# # # # # # # # # # # # # # # 
#        Per-leaf age       #
# # # # # # # # # # # # # # # 
#-----------------------------#
#  WITHOUT peak area per leaf #
#-----------------------------#
head(diversity.ps)
# Manipulate random FX structure - nothing to manipulate
rich1 <- glmmTMB(richness ~ dom_wild*leaf_stage + plant.volume_centered + repro.stage.simp, data=diversity.ps)
rich2 <- glmmTMB(richness ~ dom_wild + leaf_stage + plant.volume_centered + repro.stage.simp, data=diversity.ps)
AICctab(rich1, rich2) # rich2 is favored; no interaction
# Ancestry improve the model?
rich3 <- glmmTMB(richness ~ dom_wild + leaf_stage + plant.volume_centered + repro.stage.simp + percent_FAL + percent_CAER, data=diversity.ps)
dredge(rich3)
topmodel.rich.ps <- glmmTMB(richness ~ dom_wild + leaf_stage + percent_FAL + percent_CAER, data=diversity.ps)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.rich.ps)) # good

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.rich.ps, terms=c("dom_wild", "leaf_stage")))
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.rich.ps)
ests$model_df.residual <- df.residual(topmodel.rich.ps); ests

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.rich.ps, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- cbind(ests, sig); pred_sig

##### E) Bootstrap CI 
# Here, I'll extract the effect of domestication within each leaf stage, due to the significant interaction 
# (Used this great post as a resource: https://www.r-bloggers.com/bootstrapping-follow-up-contrasts-for-within-subject-anovas-part-2/)
# Also used this helpful page to re-level the order of comparisons in emmeans (for some reason it always considers domestic the 'control' group, so estimates are negative rather than positive). 
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
# For example - 
normal <- pairs(emmeans(topmodel.rich.ps, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.rich.ps, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.rich.ps, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
summ <- data.frame(beta = b1$t0)
summ$leaf_stage <- rownames(summ)
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
conf$leaf_stage <- rownames(conf)
ci <- merge(summ, conf, by="leaf_stage"); ci

##### F) Calculate % changes
change.rich.ps <- merge(pred_sig, ci, by="leaf_stage"); change.rich.ps
change.rich.ps$percent.change <- change.rich.ps$beta/change.rich.ps$wild.pred*100
change.rich.ps$percent.change.lower <- change.rich.ps$beta.CI_lower.boot/change.rich.ps$wild.pred*100
change.rich.ps$percent.change.upper <- change.rich.ps$beta.CI_upper.boot/change.rich.ps$wild.pred*100

change.rich.ps$response <- "richness"
change.rich.ps$scale <- "per-age-class"
change.rich.ps$model <- "status_only"

##### G) Put it all together & re-order for merging
change.rich.ps <- change.rich.ps[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.rich.ps

#--------------------------#
#  WITH peak area per leaf #
#--------------------------#
head(diversity.ps)
hist(diversity.ps$richness) # quite nice

# Manipulate random FX structure
rich1 <- glmmTMB(richness ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|dom_wild) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps) # no convergence
rich2 <- glmmTMB(richness ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|dom_wild) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps) # no convergence
rich3 <- glmmTMB(richness ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|dom_wild), data=diversity.ps) # no convergence
rich4 <- glmmTMB(richness ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
rich5 <- glmmTMB(richness ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|dom_wild), data=diversity.ps)
rich6 <- glmmTMB(richness ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
AICctab(rich4, rich6) # rich 4 best of models that converge; area-richness relationship varies within leaf age class and within plant individual!!!!!
DHARMa::testResiduals(DHARMa::simulateResiduals(rich4))

# manipulate fixed FX structure
rich4.1 <- glmmTMB(richness ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
rich4.2 <- glmmTMB(richness ~ dom_wild*leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps) 
rich4.3 <- glmmTMB(richness ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps) # failed to converge
AICctab(rich4.1, rich4.2) # rich 4.2 best; 2-way interaction
rich4.4 <- glmmTMB(richness ~ dom_wild*leaf_stage + total_peak_area_perLeaf.log.cent + plant.volume_centered + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
rich4.5 <- glmmTMB(richness ~ dom_wild*leaf_stage + total_peak_area_perLeaf.log.cent + plant.volume_centered + repro.stage.simp + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
AICctab(rich4.1, rich4.2, rich4.4, rich4.5) # rich 4.2 still best

# Check if adding ancestry improves model
rich4.6 <-  glmmTMB(richness ~ dom_wild*leaf_stage + total_peak_area_perLeaf.log.cent + percent_FAL + percent_CAER + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps) 
rich4.7 <-  glmmTMB(richness ~ dom_wild*leaf_stage + total_peak_area_perLeaf.log.cent + percent_FAL + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps) 
rich4.8 <-  glmmTMB(richness ~ dom_wild*leaf_stage + total_peak_area_perLeaf.log.cent + percent_CAER + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps) # doesn'tconverge
AICctab(rich4.2, rich4.6, rich4.7) # rich 4.2 best

# This is the top model, with overall peak area in it as well
topmodel.rich.ps_total <-  glmmTMB(richness ~ dom_wild*leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps) 
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.rich.ps_total)) # OK
summary(topmodel.rich.ps_total) # Large positive effect of overall peak area on compound richness, which makes sense (more overall = more compounds). 

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.rich.ps_total, terms=c("dom_wild", "leaf_stage")))
ests <- spread(ests0, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.rich.ps_total)
ests$model_df.residual <- df.residual(topmodel.rich.ps_total); ests

##### D) Calculate / append significance 
sig0 <- data.frame(pairs(emmeans(topmodel.rich.ps_total, ~ dom_wild|leaf_stage))) # post-hoc test 
names(sig0)[names(sig0) == 't.ratio'] <- 'statistic'
names(sig0)[names(sig0) == 'df'] <- 'statistic.df'
sig0$test.statistic <- "t.ratio"
sig <- subset(sig0, select=c("leaf_stage", "test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- merge(ests, sig, by="leaf_stage"); pred_sig

##### E) Bootstrap CI 
# Here, I'll extract the effect of domestication within each leaf stage, due to the significant interaction 
# (Used this great post as a resource: https://www.r-bloggers.com/bootstrapping-follow-up-contrasts-for-within-subject-anovas-part-2/)
# Also used this helpful page to re-level the order of comparisons in emmeans (for some reason it always considers domestic the 'control' group, so estimates are negative rather than positive). 
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
# For example - 
normal <- pairs(emmeans(topmodel.rich.ps_total, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.rich.ps_total, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.rich.ps_total, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
summ <- data.frame(beta = b1$t0) # Works better on Will's computer
summ$leaf_stage <- rownames(summ)
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
conf$leaf_stage <- rownames(conf)
ci <- merge(summ, conf, by="leaf_stage"); ci

##### F) Calculate % changes
change.rich.ps_total <- merge(pred_sig, ci, by="leaf_stage"); change.rich.ps_total
change.rich.ps_total$percent.change <- change.rich.ps_total$beta/change.rich.ps_total$wild.pred*100
change.rich.ps_total$percent.change.lower <- change.rich.ps_total$beta.CI_lower.boot/change.rich.ps_total$wild.pred*100
change.rich.ps_total$percent.change.upper <- change.rich.ps_total$beta.CI_upper.boot/change.rich.ps_total$wild.pred*100

change.rich.ps_total$response <- "richness"
change.rich.ps_total$scale <- "per-age-class"
change.rich.ps_total$model <- "status+totalpL"

##### G) Put it all together & re-order for merging
change.rich.ps_total <- change.rich.ps_total[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.rich.ps_total

#---------------------------------------------------------------------------------------------#
# Rbind the two scales + w & w/o measures of mean or individual per leaf compound production #
#---------------------------------------------------------------------------------------------#
change.rich.all <- rbind(change.rich.pp, change.rich.pp_mean, change.rich.ps, change.rich.ps_total)
change.rich.all

# Name trait
change.rich.all$trait <- "richness"

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-----------------------------#
#      SHANNON DIVERSITY      #
#-----------------------------#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# NOTE: Running analyses on Shannon diversity calculated from raw peak areas vs sqrt-transformed peak areas makes a difference!
#------------------------------------#
# Abundance - Shannon relationships #
#------------------------------------#
# For fun, check out relationships between Shannon diversity and total peak area
# Per plant
ggplot(data=diversity.pp, aes(x = total_peak_area_perLeaf_mean, y=H)) + 
  geom_point(size=2, alpha=0.5) +
  facet_grid(~dom_wild) + 
  geom_smooth(size=.5, se=FALSE, span=5) +
  theme(panel.background = element_rect(fill="white", color="black")) +
  theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
  #  theme_classic() +
  theme(legend.position = "none") +
  ylab("Shannon diversity of compounds in a plant\n(summed across N=9 leaves)") +
  xlab("Compound production\n(total compounds/leaf, avg across N=9 leaves/plant)") +
  theme(legend.position = "none") +
  theme(axis.title.x=element_text(size=13, face="bold", margin = ggplot2::margin(t = 10, r = 10, b = 0, l = 0)), axis.title.y=element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 10, vjust=.5), axis.text.y=element_text(size=10))

# Per leaf & age class
ggplot(data=diversity.ps, aes(x = total_peak_area_perLeaf, y=H)) +
  geom_point(aes(color=leaf_stage), size=2, alpha=0.5) +
  scale_color_manual(values=c("young" = "olivedrab1", "middle" = "limegreen", "old" = "darkgreen")) +
  facet_grid(~dom_wild) + 
  geom_smooth(aes(group=leaf_stage, col=leaf_stage), size=.5, se=FALSE, span=5) +
  geom_smooth(aes(group=dom_wild), size=.8, col="gray10", se=FALSE, span=5) +
  theme(panel.background = element_rect(fill="white", color="black")) +
  theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
  #  theme_classic() +
  theme(legend.position = "none") +
  ylab("Shannon Diversity of compounds in a leaf") +
  xlab("Total peak area in a leaf") +
  theme(legend.position = "none") +
  theme(axis.title.x=element_text(size=13, face="bold", margin = ggplot2::margin(t = 10, r = 10, b = 0, l = 0)), axis.title.y=element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 10, vjust=.5), axis.text.y=element_text(size=10))
# These suggest that the slope of the peak area - Shannon diversity relationship might vary within leaf age group... should try to fit models that allow this

# # # # # # # # # # # # # # # 
#        Whole-plant        # 
# # # # # # # # # # # # # # # 
#---------------------------------------#
#  WITHOUT average peak area/leaf/plant #
#---------------------------------------#
head(diversity.pp); hist(diversity.pp$H)  # actually quite normal! Bodes well!

m1 <- glmmTMB(H ~ dom_wild + plant.volume_centered + repro.stage.simp, data=diversity.pp)
m2 <- glmmTMB(H ~ dom_wild + plant.volume_centered, data=diversity.pp)
m3 <- glmmTMB(H ~ dom_wild + repro.stage.simp, data=diversity.pp)
AICctab(m1, m2, m3) # m2 best

m4 <- glmmTMB(H ~ dom_wild + plant.volume_centered + percent_FAL + percent_CAER, data=diversity.pp) 
m5 <- glmmTMB(H ~ dom_wild + plant.volume_centered + percent_FAL, data=diversity.pp) 
m6 <- glmmTMB(H ~ dom_wild + plant.volume_centered + percent_CAER, data=diversity.pp) 
AICctab(m2, m4, m5, m6) # m6 best 
topmodel.H.pp <- glmmTMB(H ~ dom_wild + plant.volume_centered + percent_CAER, data=diversity.pp) 
summary(topmodel.H.pp)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.H.pp)) # looks very good

##### Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.H.pp, terms=c("dom_wild"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted")); ests1
ests <- spread(ests1, x, predicted); ests
names(ests) <- c("wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.H.pp)
ests$model_df.residual <- df.residual(topmodel.H.pp); ests

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.H.pp, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- cbind(ests, sig); pred_sig

##### E) Bootstrap CI 
# This is from the March 15, 2020 glmmTMB vignette. See "isLMM.glmmTMB" section. Instead of "$zi", in their example, I am interested in the cond estimate for domestic - so I name that specifically within the function
fixef(topmodel.H.pp)
fixef(topmodel.H.pp)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object

b1 <- lme4::bootMer(topmodel.H.pp, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")

summ <- data.frame(beta = b1$t0); summ
conf <- data.frame(confint(b1, type = "perc")); conf
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci <- merge(summ, conf); ci

##### F) Calculate % changes
change.H.pp <- cbind(pred_sig, ci); change.H.pp
change.H.pp$percent.change <- change.H.pp$beta/change.H.pp$wild.pred*100
change.H.pp$percent.change.lower <- change.H.pp$beta.CI_lower.boot/change.H.pp$wild.pred*100
change.H.pp$percent.change.upper <- change.H.pp$beta.CI_upper.boot/change.H.pp$wild.pred*100

change.H.pp$leaf_stage <- "all"
change.H.pp$response <- "shannon"
change.H.pp$scale <- "per-plant"
change.H.pp$model <- "status_only"

##### G) Put it all together & re-order for merging
change.H.pp <- change.H.pp[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.H.pp

#------------------------------------#
#  WITH average peak area/leaf/plant #
#------------------------------------#
topmodel.H.pp_mean <- glmmTMB(H ~ dom_wild + total_peak_area_perLeaf_mean.cent + plant.volume_centered + percent_CAER, data=diversity.pp)
summary(topmodel.H.pp_mean)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.H.pp_mean)) # looks good

##### Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.H.pp_mean, terms=c("dom_wild"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted")); ests1
ests <- spread(ests1, x, predicted); ests
names(ests) <- c("wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.H.pp_mean)
ests$model_df.residual <- df.residual(topmodel.H.pp_mean); ests

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.H.pp_mean, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- cbind(ests, sig); pred_sig

##### E) Bootstrap CI 
# This is from the March 15, 2020 glmmTMB vignette. See "isLMM.glmmTMB" section. Instead of "$zi", in their example, I am interested in the cond estimate for domestic - so I name that specifically within the function
fixef(topmodel.H.pp_mean)
fixef(topmodel.H.pp_mean)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object

b1 <- lme4::bootMer(topmodel.H.pp_mean, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")

summ <- data.frame(beta = b1$t0); summ
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci <- merge(summ, conf); ci

##### F) Calculate % changes
change.H.pp_mean <- cbind(pred_sig, ci); change.H.pp_mean
change.H.pp_mean$percent.change <- change.H.pp_mean$beta/change.H.pp_mean$wild.pred*100
change.H.pp_mean$percent.change.lower <- change.H.pp_mean$beta.CI_lower.boot/change.H.pp_mean$wild.pred*100
change.H.pp_mean$percent.change.upper <- change.H.pp_mean$beta.CI_upper.boot/change.H.pp_mean$wild.pred*100

change.H.pp_mean$leaf_stage <- "all"
change.H.pp_mean$response <- "shannon"
change.H.pp_mean$scale <- "per-plant"
change.H.pp_mean$model <- "status+meanpL"

##### G) Put it all together & re-order for merging
change.H.pp_mean <- change.H.pp_mean[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.H.pp_mean

# # # # # # # # # # # # # # # 
#        Per-leaf age       # 
# # # # # # # # # # # # # # # 
#-----------------------------#
#  WITHOUT peak area per leaf # 
#-----------------------------#
head(diversity.ps)
# Manipulate random FX structure - nothing to manipulate
H1 <- glmmTMB(H ~ dom_wild + leaf_stage + plant.volume_centered + repro.stage.simp, data=diversity.ps); DHARMa::testResiduals(DHARMa::simulateResiduals(H1))
H2 <- glmmTMB(H ~ dom_wild + leaf_stage + plant.volume_centered, data=diversity.ps); DHARMa::testResiduals(DHARMa::simulateResiduals(H2))
H3 <- glmmTMB(H ~ dom_wild + leaf_stage + repro.stage.simp, data=diversity.ps); DHARMa::testResiduals(DHARMa::simulateResiduals(H3))
AICctab(H1, H2, H3) # H1 best, + best diagnnostics

H4 <- glmmTMB(H ~ dom_wild*leaf_stage + plant.volume_centered + repro.stage.simp, data=diversity.ps); DHARMa::testResiduals(DHARMa::simulateResiduals(H4))
AICctab(H1, H4) # H1 better

# Does ancestry improve the model
H5 <- glmmTMB(H ~ dom_wild + leaf_stage + plant.volume_centered + repro.stage.simp + percent_FAL + percent_CAER, data=diversity.ps)
H6 <- glmmTMB(H ~ dom_wild + leaf_stage + plant.volume_centered + repro.stage.simp + percent_FAL, data=diversity.ps)
H7 <- glmmTMB(H ~ dom_wild + leaf_stage + plant.volume_centered + repro.stage.simp + percent_CAER, data=diversity.ps)
AICctab(H2, H5, H6, H7)

topmodel.H.ps <- glmmTMB(H ~ dom_wild + leaf_stage + plant.volume_centered + percent_CAER, data=diversity.ps); summary(topmodel.H.ps)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.H.ps)) # Q-Q not perfect, but not a lot I can do (log transformation doesn't help). Interactive model actually has a better QQ, despite being the lower AIC, interestingly

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.H.ps, terms=c("dom_wild", "leaf_stage")))
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.H.ps)
ests$model_df.residual <- df.residual(topmodel.H.ps); ests

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.H.ps, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- merge(ests, sig); pred_sig # no interaction, so will be same sig across

##### E) Bootstrap CI 
# Here, I'll extract the effect of domestication within each leaf stage, due to the significant interaction 
# (Used this great post as a resource: https://www.r-bloggers.com/bootstrapping-follow-up-contrasts-for-within-subject-anovas-part-2/)
# Also used this helpful page to re-level the order of comparisons in emmeans (for some reason it always considers domestic the 'control' group, so estimates are negative rather than positive). 
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
# For example - 
normal <- pairs(emmeans(topmodel.H.ps, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.H.ps, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.H.ps, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
summ <- data.frame(beta = b1$t0)
summ$leaf_stage <- rownames(summ)
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
conf$leaf_stage <- rownames(conf)
ci <- merge(summ, conf, by="leaf_stage"); ci

##### F) Calculate % changes
change.H.ps <- merge(pred_sig, ci, by="leaf_stage"); change.H.ps
change.H.ps$percent.change <- change.H.ps$beta/change.H.ps$wild.pred*100
change.H.ps$percent.change.lower <- change.H.ps$beta.CI_lower.boot/change.H.ps$wild.pred*100
change.H.ps$percent.change.upper <- change.H.ps$beta.CI_upper.boot/change.H.ps$wild.pred*100

change.H.ps$response <- "shannon"
change.H.ps$scale <- "per-age-class"
change.H.ps$model <- "status_only"

##### G) Put it all together & re-order for merging
change.H.ps <- change.H.ps[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.H.ps

#--------------------------#
#  WITH peak area per leaf #
#--------------------------#
# Manipulate random FX structure
H1 <- glmmTMB(H ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|dom_wild) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps) # doesn't converge
H2 <- glmmTMB(H ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|dom_wild) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
H3 <- glmmTMB(H ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|dom_wild), data=diversity.ps)
H4 <- glmmTMB(H ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
H5 <- glmmTMB(H ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|dom_wild), data=diversity.ps)
H6 <- glmmTMB(H ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
# only 6 converges; area-diversity relationship varies w/in plant individual

# manipulate fixed FX structure
H6.1 <- glmmTMB(H ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
H6.2 <- glmmTMB(H ~ dom_wild*leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
H6.3 <- glmmTMB(H ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps) #
AICctab(H6.1, H6.2, H6.3) # H6.3 best; 3-way interaction
H6.4 <- glmmTMB(H ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + plant.volume_centered + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
H6.5 <- glmmTMB(H ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + plant.volume_centered + repro.stage.simp + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
AICctab(H6.3, H6.4, H6.5) # H6.3 still best

# Add ancestry
H6.6 <- glmmTMB(H ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + percent_FAL + percent_CAER + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
H6.7 <- glmmTMB(H ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + percent_FAL + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
H6.8 <- glmmTMB(H ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + percent_CAER + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
AICctab(H6.3, H6.6, H6.7, H6.8)

DHARMa::testResiduals(DHARMa::simulateResiduals(H6.3)) # not perfect, but not a lot to tweak

# This is the top model
topmodel.H.ps_total <-  glmmTMB(H ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + percent_CAER + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
summary(topmodel.H.ps_total) 

##### C) Pull out predicted means for wild & domestic + model nobs + df
mean(diversity.ps$total_peak_area_perLeaf.log.cent)
ests0 <- data.frame(ggpredict(topmodel.H.ps_total, terms=c("dom_wild", "leaf_stage", "total_peak_area_perLeaf.log.cent [-1.997946e-16]")))
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.H.ps_total)
ests$model_df.residual <- df.residual(topmodel.H.ps_total); ests

##### D) Calculate / append significance 
sig0 <- data.frame(pairs(emmeans(topmodel.H.ps_total, ~ dom_wild|leaf_stage*total_peak_area_perLeaf.log.cent))) # post-hoc test 
names(sig0)[names(sig0) == 't.ratio'] <- 'statistic'
names(sig0)[names(sig0) == 'df'] <- 'statistic.df'
sig0$test.statistic <- "t.ratio"
sig <- subset(sig0, select=c("leaf_stage", "test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- merge(ests, sig, by="leaf_stage"); pred_sig

##### E) Bootstrap CI 
# Here, I'll extract the effect of domestication within each leaf stage, due to the significant interaction 
# (Used this great post as a resource: https://www.r-bloggers.com/bootstrapping-follow-up-contrasts-for-within-subject-anovas-part-2/)
# Also used this helpful page to re-level the order of comparisons in emmeans (for some reason it always considers domestic the 'control' group, so estimates are negative rather than positive). 
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
# For example - 
normal <- pairs(emmeans(topmodel.H.ps_total, ~ dom_wild|leaf_stage*total_peak_area_perLeaf.log.cent)); normal # estimates reversed
#versus
em <- emmeans(topmodel.H.ps_total, specs = trt.vs.ctrl ~ dom_wild|leaf_stage*total_peak_area_perLeaf.log.cent); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage*total_peak_area_perLeaf.log.cent)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.H.ps_total, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
summ <- data.frame(beta = b1$t0) # Works better on Will's computer
summ$leaf_stage <- rownames(summ)
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
conf$leaf_stage <- rownames(conf)
ci <- merge(summ, conf, by="leaf_stage"); ci

##### F) Calculate % changes
change.H.ps_total <- merge(pred_sig, ci, by="leaf_stage"); change.H.ps_total
change.H.ps_total$percent.change <- change.H.ps_total$beta/change.H.ps_total$wild.pred*100
change.H.ps_total$percent.change.lower <- change.H.ps_total$beta.CI_lower.boot/change.H.ps_total$wild.pred*100
change.H.ps_total$percent.change.upper <- change.H.ps_total$beta.CI_upper.boot/change.H.ps_total$wild.pred*100

change.H.ps_total$response <- "shannon"
change.H.ps_total$scale <- "per-age-class"
change.H.ps_total$model <- "status+totalpL"

##### G) Put it all together & re-order for merging
change.H.ps_total <- change.H.ps_total[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.H.ps_total

#---------------------------------------------------------#
# Rbind the two scales + w & w/o total peak area together #
#---------------------------------------------------------#
change.H.all <- rbind(change.H.pp, change.H.pp_mean, change.H.ps, change.H.ps_total)
change.H.all 

# Name trait
change.H.all$trait <- "shannon"

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-----------------------------#
#       SIMPSON DIVERSITY     # 
#-----------------------------#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#------------------------------------#
# Abundance - Simpson relationships  #
#------------------------------------#
head(diversity.pp); head(diversity.ps) # 
# Per plant
ggplot(data=diversity.pp, aes(x = total_peak_area_perLeaf_mean, y=simp)) + 
  geom_point(size=2, alpha=0.5) +
  facet_grid(~dom_wild) + 
  geom_smooth(size=.5, se=FALSE, span=5) +
  theme(panel.background = element_rect(fill="white", color="black")) +
  theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
  #  theme_classic() +
  theme(legend.position = "none") +
  ylab("Simpson diversity of compounds in a plant\n(summed across N=9 leaves)") +
  xlab("Total peak area in a plant\n(summed across N=9 leaves)") +
  theme(legend.position = "none") +
  theme(axis.title.x=element_text(size=13, face="bold", margin = ggplot2::margin(t = 10, r = 10, b = 0, l = 0)), axis.title.y=element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 10, vjust=.5), axis.text.y=element_text(size=10))

# Per leaf & age class
ggplot(data=diversity.ps, aes(x = total_peak_area_perLeaf, y=simp)) +
  geom_point(aes(color=leaf_stage), size=2, alpha=0.5) +
  scale_color_manual(values=c("young" = "olivedrab1", "middle" = "limegreen", "old" = "darkgreen")) +
  facet_grid(~dom_wild) + 
  geom_smooth(aes(group=leaf_stage, col=leaf_stage), size=.5, se=FALSE, span=5) +
  geom_smooth(aes(group=dom_wild), size=.8, col="gray10", se=FALSE, span=5) +
  theme(panel.background = element_rect(fill="white", color="black")) +
  theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
  #  theme_classic() +
  theme(legend.position = "none") +
  ylab("Simpson Diversity of compounds in a leaf") +
  xlab("Total peak area in a leaf") +
  theme(legend.position = "none") +
  theme(axis.title.x=element_text(size=13, face="bold", margin = ggplot2::margin(t = 10, r = 10, b = 0, l = 0)), axis.title.y=element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 10, vjust=.5), axis.text.y=element_text(size=10))
# These suggest that the slope of the peak area - Shannon diversity relationship might vary within leaf age group... should try to fit models that allow this

# # # # # # # # # # # # # # # 
#        Whole-plant        # 
# # # # # # # # # # # # # # # 
#---------------------------------------#
#  WITHOUT average peak area/leaf/plant #
#---------------------------------------#
head(diversity.pp); hist(diversity.pp$simp); hist(log(diversity.pp$simp))  # hmm not great, log increases skew

m1 <- glmmTMB(simp ~ dom_wild + plant.volume_centered + repro.stage.simp, data=diversity.pp)
m2 <- glmmTMB(simp ~ dom_wild + plant.volume_centered, data=diversity.pp)
m3 <- glmmTMB(simp ~ dom_wild + repro.stage.simp, data=diversity.pp)
AICctab(m1, m2, m3) # m2 a bit better
AICctab(m2, m4, m5, m6)

# Does adding ancestry improve the model
m4 <- glmmTMB(simp ~ dom_wild + plant.volume_centered + percent_FAL + percent_CAER, data=diversity.pp)
m5 <- glmmTMB(simp ~ dom_wild + plant.volume_centered + percent_FAL, data=diversity.pp)
m6 <- glmmTMB(simp ~ dom_wild + plant.volume_centered + percent_CAER, data=diversity.pp)
AICctab(m2, m4, m5, m6)

topmodel.simp.pp <- glmmTMB(simp ~ dom_wild + plant.volume_centered + percent_CAER, data=diversity.pp)
summary(topmodel.simp.pp)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.simp.pp)) # look fine

##### Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.simp.pp, terms=c("dom_wild"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted")); ests1
ests <- spread(ests1, x, predicted); ests
names(ests) <- c("wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.simp.pp)
ests$model_df.residual <- df.residual(topmodel.simp.pp); ests

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.simp.pp, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- cbind(ests, sig); pred_sig

##### E) Bootstrap CI 
# This is from the March 15, 2020 glmmTMB vignette. See "isLMM.glmmTMB" section. Instead of "$zi", in their example, I am interested in the cond estimate for domestic - so I name that specifically within the function
fixef(topmodel.simp.pp)
fixef(topmodel.simp.pp)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object

b1 <- lme4::bootMer(topmodel.simp.pp, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")

summ <- data.frame(beta = b1$t0); summ
conf <- data.frame(confint(b1, type = "perc")); conf
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci <- merge(summ, conf); ci

##### F) Calculate % changes
change.simp.pp <- cbind(pred_sig, ci); change.simp.pp
change.simp.pp$percent.change <- change.simp.pp$beta/change.simp.pp$wild.pred*100
change.simp.pp$percent.change.lower <- change.simp.pp$beta.CI_lower.boot/change.simp.pp$wild.pred*100
change.simp.pp$percent.change.upper <- change.simp.pp$beta.CI_upper.boot/change.simp.pp$wild.pred*100

change.simp.pp$leaf_stage <- "all"
change.simp.pp$response <- "simpson"
change.simp.pp$scale <- "per-plant"
change.simp.pp$model <- "status_only"

##### G) Put it all together & re-order for merging
change.simp.pp <- change.simp.pp[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.simp.pp

#------------------------------------#
#  WITH average peak area/leaf/plant #
#------------------------------------#
topmodel.simp.pp_mean <- glmmTMB(simp ~ dom_wild + total_peak_area_perLeaf_mean.cent + plant.volume_centered + percent_CAER, data=diversity.pp)
summary(topmodel.simp.pp_mean)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.simp.pp_mean)) # looks good

##### Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.simp.pp_mean, terms=c("dom_wild"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted")); ests1
ests <- spread(ests1, x, predicted); ests
names(ests) <- c("wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.simp.pp_mean)
ests$model_df.residual <- df.residual(topmodel.simp.pp_mean); ests

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.simp.pp_mean, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- cbind(ests, sig); pred_sig

##### E) Bootstrap CI 
# This is from the March 15, 2020 glmmTMB vignette. See "isLMM.glmmTMB" section. Instead of "$zi", in their example, I am interested in the cond estimate for domestic - so I name that specifically within the function
fixef(topmodel.simp.pp_mean)
fixef(topmodel.simp.pp_mean)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object

b1 <- lme4::bootMer(topmodel.simp.pp_mean, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")

summ <- data.frame(beta = b1$t0); summ
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci <- merge(summ, conf); ci

##### F) Calculate % changes
change.simp.pp_mean <- cbind(pred_sig, ci); change.simp.pp_mean
change.simp.pp_mean$percent.change <- change.simp.pp_mean$beta/change.simp.pp_mean$wild.pred*100
change.simp.pp_mean$percent.change.lower <- change.simp.pp_mean$beta.CI_lower.boot/change.simp.pp_mean$wild.pred*100
change.simp.pp_mean$percent.change.upper <- change.simp.pp_mean$beta.CI_upper.boot/change.simp.pp_mean$wild.pred*100

change.simp.pp_mean$leaf_stage <- "all"
change.simp.pp_mean$response <- "simpson"
change.simp.pp_mean$scale <- "per-plant"
change.simp.pp_mean$model <- "status+meanpL"

##### G) Put it all together & re-order for merging
change.simp.pp_mean <- change.simp.pp_mean[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.simp.pp_mean

# # # # # # # # # # # # # # # 
#        Per-leaf age       # 
# # # # # # # # # # # # # # # 
#-----------------------------#
#  WITHOUT peak area per leaf #
#-----------------------------#
head(diversity.ps)

m1 <- glmmTMB(simp ~ dom_wild + leaf_stage + plant.volume_centered + repro.stage.simp, data=diversity.ps)
m2 <- glmmTMB(simp ~ dom_wild + plant.volume_centered + repro.stage.simp, data=diversity.ps)
m3 <- glmmTMB(simp ~ dom_wild + leaf_stage + repro.stage.simp, data=diversity.ps)
m3.2 <- glmmTMB(simp ~ dom_wild + leaf_stage, data=diversity.ps); DHARMa::testResiduals(DHARMa::simulateResiduals(m3))
AICctab(m1, m2, m3, m3.2) # m3 best

# Manipulate random FX structure - nothing to manipulate
m4 <- glmmTMB(simp ~ dom_wild*leaf_stage + repro.stage.simp, data=diversity.ps); DHARMa::testResiduals(DHARMa::simulateResiduals(m4))
AICctab(m3, m4) # m3, additive model still best

# Test ancestry
m5 <- glmmTMB(simp ~ dom_wild + leaf_stage + repro.stage.simp + percent_FAL + percent_CAER, data=diversity.ps)
m6 <- glmmTMB(simp ~ dom_wild + leaf_stage + repro.stage.simp + percent_CAER, data=diversity.ps)
m7 <- glmmTMB(simp ~ dom_wild + leaf_stage + repro.stage.simp + percent_FAL, data=diversity.ps)
AICctab(m3, m5, m6, m7)

topmodel.simp.ps <- glmmTMB(simp ~ dom_wild + leaf_stage + repro.stage.simp + percent_FAL + percent_CAER, data=diversity.ps)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.simp.ps)) # hmm not fabulous

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.simp.ps, terms=c("dom_wild", "leaf_stage")))
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.simp.ps)
ests$model_df.residual <- df.residual(topmodel.simp.ps); ests

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.simp.ps, type= "II")) # no interactions, so use type II
sig1$response <- rownames(sig1)
sig0 <- subset(sig1, response == "dom_wild")
sig0$test.statistic <- "chisq"
names(sig0)[names(sig0) == 'Df'] <- 'statistic.df'
names(sig0)[names(sig0) == 'Chisq'] <- 'statistic'
names(sig0)[names(sig0) == 'Pr..Chisq.'] <- 'p.value'
rownames(sig0) <- NULL
sig <- subset(sig0, select=c("test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- cbind(ests, sig); pred_sig

##### E) Bootstrap CI 
# Here, I'll extract the effect of domestication within each leaf stage, due to the significant interaction 
# (Used this great post as a resource: https://www.r-bloggers.com/bootstrapping-follow-up-contrasts-for-within-subject-anovas-part-2/)
# Also used this helpful page to re-level the order of comparisons in emmeans (for some reason it always considers domestic the 'control' group, so estimates are negative rather than positive). 
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
# For example - 
normal <- pairs(emmeans(topmodel.simp.ps, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.simp.ps, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.simp.ps, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
summ <- data.frame(beta = b1$t0)
summ$leaf_stage <- rownames(summ)
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
conf$leaf_stage <- rownames(conf)
ci <- merge(summ, conf); ci

##### F) Calculate % changes
change.simp.ps <- merge(pred_sig, ci); change.simp.ps
change.simp.ps$percent.change <- change.simp.ps$beta/change.simp.ps$wild.pred*100
change.simp.ps$percent.change.lower <- change.simp.ps$beta.CI_lower.boot/change.simp.ps$wild.pred*100
change.simp.ps$percent.change.upper <- change.simp.ps$beta.CI_upper.boot/change.simp.ps$wild.pred*100

change.simp.ps$response <- "simpson"
change.simp.ps$scale <- "per-age-class"
change.simp.ps$model <- "status_only"

##### G) Put it all together & re-order for merging
change.simp.ps <- change.simp.ps[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.simp.ps

#--------------------------#
#  WITH peak area per leaf #
#--------------------------#
# Manipulate random FX structure
simp1 <- glmmTMB(simp ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|dom_wild) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
simp2 <- glmmTMB(simp ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|dom_wild) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
simp3 <- glmmTMB(simp ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|dom_wild), data=diversity.ps)
simp4 <- glmmTMB(simp ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|leaf_stage) + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
simp5 <- glmmTMB(simp ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|dom_wild), data=diversity.ps)
simp6 <- glmmTMB(simp ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
AICctab(simp2, simp6) # only these converge, and simp6 best; area-diversity relationship varies w/in plant individual!!!!!

# manipulate fixed FX structure
simp6.1 <- glmmTMB(simp ~ dom_wild + leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
simp6.2 <- glmmTMB(simp ~ dom_wild*leaf_stage + total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
simp6.3 <- glmmTMB(simp ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps) #
AICctab(simp6.1, simp6.2, simp6.3) # simp6.3 best; 3-way interaction
simp6.4 <- glmmTMB(simp ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + plant.volume_centered + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
simp6.5 <- glmmTMB(simp ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + plant.volume_centered + repro.stage.simp + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
AICctab(simp6.3, simp6.4) # simp6.3 still best

# Check if ancestry improves the model
simp6.6 <- glmmTMB(simp ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + percent_FAL + percent_CAER + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps) # doesn't converge
simp6.7 <- glmmTMB(simp ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + percent_FAL + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
simp6.8 <- glmmTMB(simp ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + percent_CAER + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
AICctab(simp6.3, simp6.7,simp6.8) # simp6.7 best

topmodel.simp.ps_total <-  glmmTMB(simp ~ dom_wild*leaf_stage*total_peak_area_perLeaf.log.cent + percent_FAL + (1+total_peak_area_perLeaf.log.cent|plant_pop), data=diversity.ps)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.simp.ps_total)) # also not great, but not a lot I can tweak
summary(topmodel.simp.ps_total) 

##### C) Pull out predicted means for wild & domestic + model nobs + df
mean(diversity.ps$total_peak_area_perLeaf.log.cent)
ests0 <- data.frame(ggpredict(topmodel.simp.ps_total, terms=c("dom_wild", "leaf_stage", "total_peak_area_perLeaf.log.cent [-1.997946e-16]")))
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.simp.ps_total)
ests$model_df.residual <- df.residual(topmodel.simp.ps_total); ests

##### D) Calculate / append significance 
sig0 <- data.frame(pairs(emmeans(topmodel.simp.ps_total, ~ dom_wild|leaf_stage*total_peak_area_perLeaf.log.cent))) # post-hoc test 
names(sig0)[names(sig0) == 't.ratio'] <- 'statistic'
names(sig0)[names(sig0) == 'df'] <- 'statistic.df'
sig0$test.statistic <- "t.ratio"
sig <- subset(sig0, select=c("leaf_stage", "test.statistic", "statistic", "statistic.df", "p.value")); sig

# Put together the model-predicted means + significance
pred_sig <- merge(ests, sig, by="leaf_stage"); pred_sig

##### E) Bootstrap CI 
# Here, I'll extract the effect of domestication within each leaf stage, due to the significant interaction 
# (Used this great post as a resource: https://www.r-bloggers.com/bootstrapping-follow-up-contrasts-for-within-subject-anovas-part-2/)
# Also used this helpful page to re-level the order of comparisons in emmeans (for some reason it always considers domestic the 'control' group, so estimates are negative rather than positive). 
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
# For example - 
normal <- pairs(emmeans(topmodel.simp.ps_total, ~ dom_wild|leaf_stage*total_peak_area_perLeaf.log.cent)); normal # estimates reversed
#versus
em <- emmeans(topmodel.simp.ps_total, specs = trt.vs.ctrl ~ dom_wild|leaf_stage*total_peak_area_perLeaf.log.cent); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage*total_peak_area_perLeaf.log.cent)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.simp.ps_total, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
summ <- data.frame(beta = b1$t0) # Works better on Will's computer
summ$leaf_stage <- rownames(summ)
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
conf$leaf_stage <- rownames(conf)
ci <- merge(summ, conf, by="leaf_stage"); ci

##### F) Calculate % changes
change.simp.ps_total <- merge(pred_sig, ci, by="leaf_stage"); change.simp.ps_total
change.simp.ps_total$percent.change <- change.simp.ps_total$beta/change.simp.ps_total$wild.pred*100
change.simp.ps_total$percent.change.lower <- change.simp.ps_total$beta.CI_lower.boot/change.simp.ps_total$wild.pred*100
change.simp.ps_total$percent.change.upper <- change.simp.ps_total$beta.CI_upper.boot/change.simp.ps_total$wild.pred*100

change.simp.ps_total$response <- "simpson"
change.simp.ps_total$scale <- "per-age-class"
change.simp.ps_total$model <- "status+totalpL"

##### G) Put it all together & re-order for merging
change.simp.ps_total <- change.simp.ps_total[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.simp.ps_total

#---------------------------------------------------------#
# Rbind the two scales + w & w/o total peak area together #
#---------------------------------------------------------#
change.simp.all <- rbind(change.simp.pp, change.simp.pp_mean, change.simp.ps, change.simp.ps_total)
change.simp.all 

# Name trait
change.simp.all$trait <- "simpson"

######################################
#------------------------------------#
# Rbind the diversity metrics + plot #
#------------------------------------#
######################################
head(change.abu.all); head(change.rich.all); head(change.H.all); head(change.simp.all)
str(change.abu.all) # 12 x 19
str(change.rich.all) # 8 x 19
str(change.H.all) # 8 x 19
str(change.simp.all) # 8 x 19

# Put them all together
diversity_saponins0 <- rbind(change.abu.all, change.rich.all, change.H.all, change.simp.all); diversity_saponins0

# Add hypotheses column - for total peak area, match for other univariate traits. For diversity, new hypothesis
#-- Total peak area
tpa <- subset(diversity_saponins0, trait == "Total saponins per leaf: avg"); tpa
tpa$hypothesis <- "Var ~ Domestication + Mean"
tpa$hypothesis[which(tpa$model == "status_only" & tpa$response == "mean")] <- "Mean ~ Domestication"
tpa$hypothesis[which(tpa$model == "status_only" & tpa$response == "sd")] <- "Var ~ Domestication"
tpa$hypothesis <- factor(tpa$hypothesis, levels=c("Var ~ Domestication", "Var ~ Domestication + Mean", "Mean ~ Domestication"))

#-- Diversity
div <- subset(diversity_saponins0, trait != "Total saponins per leaf: avg"); div
div$hypothesis <- "Diversity ~ Domestication"
div$hypothesis[which(div$model == "status+meanpL" | div$model == "status+totalpL")] <- "Diversity ~ Domestication + Compound production"
div
# - re-combine
diversity_saponins <- rbind(tpa, div); diversity_saponins

# Rename / Relevel some variables for plotting
diversity_saponins$leaf_stage <- factor(diversity_saponins$leaf_stage, levels=c("young", "middle", "old", "all"))
# Add a character column for significance
# Add a character column for significance
# . = trend, 0.1 level; * = .05 level; ** = .01; *** = <.001
diversity_saponins$sig.asterisk <- ""
diversity_saponins$sig.asterisk[which(diversity_saponins$p.value <.15)] <- "."
diversity_saponins$sig.asterisk[which(diversity_saponins$p.value <.055)] <- "*"
diversity_saponins$sig.asterisk[which(diversity_saponins$p.value <.015)] <- "**"
diversity_saponins$sig.asterisk[which(diversity_saponins$p.value <.0015)] <- "***"
diversity_saponins
