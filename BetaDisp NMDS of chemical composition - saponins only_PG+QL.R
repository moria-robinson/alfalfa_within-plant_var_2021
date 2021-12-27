##############################
# NDMS analyses of chemistry # 
##############################
# Goal: quantify dissimilarity in chemical composition (using dispersion). 
# 1) Run an NMDS with each point = a leaf
#     - Calculate 1 centroid per plant (n = 9 leaves per plant)
#     - Calculate the average distance of each leaf to the centroid to that plant
#     - This will give 1 value/plant that represents multivariate dispersion in chemistry among leaves of that plant
#     - Plot wild versus domestic

# 2) Run an NMDS with each point = a leaf
#     - Calculate 1 centroid per plant X leaf age (n = 3 leaves)
#     - Calculate the average distance of each leaf to the centroid of that leaf age group 
#     - Compare wild vs domestic, including plant individual as a covariate

# Do this with different slices of data:
# - raw peak intensities
# - presence/absence of compounds
# - standardized/proportional peak intensities 
# To varying degrees, these methods avoid strong leverage by abundant compounds / give different weight to rare compounds

############ Packages ############
library(vegan)
library(Rmisc)
library(ggplot2)
library(glmmTMB)
library(MuMIn)
library(ggeffects)
library(tidyr)
library(lme4)
library(boot)
library(emmeans)
library(bbmle)
library(sjPlot)

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

# Rename
data_all <- data_saps2; nrow(data_saps2); nrow(data_all)

###############################################
############ Data manipulation ################ 
###############################################
head(data_all); nrow(data_all)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Calculate proportional abundances of each compound, in each leaf #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
head(data_all)
total_sample_id <- aggregate(peak_intensity_norm_noNoise ~ sample_id, data=data_all, sum); nrow(total_sample_id); head(total_sample_id); names(total_sample_id) <- c("sample_id", "total_peak_intensity_norm_noNoise"); head(total_sample_id); nrow(total_sample_id) # 540 rows
data_all2 <- merge(data_all, total_sample_id, by="sample_id"); nrow(data_all); nrow(data_all2); head(data_all2)
data_all2$peak_intensity_norm_noNoise.prop <- data_all2$peak_intensity_norm_noNoise/data_all2$total_peak_intensity_norm_noNoise

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Calculate presence/absence of each compound, in each leaf #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
data_all2$presence.absence <- 0
data_all2$presence.absence[which(data_all2$peak_intensity_norm_noNoise > 0)] <- 1

# Add in a per plant + leaf stage grouping variable, for using in the NMDS
data_all2$plant_pop.leaf_stage <- paste(data_all2$plant_pop, data_all2$leaf_stage, sep=".")

# Make an info DF with the relevant covariates (useful in some NMDS visualizations, I think?)
data_all2$temp <- 1
info <- aggregate(temp ~sample_id + plant_pop + dom_wild + leaf_stage + repro.stage.simp + percent_FAL + percent_HEM + percent_CAER, data=data_all2, length); info$temp <- NULL; 
data_all2$temp <- NULL
head(info); nrow(info)

# Work with subsets of this DF
head(data_all2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Calculate total compound amounts, at the appropriate scales - not for use in NMDS, but as covariates in GLMMs #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# = average TOTAL peak areas at each scale (across all 9 leaves or across all 3 leaves/stage) to use in models later as a covariate; does producing more compounds overall at a given scale >> associate w/ greater dissimilarity?

# 1) first, sum up all peaks per leaf to get total peak area/leaf, across all saponins
area_pL <- do.call(data.frame, (aggregate(peak_intensity_norm_noNoise ~ sample_id + plant_pop + dom_wild + leaf_stage + repro.stage.simp + percent_FAL + percent_HEM + percent_CAER, data=data_all, function(x) c(sum = sum(x), N = length(x))))); head(area_pL); str(area_pL); nrow(area_pL)
head(area_pL) # N columnn - good, summed across 86 compounds per leaf
area_pL$peak_intensity_norm_noNoise.N <- NULL; head(area_pL)
names(area_pL) <- c("sample_id", "plant_pop", "dom_wild", "leaf_stage", "repro.stage.simp", "percent_FAL", "percent_HEM", "percent_CAER", "total_peak_area_perLeaf") 
head(area_pL); nrow(area_pL) # 540 measures; one per leaf

# 2) Mean peak area per plant (averaged across all N=9 leaves/plant : covariate for beta diversity among all leaves/plant)
mean.peak.area_raw.pp <- do.call(data.frame, (aggregate(total_peak_area_perLeaf ~ plant_pop, data=area_pL, function(x) c(mean = mean(x), N = length(x))))); head(mean.peak.area_raw.pp); str(mean.peak.area_raw.pp); nrow(mean.peak.area_raw.pp) # 60 average peak intensities, one per plant, N=9 observations per row
mean.peak.area_raw.pp$total_peak_area_perLeaf.N <- NULL; names(mean.peak.area_raw.pp) <- c("plant_pop", "total_peak_area_perLeaf_mean"); head(mean.peak.area_raw.pp) 

# 3) Mean peak area per leaf (averaged across N=3 leaves/plant/age class : covariate for beta diversity among leaves within each age class)
mean.peak.area_raw.ps <- do.call(data.frame, (aggregate(total_peak_area_perLeaf ~ plant_pop + leaf_stage, data=area_pL, function(x) c(mean = mean(x), N = length(x))))); head(mean.peak.area_raw.ps); str(mean.peak.area_raw.ps); nrow(mean.peak.area_raw.ps) # 180 average peak intensities; one per set of 3 leaves/age class
mean.peak.area_raw.ps$total_peak_area_perLeaf.N <- NULL; names(mean.peak.area_raw.ps) <- c("plant_pop","leaf_stage", "total_peak_area_perLeaf_mean"); head(mean.peak.area_raw.ps) 

############################################################################
# Follow procedure from Betadisp to visualize dispersion around centroids ## 
# SKIP TO LINE 180 - RETURN TO MAKE A PLOT FOR ANY OF THE BETADISPS BELOW  #
############################################################################
# Calculating and plotting centroids of NMDS Result - from http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/betadisper.html
# Nice tutorial using the iris dataset https://mattsigal.github.io/eqcov_supp/betadisp-ex.html
## Plot the dispersions in their top 2 dimensions. Helps to diagnose how well the method is doing.
# TO DO THIS: RUN ONE OF THE CHUNKS BELOW FIRST TO GENERATE DISTANCE MATRICES FOR A PARTICULAR TRANSFORMATION OF THE DATA
# Can run for any of the mod.plants (distance matrices) calculated below (e.g. with raw, proportional, or presence-absence of compounds)
# % var explained from PCA 1 and 2:
mod.plant$eig[1]/sum(mod.plant$eig); mod.plant$eig[2]/sum(mod.plant$eig) 
# Raw peak intensities: 24.5 % & 19.2% 
# Proportional abundances: 32.9% & 15.7% 
# Presence/absence: 37.1% and 24.0% 
mod.plant
str(mod.plant) # one centroid per plant
permutest(mod.plant)
anova(mod.plant)
# Pull out coordinates for each sample, along the first two PCA axes (PCA1: 32.5% of var, PCA2: 15.52% of var)
temp <- data.frame(mod.plant$vectors); head(temp); str(temp)
sample_id <- rownames(temp)
PCoA1 <- temp$PCoA1
PCoA2 <- temp$PCoA2
mod.plant_PCA1 <- data.frame(sample_id, PCoA1, PCoA2); head(mod.plant_PCA1); nrow(mod.plant_PCA1)
info <- subset(compounds_wide,  select= c("dom_wild", "plant_pop", "leaf_stage", "sample_id"))
mod.plant_PCA <- merge(mod.plant_PCA1, info, by="sample_id"); head(mod.plant_PCA)
mod.plant_PCA$type <- "leaf"

# Pull out group centroids:
str(mod.plant$centroids)
temp2 <- data.frame(mod.plant$centroids)
cent0 <- subset(temp2, select=c("PCoA1", "PCoA2")); cent0$plant_pop <- rownames(cent0); cent0
# Merge in dom_wild status with each plant pop centroid
info$tally <- 1; info2 <- aggregate(tally ~ plant_pop + dom_wild, data=info, length); info2$tally <- NULL; head(info2); nrow(info2) # 60
cent <- merge(cent0, info2, by="plant_pop"); head(cent)
cent$type <- "centroid"

# Calculate segments
seg0 = mod.plant_PCA[, c("PCoA1", "PCoA2", "plant_pop")]; head(seg0)
tmp = cent0; names(tmp) <- c("PCoA1_ctr", "PCoA2_ctr", "plant_pop"); head(tmp)
seg1 = join(seg0, tmp, "plant_pop"); head(seg1)
seg <- merge(seg1, info, by="plant_pop"); head(seg)

plant_PCA <- rbind.fill(mod.plant_PCA, cent)

ggplot() +
  geom_segment(data = seg, aes(x = PCoA1, y = PCoA2, xend = PCoA1_ctr, yend = PCoA2_ctr, colour = dom_wild)) +
  geom_point(data = plant_PCA, aes(x = PCoA1, y = PCoA2, colour = dom_wild, fill = dom_wild, shape = type, size=type)) +
  scale_shape_manual(values=c("leaf"=16, "centroid"=21)) +
  scale_size_manual(values=c("leaf"=.5, "centroid"=3)) +
  scale_color_manual("Domestication status", values=c("domestic" = "aquamarine3", "wild" = "chocolate1")) +
  scale_fill_manual(values=c("domestic" = "aquamarine3", "wild" = "chocolate1")) +
  coord_equal() +  
  theme_bw() +
  theme(panel.border = element_rect(colour="black"), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), strip.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.y=element_text(size=12, margin=margin(t = 0, r = 4, b = 0, l = 0)), axis.text.x=element_text(size=12, margin=margin(t = 4, r = 0, b = 0, l = 0)), axis.title.y=element_text(size=13,face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x=element_text(size=13,face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  guides(size=FALSE, shape = FALSE, fill=FALSE)
# 8 x 8 is nice.

#-#-#-#-#-#-#
# Analyses #
#-#-#-#-#-#-#
#------------------------------------#
#------##--------##--------##-------##
# 1)     Raw peak intensities        #
#------##--------##--------##-------##
#------------------------------------#
head(data_all2)
compounds <- subset(data_all2, select=c("plant_pop", "dom_wild", "leaf_stage", "plant_pop.leaf_stage", "sample_id", "compound", "peak_intensity_norm_noNoise")); str(compounds); head(compounds)
# Reshape to wide
compounds_wide <- reshape(compounds, idvar = c("dom_wild", "plant_pop", "leaf_stage", "plant_pop.leaf_stage", "sample_id"), timevar = "compound", direction = "wide")
compounds_wide[c(1:10), c(1:10)] 
nrow(compounds_wide); length(unique(compounds_wide$sample_id)); ncol(compounds_wide) # 540 rows, one row per leaf, and 91 columns - 86 compounds plus some ID variables 
rownames(compounds_wide) <- compounds_wide$sample_id
head(compounds_wide)
## Bray-Curtis distances between samples
dst <- vegdist(compounds_wide[,6:ncol(compounds_wide)], method="bray")
## Calculate multivariate dispersions
plant <- betadisper(dst, compounds_wide$plant_pop) # group = plant ID (N=9 leaves)
stage <- betadisper(dst, compounds_wide$plant_pop.leaf_stage) # group = leaf age within plant (N = 3 leaves)

# Calculate mean distance to centroid among 9 leaves per plant
dist.plant <- data.frame(data.frame(plant$distances)); head(dist.plant); nrow(dist.plant) # 540; one per leaf
dist.plant$sample_id <- rownames(dist.plant); head(dist.plant); nrow(dist.plant)
centroid_dist.plant <- merge(dist.plant, info, by="sample_id"); head(centroid_dist.plant); nrow(centroid_dist.plant)
# calculate mean distance
dist.plant2 <- aggregate(plant.distances ~ plant_pop + dom_wild, data=centroid_dist.plant, mean); head(dist.plant2); nrow(dist.plant2) # 60; one mean distance per plant
# Merge back plant-level covariates and total peak area
head(covariates); nrow(covariates)
sap.dist.pp0 <- merge(dist.plant2, covariates, by="plant_pop"); head(sap.dist.pp0); nrow(sap.dist.pp0)
# merge in the total compound production covariate (average total leaf saponins/plant)
sap.dist.pp <- merge(sap.dist.pp0, mean.peak.area_raw.pp, by="plant_pop"); head(sap.dist.pp); nrow(sap.dist.pp) # 60

# Calculate mean distance to centroid among 3 leaves/stage/plant
dist.stage <- data.frame(data.frame(stage$distances)); str(dist.stage); nrow(dist.stage)  # 540; also one per leaf
dist.stage$sample_id <- rownames(dist.stage); head(dist.stage); nrow(dist.stage)
centroid_dist.stage <- merge(dist.stage, info, by="sample_id"); head(centroid_dist.stage); nrow(centroid_dist.stage)
#-calculate the mean distance from each leaf to the per-plant/per-stage centroid
dist.stage2 <- aggregate(stage.distances ~ plant_pop + dom_wild + leaf_stage, data=centroid_dist.stage, mean); head(dist.stage2); nrow(dist.stage2) # 180; 1 mean distance per plant + leaf stage
# Merge back plant-level covariates
sap.dist.ps0 <- merge(dist.stage2, covariates, by="plant_pop"); head(sap.dist.ps0); nrow(sap.dist.ps0)
sap.dist.ps <- merge(sap.dist.ps0, mean.peak.area_raw.ps, by=c("plant_pop", "leaf_stage")); head(sap.dist.ps); nrow(sap.dist.ps) # down to 180

sap.dist.pp$leaf_stage <- "all"
names(sap.dist.pp)[names(sap.dist.pp) == 'plant.distances'] <- 'cent.dist'
names(sap.dist.ps)[names(sap.dist.ps) == 'stage.distances'] <- 'cent.dist'

##############################
#  a)    Per - plant         #
##############################
#---------------------------------#
#   Distance from centroid alone  # 
#---------------------------------#
head(sap.dist.pp); hist(sap.dist.pp$cent.dist); hist(log(sap.dist.pp$cent.dist)) # neither great; use untransformed &  check resids
m1 <- glmmTMB(cent.dist ~ dom_wild + plant.volume_centered, data=sap.dist.pp)
m1.2 <- glmmTMB(cent.dist ~ dom_wild, data=sap.dist.pp)
AICctab(m1,m1.2)
m2 <- glmmTMB(cent.dist ~ dom_wild + percent_FAL + percent_CAER, data=sap.dist.pp)
m3 <- glmmTMB(cent.dist ~ dom_wild + percent_CAER, data=sap.dist.pp)
m4 <- glmmTMB(cent.dist ~ dom_wild + percent_FAL, data=sap.dist.pp)
AICctab(m1.2, m2, m3, m4)

topmodel.pp <- glmmTMB(cent.dist ~ dom_wild + percent_CAER, data=sap.dist.pp)
summary(topmodel.pp)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.pp)) # looks good. Outliers, but barely significantly so.

##### Pull out predicted means for wild & domestic + model nobs + df
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
pred_sig <- cbind(ests, sig); pred_sig

##### E) Bootstrap CI 
# This is from the March 15, 2020 glmmTMB vignette. See "isLMM.glmmTMB" section. Instead of "$zi", in their example, I am interested in the cond estimate for domestic - so I name that specifically within the function
fixef(topmodel.pp)
fixef(topmodel.pp)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object

b1 <- lme4::bootMer(topmodel.pp, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")

summ <- data.frame(beta = b1$t0); summ
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci <- merge(summ, conf); ci

##### F) Calculate % changes
change.pp <- cbind(pred_sig, ci); change.pp
change.pp$percent.change <- change.pp$beta/change.pp$wild.pred*100
change.pp$percent.change.lower <- change.pp$beta.CI_lower.boot/change.pp$wild.pred*100
change.pp$percent.change.upper <- change.pp$beta.CI_upper.boot/change.pp$wild.pred*100

change.pp$leaf_stage <- "all"
change.pp$response <- "distance to centroid"
change.pp$scale <- "per-plant"
change.pp$model <- "status_only"

head(change.pp)
##### G) Put it all together & re-order for merging
change.pp <- change.pp[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.pp

#------------------------------------------------------#
##  Distance from centroid + total peak area per plant # 
#------------------------------------------------------#
head(sap.dist.pp); nrow(sap.dist.pp); hist(sap.dist.pp$total_peak_area_perLeaf_mean) # actually not bad! 

topmodel.pp_total <- glmmTMB(cent.dist ~ dom_wild + total_peak_area_perLeaf_mean + percent_CAER, data=sap.dist.pp)
summary(topmodel.pp_total) # strong negative effect of total compound amount on beta div. Interesting!
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.pp_total)) # looks fine; a single outlier

##### Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.pp_total, terms=c("dom_wild"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted")); ests1
ests <- spread(ests1, x, predicted); ests
names(ests) <- c("wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.pp_total)
ests$model_df.residual <- df.residual(topmodel.pp_total); ests

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.pp_total, type= "II")) # no interactions, so use type II
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
fixef(topmodel.pp_total)
fixef(topmodel.pp_total)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object

b1 <- lme4::bootMer(topmodel.pp_total, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")

summ <- data.frame(beta = b1$t0); summ
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci <- merge(summ, conf); ci

##### F) Calculate % changes
change.pp_total <- cbind(pred_sig, ci); change.pp_total
change.pp_total$percent.change <- change.pp_total$beta/change.pp_total$wild.pred*100
change.pp_total$percent.change.lower <- change.pp_total$beta.CI_lower.boot/change.pp_total$wild.pred*100
change.pp_total$percent.change.upper <- change.pp_total$beta.CI_upper.boot/change.pp_total$wild.pred*100

change.pp_total$leaf_stage <- "all"
change.pp_total$response <- "distance to centroid"
change.pp_total$scale <- "per-plant"
change.pp_total$model <- "status + mean peak area / group"

##### G) Put it all together & re-order for merging
change.pp_total <- change.pp_total[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.pp_total

###################################
#  b)    Per - Leaf stage        #
###################################
sap.dist.ps$cent.dist_log <- log(sap.dist.ps$cent.dist)
sap.dist.ps$plant_pop_dummy <- as.factor(sample(unique(sap.dist.ps$plant_pop), replace=T, nrow(sap.dist.ps))); head(sap.dist.ps)

#---------------------------------#
##  Distance from centroid alone ## 
#---------------------------------#
head(sap.dist.ps); hist(sap.dist.ps$cent.dist); hist(log(sap.dist.ps$cent.dist)) # skewed; better with log-transformation

m1 <- glmmTMB(cent.dist ~ dom_wild*leaf_stage + repro.stage + (1|plant_pop), data=sap.dist.ps)
m2 <- glmmTMB(cent.dist ~ dom_wild*leaf_stage + plant.volume_centered + repro.stage + (1|plant_pop), data=sap.dist.ps)
m3 <- glmmTMB(cent.dist ~ dom_wild*leaf_stage + plant.volume_centered + (1|plant_pop), data=sap.dist.ps)
m3_dummy <- glmmTMB(cent.dist ~ dom_wild*leaf_stage + plant.volume_centered + (1|plant_pop_dummy), data=sap.dist.ps)
AICctab(m1, m2, m3, m3_dummy) # m3 best
m4 <- glmmTMB(cent.dist ~ dom_wild + leaf_stage + plant.volume_centered + (1|plant_pop), data=sap.dist.ps)
AICctab(m3, m4) # m3, with interaction, much better
m5 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + (1|plant_pop), data=sap.dist.ps)
DHARMa::testResiduals(DHARMa::simulateResiduals(m4)) 
DHARMa::testResiduals(DHARMa::simulateResiduals(m5)) # better w/ logged response

# Test for including ancestry
m6 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + percent_FAL + percent_CAER + (1|plant_pop), data=sap.dist.ps)
m7 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + percent_CAER + (1|plant_pop), data=sap.dist.ps)
m8 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + percent_FAL + (1|plant_pop), data=sap.dist.ps)
AICctab(m3, m6, m7, m8) # m3 much better - cool!

# With or without plant size?
m3.1 <- glmmTMB(cent.dist ~ dom_wild*leaf_stage + (1|plant_pop), data=sap.dist.ps)
AICctab(m3, m3.1)

topmodel.ps <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + (1|plant_pop), data=sap.dist.ps)
topmodel.ps2 <- glmmTMB(cent.dist_log ~ dom_wild*leaf_stage + (1|plant_pop), data=sap.dist.ps)
summary(topmodel.ps)

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.ps, terms=c("dom_wild", "leaf_stage"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.ps)
ests$model_df.residual <- df.residual(topmodel.ps); ests

##### D) Calculate / append significance
sig0 <- data.frame(pairs(emmeans(topmodel.ps, ~ dom_wild|leaf_stage))) # post-hoc test 
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
normal <- pairs(emmeans(topmodel.ps, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.ps, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.ps2, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
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
change.ps <- cbind(pred_sig, ci); change.ps
change.ps$exp.beta <- exp(change.ps$beta)
change.ps$exp.beta.lower <- exp(change.ps$beta.CI_lower.boot)
change.ps$exp.beta.upper <- exp(change.ps$beta.CI_upper.boot)
change.ps
change.ps$percent.change <- (change.ps$exp.beta-1)*100
change.ps$percent.change.lower <-(change.ps$exp.beta.lower-1)*100
change.ps$percent.change.upper <- (change.ps$exp.beta.upper-1)*100

change.ps$response <- "distance to centroid"
change.ps$scale <- "per-age-class"
change.ps$model <- "status_only"

##### G) Put it all together & re-order for merging
change.ps <- change.ps[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.ps

#------------------------------------------------------#
##  Distance from centroid + total peak area per plant # 
#------------------------------------------------------#
head(sap.dist.ps); nrow(sap.dist.ps); hist(sap.dist.ps$total_peak_area_perLeaf_mean) # actually not bad! 

m1 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + total_peak_area_perLeaf_mean + (1+total_peak_area_perLeaf_mean|leaf_stage) + (1+total_peak_area_perLeaf_mean|plant_pop), data=sap.dist.ps) # doesn't converge
m2 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + total_peak_area_perLeaf_mean + (1+total_peak_area_perLeaf_mean|leaf_stage), data=sap.dist.ps) # doesn't converge

# Check for role of ancestry
m3 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + total_peak_area_perLeaf_mean + (1|plant_pop), data=sap.dist.ps)
m4 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + total_peak_area_perLeaf_mean + percent_FAL + percent_CAER + (1|plant_pop), data=sap.dist.ps)
m5 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + total_peak_area_perLeaf_mean + percent_FAL + (1|plant_pop), data=sap.dist.ps)
m6 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + total_peak_area_perLeaf_mean + percent_CAER + (1|plant_pop), data=sap.dist.ps)
AICctab(m3, m4, m5, m6)

#remove plant size?
m5.1 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + percent_FAL + (1|plant_pop), data=sap.dist.ps)
AICctab(m5, m5.1)

topmodel.ps_total <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + percent_FAL + (1|plant_pop), data=sap.dist.ps); car::Anova(topmodel.ps_total) # sig interaction - cool!
topmodel.ps_total2 <- glmmTMB(cent.dist_log ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + percent_FAL + (1|plant_pop), data=sap.dist.ps)
summary(topmodel.ps_total) 
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.ps_total)) # looks good

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.ps_total, terms=c("dom_wild", "leaf_stage"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.ps_total)
ests$model_df.residual <- df.residual(topmodel.ps_total); ests

##### D) Calculate / append significance
sig0 <- data.frame(pairs(emmeans(topmodel.ps_total, ~ dom_wild|leaf_stage))) # post-hoc test 
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
normal <- pairs(emmeans(topmodel.ps_total, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.ps_total, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.ps_total2, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
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
change.ps_total <- cbind(pred_sig, ci); change.ps_total
change.ps_total$exp.beta <- exp(change.ps_total$beta)
change.ps_total$exp.beta.lower <- exp(change.ps_total$beta.CI_lower.boot)
change.ps_total$exp.beta.upper <- exp(change.ps_total$beta.CI_upper.boot)
change.ps_total
change.ps_total$percent.change <- (change.ps_total$exp.beta-1)*100
change.ps_total$percent.change.lower <-(change.ps_total$exp.beta.lower-1)*100
change.ps_total$percent.change.upper <- (change.ps_total$exp.beta.upper-1)*100

change.ps_total$response <- "distance to centroid"
change.ps_total$scale <- "per-age-class"
change.ps_total$model <- "status + mean peak area / group"

##### G) Put it all together & re-order for merging
change.ps_total <- change.ps_total[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.ps_total

#------------------------------------------------------------------#
# Rbind the two scales together : whole-plant + across all leaves  #
#------------------------------------------------------------------#
change.all_raw <- rbind(change.pp, change.pp_total, change.ps, change.ps_total) 
change.all_raw$trait <- "BetaDisp: raw peak areas"
head(change.all_raw)

#------------------------------------------------------#
#------##--------##--------##-------##-------##-------##
# 3)   Proportional / standardized peak intensities    #
#------##--------##--------##-------##-------##-------##
#------------------------------------------------------#
# = raw peak areas divided by total peak area per sample; so, proportional abundance of each compound, in each sample
compounds <- subset(data_all2, select=c("plant_pop", "dom_wild", "leaf_stage", "plant_pop.leaf_stage", "sample_id", "compound", "peak_intensity_norm_noNoise.prop")); str(compounds); head(compounds)
# Reshape to wide
compounds_wide <- reshape(compounds, idvar = c("dom_wild", "plant_pop", "leaf_stage", "plant_pop.leaf_stage", "sample_id"), timevar = "compound", direction = "wide")
compounds_wide[c(1:10), c(1:10)]
nrow(compounds_wide); length(unique(compounds_wide$sample_id)) # 540! One per leaf
rownames(compounds_wide) <- compounds_wide$sample_id
head(compounds_wide); ncol(compounds_wide) # should be 91
## Bray-Curtis distances between samples
dst <- vegdist(compounds_wide[,6:ncol(compounds_wide)], method="bray")
## Calculate multivariate dispersions
plant <- betadisper(dst, compounds_wide$plant_pop)
stage <- betadisper(dst, compounds_wide$plant_pop.leaf_stage)

# Calculate mean distance to centroid among 9 leaves per plant; centroid is calculated among those 9 leaves per-plant
dist.plant <- data.frame(data.frame(plant$distances)); head(dist.plant); nrow(dist.plant) # 540; one per leaf
dist.plant$sample_id <- rownames(dist.plant); head(dist.plant); nrow(dist.plant)
centroid_dist.plant <- merge(dist.plant, info, by="sample_id"); head(centroid_dist.plant); nrow(centroid_dist.plant)
#-calculate the mean distance from each leaf to the per-plant centroid
dist.plant2 <- aggregate(plant.distances ~ plant_pop + dom_wild, data=centroid_dist.plant, mean); head(dist.plant2); nrow(dist.plant2) # 60; one mean distance per plant
# Merge back plant-level covariates and total peak area
head(covariates); nrow(covariates)
sap.dist.pp0 <- merge(dist.plant2, covariates, by="plant_pop"); head(sap.dist.pp0); nrow(sap.dist.pp0)
sap.dist.pp <- merge(sap.dist.pp0, mean.peak.area_raw.pp, by="plant_pop"); head(sap.dist.pp); nrow(sap.dist.pp) # down to 60

# Calculate mean distance to centroid among 3 leaves/stage/plant; centroid is calculated among those 3 leaves per-plant & stage
dist.stage <- data.frame(data.frame(stage$distances)); str(dist.stage); nrow(dist.stage)  # 540; also one per leaf
dist.stage$sample_id <- rownames(dist.stage); head(dist.stage); nrow(dist.stage)
centroid_dist.stage <- merge(dist.stage, info, by="sample_id"); head(centroid_dist.stage); nrow(centroid_dist.stage)
#-calculate the mean distance from each leaf to the per-plant/per-stage centroid
dist.stage2 <- aggregate(stage.distances ~ plant_pop + dom_wild + leaf_stage, data=centroid_dist.stage, mean); head(dist.stage2); nrow(dist.stage2) # 180; 1 mean distance per plant + leaf stage
# Merge back plant-level covariates
sap.dist.ps0 <- merge(dist.stage2, covariates, by="plant_pop"); head(sap.dist.ps0); nrow(sap.dist.ps0)
sap.dist.ps <- merge(sap.dist.ps0, mean.peak.area_raw.ps, by=c("plant_pop", "leaf_stage")); head(sap.dist.ps); nrow(sap.dist.ps) # down to 180

sap.dist.pp$leaf_stage <- "all"
names(sap.dist.pp)[names(sap.dist.pp) == 'plant.distances'] <- 'cent.dist'
names(sap.dist.ps)[names(sap.dist.ps) == 'stage.distances'] <- 'cent.dist'

##############################
#  a)    Per - plant         #
##############################
sap.dist.pp$cent.dist.log <- log(sap.dist.pp$cent.dist)

#---------------------------------#
#   Distance from centroid alone  # 
#---------------------------------#
head(sap.dist.pp); hist(sap.dist.pp$cent.dist); hist(sap.dist.pp$cent.dist.log) # log spreads the distributoun out a bit

m1 <- glmmTMB(cent.dist ~ dom_wild + plant.volume_centered + repro.stage, data=sap.dist.pp)
m1.2 <- glmmTMB(cent.dist ~ dom_wild + repro.stage, data=sap.dist.pp)
m1.3 <- glmmTMB(cent.dist ~ dom_wild + plant.volume_centered, data=sap.dist.pp); DHARMa::testResiduals(DHARMa::simulateResiduals(m1.3)) 
AICctab(m1, m1.2, m1.3)
m1.3_log <- glmmTMB(cent.dist.log ~ dom_wild + plant.volume_centered, data=sap.dist.pp); DHARMa::testResiduals(DHARMa::simulateResiduals(m1.3_log)) # QQ a bit tighter with the log transformation

# add ancestry
m1.4 <- glmmTMB(log(cent.dist) ~ dom_wild + plant.volume_centered + percent_FAL + percent_CAER, data=sap.dist.pp)
m1.5 <- glmmTMB(log(cent.dist) ~ dom_wild + plant.volume_centered + percent_FAL, data=sap.dist.pp)
m1.6 <- glmmTMB(log(cent.dist) ~ dom_wild + plant.volume_centered + percent_CAER, data=sap.dist.pp)
m1.7 <- glmmTMB(log(cent.dist) ~ dom_wild + plant.volume_centered, data=sap.dist.pp)
AICctab(m1.4, m1.5, m1.6, m1.7)

# With or without plant size?
m1.8 <- glmmTMB(log(cent.dist) ~ dom_wild + percent_CAER, data=sap.dist.pp)
AICctab(m1.6, m1.8)

topmodel.pp <- glmmTMB(log(cent.dist) ~ dom_wild + percent_CAER, data=sap.dist.pp)
topmodel.pp2 <- glmmTMB(cent.dist.log ~ dom_wild + percent_CAER, data=sap.dist.pp)
summary(topmodel.pp)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.pp)) # looks good

##### Pull out predicted means for wild & domestic + model nobs + df
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
pred_sig <- cbind(ests, sig); pred_sig

##### E) Bootstrap CI 
# This is from the March 15, 2020 glmmTMB vignette. See "isLMM.glmmTMB" section. Instead of "$zi", in their example, I am interested in the cond estimate for domestic - so I name that specifically within the function
fixef(topmodel.pp)
fixef(topmodel.pp)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object

b1 <- lme4::bootMer(topmodel.pp2, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")

summ <- data.frame(beta = b1$t0); summ
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci <- merge(summ, conf); ci

##### F) Calculate % changes - tweaked for log transformation
change.pp <- cbind(pred_sig, ci); change.pp
change.pp$exp.beta <- exp(change.pp$beta)
change.pp$exp.beta.lower <- exp(change.pp$beta.CI_lower.boot)
change.pp$exp.beta.upper <- exp(change.pp$beta.CI_upper.boot)
change.pp
change.pp$percent.change <- (change.pp$exp.beta-1)*100
change.pp$percent.change.lower <-(change.pp$exp.beta.lower-1)*100
change.pp$percent.change.upper <- (change.pp$exp.beta.upper-1)*100

change.pp$leaf_stage <- "all"
change.pp$response <- "distance to centroid"
change.pp$scale <- "per-plant"
change.pp$model <- "status_only"

head(change.pp)
##### G) Put it all together & re-order for merging
change.pp <- change.pp[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.pp

#------------------------------------------------------#
##  Distance from centroid + total peak area per plant # 
#------------------------------------------------------#
head(sap.dist.pp); nrow(sap.dist.pp); hist(sap.dist.pp$total_peak_area_perLeaf_mean) # actually not bad! 

# Check if plant size improves model
m1 <- glmmTMB(log(cent.dist) ~ dom_wild + total_peak_area_perLeaf_mean + plant.volume_centered, data=sap.dist.pp)
m2 <- glmmTMB(log(cent.dist) ~ dom_wild + total_peak_area_perLeaf_mean, data=sap.dist.pp)
AICctab(m1, m2) # no - remove

# Check if ancestry improves model
m3 <- glmmTMB(log(cent.dist) ~ dom_wild + total_peak_area_perLeaf_mean + percent_FAL + percent_CAER, data=sap.dist.pp)
m4 <- glmmTMB(log(cent.dist) ~ dom_wild + total_peak_area_perLeaf_mean + percent_FAL, data=sap.dist.pp)
m5 <- glmmTMB(log(cent.dist) ~ dom_wild + total_peak_area_perLeaf_mean + percent_FAL, data=sap.dist.pp)
AICctab(m2, m3, m4, m5) # m2 best - no ancestry

topmodel.pp_total <- glmmTMB(log(cent.dist) ~ dom_wild + total_peak_area_perLeaf_mean, data=sap.dist.pp)
topmodel.pp_total2 <- glmmTMB(cent.dist.log ~ dom_wild + total_peak_area_perLeaf_mean, data=sap.dist.pp)
summary(topmodel.pp_total) # interesting - here, total peak area doesn't matter. Cool!
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.pp_total)) # looks good

##### Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.pp_total, terms=c("dom_wild"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted")); ests1
ests <- spread(ests1, x, predicted); ests
names(ests) <- c("wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.pp_total)
ests$model_df.residual <- df.residual(topmodel.pp_total); ests

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.pp_total, type= "II")) # no interactions, so use type II
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
fixef(topmodel.pp_total)
fixef(topmodel.pp_total)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object

b1 <- lme4::bootMer(topmodel.pp_total2, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")

summ <- data.frame(beta = b1$t0); summ
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci <- merge(summ, conf); ci

##### F) Calculate % changes
change.pp_total <- cbind(pred_sig, ci); change.pp_total
change.pp_total$exp.beta <- exp(change.pp_total$beta)
change.pp_total$exp.beta.lower <- exp(change.pp_total$beta.CI_lower.boot)
change.pp_total$exp.beta.upper <- exp(change.pp_total$beta.CI_upper.boot)
change.pp_total
change.pp_total$percent.change <- (change.pp_total$exp.beta-1)*100
change.pp_total$percent.change.lower <-(change.pp_total$exp.beta.lower-1)*100
change.pp_total$percent.change.upper <- (change.pp_total$exp.beta.upper-1)*100

change.pp_total$leaf_stage <- "all"
change.pp_total$response <- "distance to centroid"
change.pp_total$scale <- "per-plant"
change.pp_total$model <- "status + mean peak area / group"

##### G) Put it all together & re-order for merging
change.pp_total <- change.pp_total[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.pp_total

###################################
#  b)    Per - Leaf stage         #
###################################
sap.dist.ps$cent.dist_log <- log(sap.dist.ps$cent.dist)

#---------------------------------#
##  Distance from centroid alone ## 
#---------------------------------#
head(sap.dist.ps); hist(sap.dist.ps$cent.dist); hist(log(sap.dist.ps$cent.dist)) # skewed; better with log-transformation

m1 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + repro.stage + (1|plant_pop), data=sap.dist.ps)
m2 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + (1|plant_pop), data=sap.dist.ps)
m3 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + repro.stage + (1|plant_pop), data=sap.dist.ps)
AICctab(m1, m2, m3) # m2 best
m4 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage +  plant.volume_centered + (1|plant_pop), data=sap.dist.ps)
m4.1 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + (1|plant_pop), data=sap.dist.ps)
AICctab(m2, m4, m4.1) # m4.1 best 

# Check role of ancestry...
m4.2 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + percent_FAL + percent_CAER + (1|plant_pop), data=sap.dist.ps)
m4.3 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + percent_FAL + (1|plant_pop), data=sap.dist.ps)
m4.4 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + percent_CAER + (1|plant_pop), data=sap.dist.ps)
AICctab(m4.1, m4.2, m4.3, m4.4) # m4.1 best

topmodel.ps <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage +  (1|plant_pop), data=sap.dist.ps)
topmodel.ps2 <- glmmTMB(cent.dist_log ~ dom_wild + leaf_stage +  (1|plant_pop), data=sap.dist.ps)
summary(topmodel.ps)

DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.ps)) # better w/ logged response

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.ps, terms=c("dom_wild", "leaf_stage"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.ps)
ests$model_df.residual <- df.residual(topmodel.ps); ests

##### D) Calculate / append significance
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

##### E) Bootstrap CI 
# Here, I'll extract the effect of domestication within each leaf stage
# Also used this helpful page to re-level the order of comparisons in emmeans (for some reason it always considers domestic the 'control' group, so estimates are negative rather than positive). 
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
# For example - 
normal <- pairs(emmeans(topmodel.ps, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.ps, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.ps2, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
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
change.ps <- cbind(pred_sig, ci); change.raw.ps
change.ps$exp.beta <- exp(change.ps$beta)
change.ps$exp.beta.lower <- exp(change.ps$beta.CI_lower.boot)
change.ps$exp.beta.upper <- exp(change.ps$beta.CI_upper.boot)
change.ps
change.ps$percent.change <- (change.ps$exp.beta-1)*100
change.ps$percent.change.lower <-(change.ps$exp.beta.lower-1)*100
change.ps$percent.change.upper <- (change.ps$exp.beta.upper-1)*100

change.ps$response <- "distance to centroid"
change.ps$scale <- "per-age-class"
change.ps$model <- "status_only"

##### G) Put it all together & re-order for merging
change.ps <- change.ps[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.ps

#------------------------------------------------------#
##  Distance from centroid + total peak area per plant #  
#------------------------------------------------------#
head(sap.dist.ps); nrow(sap.dist.ps); hist(sap.dist.ps$total_peak_area_perLeaf_mean) # actually not bad! 

m1 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + total_peak_area_perLeaf_mean + plant.volume_centered + (1|plant_pop), data=sap.dist.ps)
m1.1 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + total_peak_area_perLeaf_mean + (1|plant_pop), data=sap.dist.ps)
AICctab(m1, m1.1)
m2 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + total_peak_area_perLeaf_mean + (1 + total_peak_area_perLeaf_mean|leaf_stage) + (1 + total_peak_area_perLeaf_mean|plant_pop) + (1|plant_pop), data=sap.dist.ps)
m3 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + total_peak_area_perLeaf_mean + (1 + total_peak_area_perLeaf_mean|leaf_stage) + (1 + total_peak_area_perLeaf_mean|plant_pop), data=sap.dist.ps)
m4 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + total_peak_area_perLeaf_mean + (1 + total_peak_area_perLeaf_mean|leaf_stage), data=sap.dist.ps)
m4 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + total_peak_area_perLeaf_mean + (1 + total_peak_area_perLeaf_mean|plant_pop), data=sap.dist.ps)
# Only m1/m1.1 converges, simplest model

m5 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + (1|plant_pop), data=sap.dist.ps)
AICctab(m1.1, m5) # m1.1, additive better

# Ancestry....
m1.2 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + total_peak_area_perLeaf_mean + percent_FAL + percent_CAER + (1|plant_pop), data=sap.dist.ps)
m1.3 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + total_peak_area_perLeaf_mean + percent_FAL + (1|plant_pop), data=sap.dist.ps)
m1.4 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + total_peak_area_perLeaf_mean + percent_CAER + (1|plant_pop), data=sap.dist.ps)
AICctab(m1.1, m1.2, m1.3, m1.4) # m1.1 best, no ancestry

topmodel.ps_total <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + total_peak_area_perLeaf_mean + (1|plant_pop), data=sap.dist.ps); car::Anova(topmodel.ps_total)
topmodel.ps_total2 <- glmmTMB(cent.dist_log ~ dom_wild + leaf_stage + total_peak_area_perLeaf_mean + (1|plant_pop), data=sap.dist.ps)
summary(topmodel.ps_total)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.ps_total)) # looks good

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.ps_total, terms=c("dom_wild", "leaf_stage"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.ps_total)
ests$model_df.residual <- df.residual(topmodel.ps_total); ests

##### D) Calculate / append significance
sig1 <- data.frame(car::Anova(topmodel.ps_total, type= "II")) # no interactions, so use type II
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

##### E) Bootstrap CI 
# Here, I'll extract the effect of domestication within each leaf stage
# Also used this helpful page to re-level the order of comparisons in emmeans (for some reason it always considers domestic the 'control' group, so estimates are negative rather than positive). 
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
# For example - 
normal <- pairs(emmeans(topmodel.ps_total, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.ps_total, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.ps_total2, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
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
change.ps_total <- cbind(pred_sig, ci); change.ps_total
change.ps_total$exp.beta <- exp(change.ps_total$beta)
change.ps_total$exp.beta.lower <- exp(change.ps_total$beta.CI_lower.boot)
change.ps_total$exp.beta.upper <- exp(change.ps_total$beta.CI_upper.boot)
change.ps_total
change.ps_total$percent.change <- (change.ps_total$exp.beta-1)*100
change.ps_total$percent.change.lower <-(change.ps_total$exp.beta.lower-1)*100
change.ps_total$percent.change.upper <- (change.ps_total$exp.beta.upper-1)*100

change.ps_total$response <- "distance to centroid"
change.ps_total$scale <- "per-age-class"
change.ps_total$model <- "status + mean peak area / group"

##### G) Put it all together & re-order for merging
change.ps_total <- change.ps_total[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.ps_total

#------------------------------------------------------------------#
# Rbind the two scales together : whole-plant + across all leaves  #
#------------------------------------------------------------------#
change.all_prop <- rbind(change.pp, change.pp_total, change.ps, change.ps_total) 
change.all_prop$trait <- "BetaDisp: prop peak areas"
head(change.all_prop)

#------------------------------------#
#------##--------##--------##-------##
# 4)  Compound presence/absence      # 
#------##--------##--------##-------##
#------------------------------------#
compounds <- subset(data_all2, select=c("plant_pop", "dom_wild", "leaf_stage", "plant_pop.leaf_stage", "sample_id", "compound", "presence.absence")); str(compounds); head(compounds)
# Reshape to wide
compounds_wide <- reshape(compounds, idvar = c("dom_wild", "plant_pop", "leaf_stage", "plant_pop.leaf_stage", "sample_id"), timevar = "compound", direction = "wide")
compounds_wide[c(1:10), c(1:10)]
nrow(compounds_wide); length(unique(compounds_wide$sample_id)) # 540! One per leaf
rownames(compounds_wide) <- compounds_wide$sample_id
head(compounds_wide); ncol(compounds_wide) # should be 91
## Bray-Curtis distances between samples
dst <- vegdist(compounds_wide[,6:ncol(compounds_wide)], method="bray")
## Calculate multivariate dispersions
plant <- betadisper(dst, compounds_wide$plant_pop)
stage <- betadisper(dst, compounds_wide$plant_pop.leaf_stage)

# Calculate mean distance to centroid among 9 leaves per plant; centroid is calculated among those 9 leaves per-plant
dist.plant <- data.frame(data.frame(plant$distances)); head(dist.plant); nrow(dist.plant) # 540; one per leaf
dist.plant$sample_id <- rownames(dist.plant); head(dist.plant); nrow(dist.plant)
centroid_dist.plant <- merge(dist.plant, info, by="sample_id"); head(centroid_dist.plant); nrow(centroid_dist.plant)
#-calculate the mean distance from each leaf to the per-plant centroid
dist.plant2 <- aggregate(plant.distances ~ plant_pop + dom_wild, data=centroid_dist.plant, mean); head(dist.plant2); nrow(dist.plant2) # 60; one mean distance per plant
# Merge back plant-level covariates and total peak area
head(covariates); nrow(covariates)
sap.dist.pp0 <- merge(dist.plant2, covariates, by="plant_pop"); head(sap.dist.pp0); nrow(sap.dist.pp0)
sap.dist.pp <- merge(sap.dist.pp0, mean.peak.area_raw.pp, by="plant_pop"); head(sap.dist.pp); nrow(sap.dist.pp) # down to 60

# Calculate mean distance to centroid among 3 leaves/stage/plant; centroid is calculated among those 3 leaves per-plant & stage
dist.stage <- data.frame(data.frame(stage$distances)); str(dist.stage); nrow(dist.stage)  # 540; also one per leaf
dist.stage$sample_id <- rownames(dist.stage); head(dist.stage); nrow(dist.stage)
centroid_dist.stage <- merge(dist.stage, info, by="sample_id"); head(centroid_dist.stage); nrow(centroid_dist.stage)
#-calculate the mean distance from each leaf to the per-plant/per-stage centroid
dist.stage2 <- aggregate(stage.distances ~ plant_pop + dom_wild + leaf_stage, data=centroid_dist.stage, mean); head(dist.stage2); nrow(dist.stage2) # 180; 1 mean distance per plant + leaf stage
# Merge back plant-level covariates
sap.dist.ps0 <- merge(dist.stage2, covariates, by="plant_pop"); head(sap.dist.ps0); nrow(sap.dist.ps0)
sap.dist.ps <- merge(sap.dist.ps0, mean.peak.area_raw.ps, by=c("plant_pop", "leaf_stage")); head(sap.dist.ps); nrow(sap.dist.ps) # down to 180

##############################
#   a)   Per - plant         # 
##############################
head(sap.dist.pp); hist(sap.dist.pp$cent.dist) # wow very nice
#---------------------------------#
#   Distance from centroid alone  # 
#---------------------------------#
hist(sap.dist.pp$cent.dist)
ggplot(data=sap.dist.pp, aes(x=dom_wild, y=cent.dist)) +
         geom_boxplot() 

m1 <- glmmTMB(cent.dist ~ dom_wild + plant.volume_centered + repro.stage, data=sap.dist.pp)
m1.2 <- glmmTMB(cent.dist ~ dom_wild + plant.volume_centered, data=sap.dist.pp)
m1.3 <- glmmTMB(cent.dist ~ dom_wild + repro.stage, data=sap.dist.pp)
m1.4 <- glmmTMB(cent.dist ~ dom_wild, data=sap.dist.pp)
AICctab(m1, m1.2, m1.3, m1.4)

# Ancestry
m1.5 <- glmmTMB(cent.dist ~ dom_wild + percent_FAL + percent_CAER, data=sap.dist.pp)
m1.6 <- glmmTMB(cent.dist ~ dom_wild + percent_CAER, data=sap.dist.pp)
m1.7 <- glmmTMB(cent.dist ~ dom_wild + percent_FAL, data=sap.dist.pp)
AICctab(m1.4, m1.5, m1.6, m1.7)

topmodel.pp <- glmmTMB(cent.dist ~ dom_wild, data=sap.dist.pp)
summary(topmodel.pp) # nice!
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.pp)) # Very nice!

##### Pull out predicted means for wild & domestic + model nobs + df
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
pred_sig <- cbind(ests, sig); pred_sig

##### E) Bootstrap CI 
# This is from the March 15, 2020 glmmTMB vignette. See "isLMM.glmmTMB" section. Instead of "$zi", in their example, I am interested in the cond estimate for domestic - so I name that specifically within the function
fixef(topmodel.pp)
fixef(topmodel.pp)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object

b1 <- lme4::bootMer(topmodel.pp, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")

summ <- data.frame(beta = b1$t0); summ
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci <- merge(summ, conf); ci

##### F) Calculate % changes
change.pp <- cbind(pred_sig, ci); change.pp
change.pp$percent.change <- change.pp$beta/change.pp$wild.pred*100
change.pp$percent.change.lower <- change.pp$beta.CI_lower.boot/change.pp$wild.pred*100
change.pp$percent.change.upper <- change.pp$beta.CI_upper.boot/change.pp$wild.pred*100

change.pp$leaf_stage <- "all"
change.pp$response <- "distance to centroid"
change.pp$scale <- "per-plant"
change.pp$model <- "status_only"

head(change.pp)
##### G) Put it all together & re-order for merging
change.pp <- change.pp[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.pp

#------------------------------------------------------#
##  Distance from centroid + total peak area per plant # 
#------------------------------------------------------#
head(sap.dist.pp); nrow(sap.dist.pp); hist(sap.dist.pp$total_peak_area_perLeaf_mean) #

# Ancestry
m1 <- glmmTMB(cent.dist ~ dom_wild + total_peak_area_perLeaf_mean, data=sap.dist.pp)
m2 <- glmmTMB(cent.dist ~ dom_wild + total_peak_area_perLeaf_mean + percent_FAL + percent_CAER, data=sap.dist.pp)
m3 <- glmmTMB(cent.dist ~ dom_wild + total_peak_area_perLeaf_mean + percent_CAER, data=sap.dist.pp)
m4 <- glmmTMB(cent.dist ~ dom_wild + total_peak_area_perLeaf_mean + percent_FAL, data=sap.dist.pp)
AICctab(m1, m2, m3, m4)

topmodel.pp_total <- glmmTMB(cent.dist ~ dom_wild + total_peak_area_perLeaf_mean, data=sap.dist.pp)
summary(topmodel.pp_total)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.pp_total)) 

##### Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.pp_total, terms=c("dom_wild"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted")); ests1
ests <- spread(ests1, x, predicted); ests
names(ests) <- c("wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.pp_total)
ests$model_df.residual <- df.residual(topmodel.pp_total); ests

##### D) Calculate / append significance 
sig1 <- data.frame(car::Anova(topmodel.pp_total, type= "II")) # no interactions, so use type II
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
fixef(topmodel.pp_total)
fixef(topmodel.pp_total)$cond["dom_wilddomestic"] # confirm how to call the right component of the model object

b1 <- lme4::bootMer(topmodel.pp_total, FUN=function(x) fixef(x)$cond["dom_wilddomestic"], nsim=500, .progress="txt")

summ <- data.frame(beta = b1$t0); summ
conf <- data.frame(confint(b1, type = "perc"))
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
ci <- merge(summ, conf); ci

##### F) Calculate % changes
change.pp_total <- cbind(pred_sig, ci); change.pp_total
change.pp_total$percent.change <- change.pp_total$beta/change.pp_total$wild.pred*100
change.pp_total$percent.change.lower <- change.pp_total$beta.CI_lower.boot/change.pp_total$wild.pred*100
change.pp_total$percent.change.upper <- change.pp_total$beta.CI_upper.boot/change.pp_total$wild.pred*100

change.pp_total$leaf_stage <- "all"
change.pp_total$response <- "distance to centroid"
change.pp_total$scale <- "per-plant"
change.pp_total$model <- "status + mean peak area / group"

##### G) Put it all together & re-order for merging
change.pp_total <- change.pp_total[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.pp_total

###################################
#   b)   Per - Leaf stage         #
###################################
hist(sap.dist.ps$cent.dist)
hist(log(sap.dist.ps$cent.dist))
sap.dist.ps$cent.dist.log <- log(sap.dist.ps$cent.dist)

#---------------------------------#
##  Distance from centroid alone ## 
#---------------------------------#
head(sap.dist.ps); hist(sap.dist.ps$cent.dist); hist(log(sap.dist.ps$cent.dist)) # fine w/o transformation

m1 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + repro.stage + (1|plant_pop), data=sap.dist.ps)
m2 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + plant.volume_centered + (1|plant_pop), data=sap.dist.ps)
m3 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + repro.stage + (1|plant_pop), data=sap.dist.ps)
m4 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + (1|plant_pop), data=sap.dist.ps)
AICctab(m1, m2, m3, m4) # m4 best

m5 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + (1|plant_pop), data=sap.dist.ps)
AICctab(m4, m5) # m4, interactive model much better

# Ancestry?
m6 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + percent_FAL + percent_CAER + (1|plant_pop), data=sap.dist.ps)
m7 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + percent_FAL + (1|plant_pop), data=sap.dist.ps)
m8 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + percent_CAER + (1|plant_pop), data=sap.dist.ps)
AICctab(m4, m6, m7, m8) # m4 still best

topmodel.ps <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + (1|plant_pop), data=sap.dist.ps)
topmodel.ps2 <- glmmTMB(cent.dist.log ~ dom_wild*leaf_stage + (1|plant_pop), data=sap.dist.ps)
summary(topmodel.ps)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.ps)) # OK here! One outlier

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.ps, terms=c("dom_wild", "leaf_stage"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.ps)
ests$model_df.residual <- df.residual(topmodel.ps); ests

##### D) Calculate / append significance
sig0 <- data.frame(pairs(emmeans(topmodel.ps, ~ dom_wild|leaf_stage))) # post-hoc test 
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
normal <- pairs(emmeans(topmodel.ps, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.ps, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.ps2, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
summ <- data.frame(beta = b1$t0); summ
summ$leaf_stage <- rownames(summ)
conf <- data.frame(confint(b1, type = "perc")); conf
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
conf$leaf_stage <- rownames(conf)
ci <- merge(summ, conf, by="leaf_stage"); ci

##### F) Calculate % changes
change.ps <- cbind(pred_sig, ci); change.ps
change.ps$exp.beta <- exp(change.ps$beta)
change.ps$exp.beta.lower <- exp(change.ps$beta.CI_lower.boot)
change.ps$exp.beta.upper <- exp(change.ps$beta.CI_upper.boot)
change.ps
change.ps$percent.change <- (change.ps$exp.beta-1)*100
change.ps$percent.change.lower <-(change.ps$exp.beta.lower-1)*100
change.ps$percent.change.upper <- (change.ps$exp.beta.upper-1)*100

change.ps$response <- "distance to centroid"
change.ps$scale <- "per-age-class"
change.ps$model <- "status_only"

##### G) Put it all together & re-order for merging
change.ps <- change.ps[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.ps

#------------------------------------------------------#
##  Distance from centroid + total peak area per plant # 
#------------------------------------------------------#
head(sap.dist.ps); nrow(sap.dist.ps); hist(sap.dist.ps$total_peak_area_perLeaf_mean) # actually not bad! 

m1 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + plant.volume_centered + (1|plant_pop), data=sap.dist.ps)
m1.1 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + (1|plant_pop), data=sap.dist.ps)
m2 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + plant.volume_centered + (1 + total_peak_area_perLeaf_mean|leaf_stage) + (1 + total_peak_area_perLeaf_mean|plant_pop) + (1|plant_pop), data=sap.dist.ps)
m3 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + plant.volume_centered + (1 + total_peak_area_perLeaf_mean|leaf_stage) + (1 + total_peak_area_perLeaf_mean|plant_pop), data=sap.dist.ps)
m4 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + plant.volume_centered + (1 + total_peak_area_perLeaf_mean|leaf_stage), data=sap.dist.ps)
m5 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + plant.volume_centered + (1 + total_peak_area_perLeaf_mean|plant_pop), data=sap.dist.ps)
AICctab(m1, m1.1, m5) # m1/m1.1 best

m6 <- glmmTMB(log(cent.dist) ~ dom_wild + leaf_stage + total_peak_area_perLeaf_mean + (1|plant_pop), data=sap.dist.ps)
AICctab(m1.1, m6) # m1.1, interactive model better

# Ancestry
m1.2 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + percent_FAL + percent_CAER + (1|plant_pop), data=sap.dist.ps)
m1.3 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + percent_FAL + (1|plant_pop), data=sap.dist.ps)
m1.4 <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + percent_CAER + (1|plant_pop), data=sap.dist.ps)
AICctab(m1.1, m1.2, m1.3, m1.4) #m1.1 best

topmodel.ps_total <- glmmTMB(log(cent.dist) ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + (1|plant_pop), data=sap.dist.ps)
topmodel.ps_total2 <- glmmTMB(cent.dist.log ~ dom_wild*leaf_stage + total_peak_area_perLeaf_mean + (1|plant_pop), data=sap.dist.ps)
summary(topmodel.ps_total); Anova(topmodel.ps_total)
DHARMa::testResiduals(DHARMa::simulateResiduals(topmodel.ps_total)) # looks good

##### C) Pull out predicted means for wild & domestic + model nobs + df
ests0 <- data.frame(ggpredict(topmodel.ps_total, terms=c("dom_wild", "leaf_stage"))); ests0
ests1 <- subset(ests0, select=c("x", "predicted", "group"))
ests <- spread(ests1, x, predicted)
names(ests) <- c("leaf_stage", "wild.pred", "dom.pred"); ests

# Add in the model nobs/df
ests$model_nobs <- nobs(topmodel.ps_total)
ests$model_df.residual <- df.residual(topmodel.ps_total); ests

##### D) Calculate / append significance
sig0 <- data.frame(pairs(emmeans(topmodel.ps_total, ~ dom_wild|leaf_stage))) # post-hoc test 
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
normal <- pairs(emmeans(topmodel.ps_total, ~ dom_wild|leaf_stage)); normal # estimates reversed
#versus
em <- emmeans(topmodel.ps_total, specs = trt.vs.ctrl ~ dom_wild|leaf_stage); em$contrasts # estimates are in the correct order, positive when an increase with domestication - yay!

contrasts <- function(x){
  em <- emmeans(x, specs = trt.vs.ctrl ~ dom_wild|leaf_stage)
  con <- data.frame(em$contrasts)
  ests <- con$estimate
  names(ests) <- con$leaf_stage
  ests
}

# Bootstrapping...
b1 <- lme4::bootMer(topmodel.ps_total2, FUN=function(x) contrasts(x), nsim=500, .progress="txt")
summ <- data.frame(beta = b1$t0); summ
summ$leaf_stage <- rownames(summ)
conf <- data.frame(confint(b1, type = "perc")); conf
names(conf)[names(conf) == 'X2.5..'] <- 'beta.CI_lower.boot'
names(conf)[names(conf) == 'X97.5..'] <- 'beta.CI_upper.boot'
conf$leaf_stage <- rownames(conf)
ci <- merge(summ, conf, by="leaf_stage"); ci

##### F) Calculate % changes 
change.ps_total <- cbind(pred_sig, ci); change.ps_total
change.ps_total$exp.beta <- exp(change.ps_total$beta)
change.ps_total$exp.beta.lower <- exp(change.ps_total$beta.CI_lower.boot)
change.ps_total$exp.beta.upper <- exp(change.ps_total$beta.CI_upper.boot)
change.ps_total
change.ps_total$percent.change <- (change.ps_total$exp.beta-1)*100
change.ps_total$percent.change.lower <-(change.ps_total$exp.beta.lower-1)*100
change.ps_total$percent.change.upper <- (change.ps_total$exp.beta.upper-1)*100

change.ps_total$response <- "distance to centroid"
change.ps_total$scale <- "per-age-class"
change.ps_total$model <- "status + mean peak area / group"

##### G) Put it all together & re-order for merging
change.ps_total <- change.ps_total[c("model", "scale", "leaf_stage", "response", "model_nobs", "model_df.residual", "wild.pred", "dom.pred", "beta", "beta.CI_lower.boot", "beta.CI_upper.boot", "percent.change", "percent.change.lower", "percent.change.upper", "test.statistic", "statistic", "statistic.df", "p.value")]
change.ps_total

#------------------------------------------------------------------#
# Rbind the two scales together : whole-plant + across all leaves  #
#------------------------------------------------------------------#
change.all_pa <- rbind(change.pp, change.pp_total, change.ps, change.ps_total) 
change.all_pa$trait <- "BetaDisp: presence.absence"
head(change.all_pa)

#----------------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------#
# Rbind the three methods together: raw peak intensities, standardized peak intensities, + presence/absence  #
#----------------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------#
head(change.all_raw)
head(change.all_prop)
head(change.all_pa)

disp_saponins <- rbind(change.all_raw, change.all_prop, change.all_pa); disp_saponins
# Similarly to the diversity profiles, these three methods give different weight to abundant vs. rare compounds

# Rename / Relevel some variables for plotting
disp_saponins$hypothesis <- "Diversity ~ Domestication"
disp_saponins$hypothesis[which(disp_saponins$model == "status + mean peak area / group")] <- "Diversity ~ Domestication + Mean peak area / group"
disp_saponins$leaf_stage <- factor(disp_saponins$leaf_stage, levels=c("young", "middle", "old", "all"))
# Add a character column for significance
# Add a character column for significance
# . = trend, 0.1 level; * = .05 level; ** = .01; *** = <.001
disp_saponins$sig.asterisk <- ""
disp_saponins$sig.asterisk[which(disp_saponins$p.value <.15)] <- "."
disp_saponins$sig.asterisk[which(disp_saponins$p.value <.055)] <- "*"
disp_saponins$sig.asterisk[which(disp_saponins$p.value <.015)] <- "**"
disp_saponins$sig.asterisk[which(disp_saponins$p.value <.0015)] <- "***"
