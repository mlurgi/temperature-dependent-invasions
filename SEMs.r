## ---------------------------
##
## Script name: SEMs.r
##
## Purpose of script: This script implements structural equation models looking into the complex relationships
## between temperature and different community / network properties on simulated complex food webs to determine
## their influence on invasion success.
## Input data for the models includes a series of measurements perfomed over in-silico food webs and communities
## which dynamics are governed by these species interactions. Invasions experiments were performed on these 
## communities along a temperature gradient, and the outcome of the invasion (either successful or unsuccessful) is recorded.
##
## Author: Dr Miguel Lurgi
## Lecturer in Biosciences (Computational Ecology)
## Computational Ecology Lab - Department of Biosciences
## Swansea University, UK
## 
## and
##
## Centre for Biodiversity Theory and Modelling
## Theoretical and Experimental Ecology Station, CNRS, France
##
## Date Created: November-2018
##
## Copyright (c) Miguel Lurgi, 2018-2020
## Email: miguel.lurgi@swansea.ac.uk
##
## ---------------------------
##
## Notes:
##
## This script is provided as supplementary material for the paper:
## Sentis, Montoya & Lurgi (2020) Warming indirectly incrases invasion success in food webs. Uploaded to BioRXiv. https://doi.org/10.1101/2020.07.20.211516
##
## ---------------------------

require(lme4)
require(piecewiseSEM)

load('data-for-sems.rda')
#### After receiving comments from the reviewers from Proceedings B, we decided to extend the SEMs analysis and
#### incorporate 6 separate model specifications that allowed us to test specific hypotheses in the manuscript.
#### These models were then compared using AIC.

#### Model 1: Only indirect effects from temperature on invasion success via network / community properties,
#### with no indirect effects going through S (Temperature influences all variables and then all variables 
#### influence invasion success)

model_1_spec <- psem(
  glmer(invaded ~ S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_orig, nAGQ = 0),
  glmer(S ~ temperature_c + (1|id), data=cur_data_orig, family='poisson', nAGQ = 0),
  glmer(L ~ temperature_c + (1|id), data=cur_data_orig, nAGQ=0, family=MASS::negative.binomial(theta=1e5)),
  lmer(L.S ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(mfcl ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(modularity ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(gensd ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(vulsd ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(basal ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(intermediate ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(max_sim ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(log_com_cv ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(log_mean_cv ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(avg_ppmr ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(total_biomass ~ temperature_c + (1|id), data=cur_data_orig),
  invaded %~~% temperature_c,
  L %~~% S,
  L.S %~~% S,
  mfcl %~~% S,
  modularity %~~% S,
  gensd %~~% S,
  vulsd %~~% S,
  basal %~~% S,
  intermediate %~~% S,
  max_sim %~~% S,
  log_com_cv %~~% S,
  log_mean_cv %~~% S,
  avg_ppmr %~~% S,
  total_biomass %~~% S,
  L.S %~~% L,
  mfcl %~~% L,
  log_com_cv %~~% L,
  log_mean_cv %~~% L,
  basal %~~% L,
  intermediate %~~% L,
  modularity %~~% L,
  gensd %~~% L,
  vulsd %~~% L,
  avg_ppmr %~~% L,
  total_biomass %~~% L,
  max_sim %~~% L,
  basal %~~% L.S,
  intermediate %~~% L.S,
  modularity %~~% L.S,
  mfcl %~~% L.S,
  log_com_cv %~~% L.S,
  log_mean_cv %~~% L.S,
  gensd %~~% L.S,
  vulsd %~~% L.S,
  avg_ppmr %~~% L.S,
  total_biomass %~~% L.S,
  max_sim %~~% L.S,
  basal %~~% mfcl,
  intermediate %~~% mfcl,
  modularity %~~% mfcl,
  log_mean_cv %~~% mfcl,
  gensd %~~% mfcl,
  vulsd %~~% mfcl,
  avg_ppmr %~~% mfcl,
  max_sim %~~% mfcl,
  modularity %~~% avg_ppmr,
  modularity %~~% log_mean_cv,
  modularity %~~% basal,
  modularity %~~% intermediate,
  modularity %~~% max_sim,
  modularity %~~% gensd,
  modularity %~~% vulsd,
  modularity %~~% total_biomass,
  gensd %~~% vulsd,
  gensd %~~% avg_ppmr,
  gensd %~~% intermediate,
  gensd %~~% basal,
  gensd %~~% max_sim,
  gensd %~~% total_biomass,
  vulsd %~~% avg_ppmr,
  vulsd %~~% basal,
  vulsd %~~% intermediate,
  vulsd %~~% max_sim,
  vulsd %~~% log_com_cv,
  vulsd %~~% log_mean_cv,
  log_com_cv %~~% log_mean_cv,
  log_com_cv %~~% intermediate,
  log_com_cv %~~% basal,
  log_com_cv %~~% mfcl,
  log_com_cv %~~% total_biomass,
  log_mean_cv %~~% intermediate,
  log_mean_cv %~~% avg_ppmr,
  log_mean_cv %~~% total_biomass,
  log_mean_cv %~~% max_sim,
  avg_ppmr %~~% intermediate,
  avg_ppmr %~~% max_sim,
  avg_ppmr %~~% basal,
  intermediate %~~% total_biomass,
  max_sim %~~% total_biomass,
  max_sim %~~% basal,
  basal %~~% total_biomass,
  intermediate %~~% basal,
  intermediate %~~% max_sim,
  data=cur_data_orig
)

summary_model_1 <- summary(model_1_spec, conserve=TRUE)

### because standardised coefficients are for some reason not obtained automatedly from psem for responses 'invaded', 'S', and 'L'
### we calculate these a posteriori
### To calculate standardised coefficients for binomial responses we use the following approach developed in :
### https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#standardized-coefficients
### section 2.2.4

glm.model <- glmer(invaded ~ S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_orig, nAGQ = 0)

Beta.glm <- summary(glm.model)$coefficients[2:15, 1]
preds <- predict(glm.model, type = "link")

# Compute sd of error variance using theoretical variances
sd.y.LT <- sqrt(var(preds) + pi^2/3)
# Compute sd of x
sd.x <- apply(cur_data_orig[names(Beta.glm)], 2, sd)

summary_model_1$coefficients[1:14,8] <- round((Beta.glm * sd.x / sd.y.LT), 4)

### for S and L we follow the common standardisation of beta * (sdx/sdy)

sd.y <- apply(cur_data_orig['S'], 2, sd)
predictors <- summary_model_1$coefficients[which(summary_model_1$coefficients$Response == 'S'),]$Predictor
sd.x <- apply(cur_data_orig[as.character(predictors)], 2, sd)

summary_model_1$coefficients[which(summary_model_1$coefficients$Response == 'S'),]$Std.Estimate <- round(summary_model_1$coefficients[which(summary_model_1$coefficients$Response == 'S'),]$Estimate * (sd.x / sd.y), 4)

sd.y <- apply(cur_data_orig['L'], 2, sd)
predictors <- summary_model_1$coefficients[which(summary_model_1$coefficients$Response == 'L'),]$Predictor
sd.x <- apply(cur_data_orig[as.character(predictors)], 2, sd)

summary_model_1$coefficients[which(summary_model_1$coefficients$Response == 'L'),]$Std.Estimate <- round(summary_model_1$coefficients[which(summary_model_1$coefficients$Response == 'L'),]$Estimate * (sd.x / sd.y), 4)

save(summary_model_1, file='sem-model-1.RData')
write.csv(summary_model_1$coefficients, file='sem-model-1-coeffs.csv')

#### Model 2: Direct and indirect effects from temperature on invasion success
#### The same as model 1 but including the direct effect of temperature on invasion success

model_2_spec <- psem(
  glmer(invaded ~ temperature_c + S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_orig, nAGQ = 0),
  glmer(S ~ temperature_c + (1|id), data=cur_data_orig, family='poisson', nAGQ = 0),
  glmer(L ~ temperature_c + (1|id), data=cur_data_orig, nAGQ=0, family=MASS::negative.binomial(theta=1e5)),
  lmer(L.S ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(mfcl ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(modularity ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(gensd ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(vulsd ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(basal ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(intermediate ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(max_sim ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(log_com_cv ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(log_mean_cv ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(avg_ppmr ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(total_biomass ~ temperature_c + (1|id), data=cur_data_orig),
  L %~~% S,
  L.S %~~% S,
  mfcl %~~% S,
  modularity %~~% S,
  gensd %~~% S,
  vulsd %~~% S,
  basal %~~% S,
  intermediate %~~% S,
  max_sim %~~% S,
  log_com_cv %~~% S,
  log_mean_cv %~~% S,
  avg_ppmr %~~% S,
  total_biomass %~~% S,
  L.S %~~% L,
  mfcl %~~% L,
  log_com_cv %~~% L,
  log_mean_cv %~~% L,
  basal %~~% L,
  intermediate %~~% L,
  modularity %~~% L,
  gensd %~~% L,
  vulsd %~~% L,
  avg_ppmr %~~% L,
  total_biomass %~~% L,
  max_sim %~~% L,
  basal %~~% L.S,
  intermediate %~~% L.S,
  modularity %~~% L.S,
  mfcl %~~% L.S,
  log_com_cv %~~% L.S,
  log_mean_cv %~~% L.S,
  gensd %~~% L.S,
  vulsd %~~% L.S,
  avg_ppmr %~~% L.S,
  total_biomass %~~% L.S,
  max_sim %~~% L.S,
  basal %~~% mfcl,
  intermediate %~~% mfcl,
  modularity %~~% mfcl,
  log_mean_cv %~~% mfcl,
  gensd %~~% mfcl,
  vulsd %~~% mfcl,
  avg_ppmr %~~% mfcl,
  max_sim %~~% mfcl,
  modularity %~~% avg_ppmr,
  modularity %~~% log_mean_cv,
  modularity %~~% basal,
  modularity %~~% intermediate,
  modularity %~~% max_sim,
  modularity %~~% gensd,
  modularity %~~% vulsd,
  modularity %~~% total_biomass,
  gensd %~~% vulsd,
  gensd %~~% avg_ppmr,
  gensd %~~% intermediate,
  gensd %~~% basal,
  gensd %~~% max_sim,
  gensd %~~% total_biomass,
  vulsd %~~% avg_ppmr,
  vulsd %~~% basal,
  vulsd %~~% intermediate,
  vulsd %~~% max_sim,
  vulsd %~~% log_com_cv,
  vulsd %~~% log_mean_cv,
  log_com_cv %~~% log_mean_cv,
  log_com_cv %~~% intermediate,
  log_com_cv %~~% basal,
  log_com_cv %~~% mfcl,
  log_com_cv %~~% total_biomass,
  log_mean_cv %~~% intermediate,
  log_mean_cv %~~% avg_ppmr,
  log_mean_cv %~~% total_biomass,
  log_mean_cv %~~% max_sim,
  avg_ppmr %~~% intermediate,
  avg_ppmr %~~% max_sim,
  avg_ppmr %~~% basal,
  intermediate %~~% total_biomass,
  max_sim %~~% total_biomass,
  max_sim %~~% basal,
  basal %~~% total_biomass,
  intermediate %~~% basal,
  intermediate %~~% max_sim,
  data=cur_data_orig
)

summary_model_2 <- summary(model_2_spec, conserve=TRUE)

### because standardised coefficients are for some reason not obtained automatedly from psem for responses 'invaded', 'S', and 'L'
### we calculate these a posteriori
### To calculate standardised coefficients for binomial responses we use the following approach developed in :
### https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#standardized-coefficients
### section 2.2.4

glm.model <- glmer(invaded ~ temperature_c + S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_orig, nAGQ = 0)

Beta.glm <- summary(glm.model)$coefficients[2:16, 1]
preds <- predict(glm.model, type = "link")

# Compute sd of error variance using theoretical variances
sd.y.LT <- sqrt(var(preds) + pi^2/3)
# Compute sd of x
sd.x <- apply(cur_data_orig[names(Beta.glm)], 2, sd)

summary_model_2$coefficients[1:15,8] <- round((Beta.glm * sd.x / sd.y.LT), 4)

### for S and L we follow the common standardisation of beta * (sdx/sdy)

sd.y <- apply(cur_data_orig['S'], 2, sd)
predictors <- summary_model_2$coefficients[which(summary_model_2$coefficients$Response == 'S'),]$Predictor
sd.x <- apply(cur_data_orig[as.character(predictors)], 2, sd)

summary_model_2$coefficients[which(summary_model_2$coefficients$Response == 'S'),]$Std.Estimate <- round(summary_model_2$coefficients[which(summary_model_2$coefficients$Response == 'S'),]$Estimate * (sd.x / sd.y), 4)

sd.y <- apply(cur_data_orig['L'], 2, sd)
predictors <- summary_model_2$coefficients[which(summary_model_2$coefficients$Response == 'L'),]$Predictor
sd.x <- apply(cur_data_orig[as.character(predictors)], 2, sd)

summary_model_2$coefficients[which(summary_model_2$coefficients$Response == 'L'),]$Std.Estimate <- round(summary_model_2$coefficients[which(summary_model_2$coefficients$Response == 'L'),]$Estimate * (sd.x / sd.y), 4)

save(summary_model_2, file='sem-model-2.RData')
write.csv(summary_model_2$coefficients, file='sem-model-2-coeffs.csv')

#### Model 3: Deterministic interactions are most important. In this SEM, temperature effects
#### are exclusively modelled through the action of S. Temperature is not assumed to affect either
#### other community properties nor invasion success

model_3_spec <- psem(
  glmer(invaded ~ S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_orig, nAGQ = 0),
  glmer(S ~ temperature_c + (1|id), data=cur_data_orig, family='poisson', nAGQ = 0),
  glmer(L ~ S + (1|id), data=cur_data_orig, nAGQ=0, family=MASS::negative.binomial(theta=1e5)),
  lmer(L.S ~ S + (1|id), data=cur_data_orig),
  lmer(mfcl ~ S + (1|id), data=cur_data_orig),
  lmer(modularity ~ S + (1|id), data=cur_data_orig),
  lmer(gensd ~ S + (1|id), data=cur_data_orig),
  lmer(vulsd ~ S + (1|id), data=cur_data_orig),
  lmer(basal ~ S + (1|id), data=cur_data_orig),
  lmer(intermediate ~ S + (1|id), data=cur_data_orig),
  lmer(max_sim ~ S + (1|id), data=cur_data_orig),
  lmer(log_com_cv ~ S + (1|id), data=cur_data_orig),
  lmer(log_mean_cv ~ S + (1|id), data=cur_data_orig),
  lmer(avg_ppmr ~ S + (1|id), data=cur_data_orig),
  lmer(total_biomass ~ S + (1|id), data=cur_data_orig),
  invaded %~~% temperature_c,
  L %~~% temperature_c,
  L.S %~~% temperature_c,
  mfcl %~~% temperature_c,
  modularity %~~% temperature_c,
  gensd %~~% temperature_c,
  vulsd %~~% temperature_c,
  intermediate %~~% temperature_c,
  max_sim %~~% temperature_c,
  log_com_cv %~~% temperature_c,
  log_mean_cv %~~% temperature_c,
  total_biomass %~~% temperature_c,
  avg_ppmr %~~% temperature_c,
  L.S %~~% L,
  mfcl %~~% L,
  log_mean_cv %~~% L,
  basal %~~% L,
  intermediate %~~% L,
  modularity %~~% L,
  gensd %~~% L,
  vulsd %~~% L,
  avg_ppmr %~~% L,
  total_biomass %~~% L,
  max_sim %~~% L,
  basal %~~% L.S,
  intermediate %~~% L.S,
  modularity %~~% L.S,
  mfcl %~~% L.S,
  log_com_cv %~~% L.S,
  log_mean_cv %~~% L.S,
  gensd %~~% L.S,
  vulsd %~~% L.S,
  avg_ppmr %~~% L.S,
  total_biomass %~~% L.S,
  max_sim %~~% L.S,
  basal %~~% mfcl,
  intermediate %~~% mfcl,
  modularity %~~% mfcl,
  log_mean_cv %~~% mfcl,
  gensd %~~% mfcl,
  vulsd %~~% mfcl,
  avg_ppmr %~~% mfcl,
  max_sim %~~% mfcl,
  modularity %~~% avg_ppmr,
  modularity %~~% basal,
  modularity %~~% intermediate,
  modularity %~~% max_sim,
  modularity %~~% gensd,
  modularity %~~% vulsd,
  gensd %~~% vulsd,
  gensd %~~% avg_ppmr,
  gensd %~~% intermediate,
  gensd %~~% basal,
  gensd %~~% max_sim,
  gensd %~~% total_biomass,
  gensd %~~% log_mean_cv,
  vulsd %~~% total_biomass,
  vulsd %~~% avg_ppmr,
  vulsd %~~% basal,
  vulsd %~~% intermediate,
  vulsd %~~% log_mean_cv,
  log_com_cv %~~% log_mean_cv,
  log_com_cv %~~% basal,
  log_com_cv %~~% total_biomass,
  log_com_cv %~~% max_sim,
  log_com_cv %~~% modularity,
  log_mean_cv %~~% intermediate,
  log_mean_cv %~~% avg_ppmr,
  log_mean_cv %~~% total_biomass,
  log_mean_cv %~~% max_sim,
  avg_ppmr %~~% intermediate,
  avg_ppmr %~~% max_sim,
  avg_ppmr %~~% basal,
  intermediate %~~% total_biomass,
  max_sim %~~% total_biomass,
  max_sim %~~% basal,
  intermediate %~~% basal,
  intermediate %~~% max_sim,
  data=cur_data_orig
)

summary_model_3 <- summary(model_3_spec, conserve=TRUE)

### because standardised coefficients are for some reason not obtained automatedly from psem for responses 'invaded', 'S', and 'L'
### we calculate these a posteriori
### To calculate standardised coefficients for binomial responses we use the following approach developed in :
### https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#standardized-coefficients
### section 2.2.4

glm.model <- glmer(invaded ~ S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_orig, nAGQ = 0)

Beta.glm <- summary(glm.model)$coefficients[2:15, 1]
preds <- predict(glm.model, type = "link")

# Compute sd of error variance using theoretical variances
sd.y.LT <- sqrt(var(preds) + pi^2/3)
# Compute sd of x
sd.x <- apply(cur_data_orig[names(Beta.glm)], 2, sd)

summary_model_3$coefficients[1:14,8] <- round((Beta.glm * sd.x / sd.y.LT), 4)

### for S and L we follow the common standardisation of beta * (sdx/sdy)

sd.y <- apply(cur_data_orig['S'], 2, sd)
predictors <- summary_model_3$coefficients[which(summary_model_3$coefficients$Response == 'S'),]$Predictor
sd.x <- apply(cur_data_orig[as.character(predictors)], 2, sd)

summary_model_3$coefficients[which(summary_model_3$coefficients$Response == 'S'),]$Std.Estimate <- round(summary_model_3$coefficients[which(summary_model_3$coefficients$Response == 'S'),]$Estimate * (sd.x / sd.y), 4)

sd.y <- apply(cur_data_orig['L'], 2, sd)
predictors <- summary_model_3$coefficients[which(summary_model_3$coefficients$Response == 'L'),]$Predictor
sd.x <- apply(cur_data_orig[as.character(predictors)], 2, sd)

summary_model_3$coefficients[which(summary_model_3$coefficients$Response == 'L'),]$Std.Estimate <- round(summary_model_3$coefficients[which(summary_model_3$coefficients$Response == 'L'),]$Estimate * (sd.x / sd.y), 4)

save(summary_model_3, file='sem-model-3.RData')
write.csv(summary_model_3$coefficients, file='sem-model-3-coeffs.csv')



#### Model 4: Deterministic and direct effects of temperature on network / community properties
#### The same as model 3 but including the direct effects of temperature on network properties
#### but not on invasion success

model_4_spec <- psem(
  glmer(invaded ~ S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_orig, nAGQ = 0),
  glmer(S ~ temperature_c + (1|id), data=cur_data_orig, family='poisson', nAGQ = 0),
  glmer(L ~ S + temperature_c + (1|id), data=cur_data_orig, nAGQ=0, family=MASS::negative.binomial(theta=1e5)),
  lmer(L.S ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(mfcl ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(modularity ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(gensd ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(vulsd ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(basal ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(intermediate ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(max_sim ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(log_com_cv ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(log_mean_cv ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(avg_ppmr ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(total_biomass ~ S + temperature_c + (1|id), data=cur_data_orig),
  invaded %~~% temperature_c,
  L.S %~~% L,
  mfcl %~~% L,
  log_com_cv %~~% L,
  basal %~~% L,
  intermediate %~~% L,
  modularity %~~% L,
  gensd %~~% L,
  vulsd %~~% L,
  avg_ppmr %~~% L,
  total_biomass %~~% L,
  max_sim %~~% L,
  basal %~~% L.S,
  intermediate %~~% L.S,
  modularity %~~% L.S,
  mfcl %~~% L.S,
  log_com_cv %~~% L.S,
  gensd %~~% L.S,
  vulsd %~~% L.S,
  avg_ppmr %~~% L.S,
  total_biomass %~~% L.S,
  max_sim %~~% L.S,
  basal %~~% mfcl,
  intermediate %~~% mfcl,
  modularity %~~% mfcl,
  log_com_cv %~~% mfcl,
  total_biomass %~~% mfcl,
  gensd %~~% mfcl,
  vulsd %~~% mfcl,
  avg_ppmr %~~% mfcl,
  max_sim %~~% mfcl,
  modularity %~~% avg_ppmr,
  modularity %~~% basal,
  modularity %~~% intermediate,
  modularity %~~% max_sim,
  modularity %~~% gensd,
  modularity %~~% vulsd,
  modularity %~~% log_mean_cv,
  modularity %~~% total_biomass,
  gensd %~~% vulsd,
  gensd %~~% avg_ppmr,
  gensd %~~% intermediate,
  gensd %~~% basal,
  gensd %~~% max_sim,
  gensd %~~% total_biomass,
  gensd %~~% log_mean_cv,
  vulsd %~~% total_biomass,
  vulsd %~~% avg_ppmr,
  vulsd %~~% basal,
  vulsd %~~% intermediate,
  vulsd %~~% log_mean_cv,
  vulsd %~~% log_com_cv,
  vulsd %~~% max_sim,
  log_com_cv %~~% log_mean_cv,
  log_com_cv %~~% basal,
  log_com_cv %~~% max_sim,
  log_com_cv %~~% intermediate,
  log_mean_cv %~~% intermediate,
  log_mean_cv %~~% avg_ppmr,
  log_mean_cv %~~% total_biomass,
  avg_ppmr %~~% intermediate,
  avg_ppmr %~~% max_sim,
  avg_ppmr %~~% basal,
  avg_ppmr %~~% total_biomass,
  intermediate %~~% total_biomass,
  max_sim %~~% total_biomass,
  max_sim %~~% basal,
  intermediate %~~% basal,
  intermediate %~~% max_sim,
  data=cur_data_orig
)

summary_model_4 <- summary(model_4_spec, conserve=TRUE)

### because standardised coefficients are for some reason not obtained automatedly from psem for responses 'invaded', 'S', and 'L'
### we calculate these a posteriori
### To calculate standardised coefficients for binomial responses we use the following approach developed in :
### https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#standardized-coefficients
### section 2.2.4

glm.model <- glmer(invaded ~ S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_orig, nAGQ = 0)

Beta.glm <- summary(glm.model)$coefficients[2:15, 1]
preds <- predict(glm.model, type = "link")

# Compute sd of error variance using theoretical variances
sd.y.LT <- sqrt(var(preds) + pi^2/3)
# Compute sd of x
sd.x <- apply(cur_data_orig[names(Beta.glm)], 2, sd)

summary_model_4$coefficients[1:14,8] <- round((Beta.glm * sd.x / sd.y.LT), 4)

### for S and L we follow the common standardisation of beta * (sdx/sdy)

sd.y <- apply(cur_data_orig['S'], 2, sd)
predictors <- summary_model_4$coefficients[which(summary_model_4$coefficients$Response == 'S'),]$Predictor
sd.x <- apply(cur_data_orig[as.character(predictors)], 2, sd)

summary_model_4$coefficients[which(summary_model_4$coefficients$Response == 'S'),]$Std.Estimate <- round(summary_model_4$coefficients[which(summary_model_4$coefficients$Response == 'S'),]$Estimate * (sd.x / sd.y), 4)

sd.y <- apply(cur_data_orig['L'], 2, sd)
predictors <- summary_model_4$coefficients[which(summary_model_4$coefficients$Response == 'L'),]$Predictor
sd.x <- apply(cur_data_orig[as.character(predictors)], 2, sd)

summary_model_4$coefficients[which(summary_model_4$coefficients$Response == 'L'),]$Std.Estimate <- round(summary_model_4$coefficients[which(summary_model_4$coefficients$Response == 'L'),]$Estimate * (sd.x / sd.y), 4)

save(summary_model_4, file='sem-model-4.RData')
write.csv(summary_model_4$coefficients, file='sem-model-4-coeffs.csv')




#### Model 5: Instability and productivity are most important.
#### The premise behind this model is that temperature influences invasion success through
#### its effects on ecosystem dynamics and productivity. Thus we model the effect of temperature
#### on community and population stability and on total biomass and then the effects of these
#### on community / network structure and finally on invasion success

model_5_spec <- psem(
  glmer(invaded ~ S + L + L.S + mfcl + modularity + gensd + vulsd + avg_ppmr + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_orig, nAGQ = 0),
  glmer(S ~ log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig, family='poisson', nAGQ = 0),
  glmer(L ~ log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig, nAGQ=0, family=MASS::negative.binomial(theta=1e5)),
  lmer(L.S ~ log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(mfcl ~ log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(modularity ~ log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(gensd ~ log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(vulsd ~ log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(basal ~ log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(intermediate ~ log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(max_sim ~ log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(log_com_cv ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(log_mean_cv ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(avg_ppmr ~ log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(total_biomass ~ temperature_c + (1|id), data=cur_data_orig),
  invaded %~~% temperature_c,
  S %~~% temperature_c,
  L %~~% temperature_c,
  L.S %~~% temperature_c,
  mfcl %~~% temperature_c,
  gensd %~~% temperature_c,
  vulsd %~~% temperature_c,
  modularity %~~% temperature_c,
  basal %~~% temperature_c,
  avg_ppmr %~~% temperature_c,
  log_mean_cv %~~% log_com_cv,
  invaded %~~% log_mean_cv,
  invaded %~~% total_biomass,
  L %~~% S,
  L.S %~~% S,
  mfcl %~~% S,
  modularity %~~% S,
  gensd %~~% S,
  vulsd %~~% S,
  basal %~~% S,
  intermediate %~~% S,
  max_sim %~~% S,
  avg_ppmr %~~% S,
  L.S %~~% L,
  mfcl %~~% L,
  basal %~~% L,
  intermediate %~~% L,
  modularity %~~% L,
  gensd %~~% L,
  vulsd %~~% L,
  avg_ppmr %~~% L,
  max_sim %~~% L,
  basal %~~% L.S,
  intermediate %~~% L.S,
  modularity %~~% L.S,
  mfcl %~~% L.S,
  gensd %~~% L.S,
  vulsd %~~% L.S,
  avg_ppmr %~~% L.S,
  max_sim %~~% L.S,
  basal %~~% mfcl,
  intermediate %~~% mfcl,
  modularity %~~% mfcl,
  gensd %~~% mfcl,
  vulsd %~~% mfcl,
  avg_ppmr %~~% mfcl,
  max_sim %~~% mfcl,
  modularity %~~% avg_ppmr,
  modularity %~~% basal,
  modularity %~~% intermediate,
  modularity %~~% max_sim,
  modularity %~~% gensd,
  modularity %~~% vulsd,
  gensd %~~% vulsd,
  gensd %~~% avg_ppmr,
  gensd %~~% intermediate,
  gensd %~~% basal,
  gensd %~~% max_sim,
  vulsd %~~% avg_ppmr,
  vulsd %~~% basal,
  vulsd %~~% intermediate,
  vulsd %~~% max_sim,
  log_mean_cv %~~% total_biomass,
  avg_ppmr %~~% intermediate,
  avg_ppmr %~~% max_sim,
  avg_ppmr %~~% basal,
  max_sim %~~% basal,
  intermediate %~~% basal,
  intermediate %~~% max_sim,
  data=cur_data_orig
)
 
summary_model_5 <- summary(model_5_spec, conserve=TRUE)

### because standardised coefficients are for some reason not obtained automatedly from psem for responses 'invaded', 'S', and 'L'
### we calculate these a posteriori
### To calculate standardised coefficients for binomial responses we use the following approach developed in :
### https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#standardized-coefficients
### section 2.2.4

glm.model <- glmer(invaded ~ S + L + L.S + mfcl + modularity + gensd + vulsd + avg_ppmr + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_orig, nAGQ = 0)

Beta.glm <- summary(glm.model)$coefficients[2:12, 1]
preds <- predict(glm.model, type = "link")

# Compute sd of error variance using theoretical variances
sd.y.LT <- sqrt(var(preds) + pi^2/3)
# Compute sd of x
sd.x <- apply(cur_data_orig[names(Beta.glm)], 2, sd)

summary_model_5$coefficients[1:11,8] <- round((Beta.glm * sd.x / sd.y.LT), 4)

### for S and L we follow the common standardisation of beta * (sdx/sdy)

sd.y <- apply(cur_data_orig['S'], 2, sd)
predictors <- summary_model_5$coefficients[which(summary_model_5$coefficients$Response == 'S'),]$Predictor
sd.x <- apply(cur_data_orig[as.character(predictors)], 2, sd)

summary_model_5$coefficients[which(summary_model_5$coefficients$Response == 'S'),]$Std.Estimate <- round(summary_model_5$coefficients[which(summary_model_5$coefficients$Response == 'S'),]$Estimate * (sd.x / sd.y), 4)

sd.y <- apply(cur_data_orig['L'], 2, sd)
predictors <- summary_model_5$coefficients[which(summary_model_5$coefficients$Response == 'L'),]$Predictor
sd.x <- apply(cur_data_orig[as.character(predictors)], 2, sd)

summary_model_5$coefficients[which(summary_model_5$coefficients$Response == 'L'),]$Std.Estimate <- round(summary_model_5$coefficients[which(summary_model_5$coefficients$Response == 'L'),]$Estimate * (sd.x / sd.y), 4)

save(summary_model_5, file='sem-model-5.RData')
write.csv(summary_model_5$coefficients, file='sem-model-5-coeffs.csv')


#### Model 6: Community / network properties are affected by biomass, stability and temperature
#### The same as the model above (5) but incorporating the direct effect of temperature on all 
#### community / network variables and on invasion success.

model_6_spec <- psem(
  glmer(invaded ~ temperature_c + S + L + L.S + mfcl + modularity + gensd + vulsd + avg_ppmr + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_orig, nAGQ = 0),
  glmer(S ~ temperature_c + log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig, family='poisson', nAGQ = 0),
  glmer(L ~ temperature_c + log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig, nAGQ=0, family=MASS::negative.binomial(theta=1e5)),
  lmer(L.S ~ temperature_c + log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(mfcl ~ temperature_c + log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(modularity ~ temperature_c + log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(gensd ~ temperature_c + log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(vulsd ~ temperature_c + log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(basal ~ temperature_c + log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(intermediate ~ temperature_c + log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(max_sim ~ temperature_c + log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(log_com_cv ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(log_mean_cv ~ temperature_c + (1|id), data=cur_data_orig),
  lmer(avg_ppmr ~ temperature_c + log_com_cv + log_mean_cv + total_biomass + (1|id), data=cur_data_orig),
  lmer(total_biomass ~ temperature_c + (1|id), data=cur_data_orig),
  log_mean_cv %~~% log_com_cv,
  invaded %~~% log_mean_cv,
  invaded %~~% total_biomass,
  L %~~% S,
  L.S %~~% S,
  mfcl %~~% S,
  modularity %~~% S,
  gensd %~~% S,
  vulsd %~~% S,
  basal %~~% S,
  intermediate %~~% S,
  max_sim %~~% S,
  avg_ppmr %~~% S,
  L.S %~~% L,
  mfcl %~~% L,
  basal %~~% L,
  intermediate %~~% L,
  modularity %~~% L,
  gensd %~~% L,
  vulsd %~~% L,
  avg_ppmr %~~% L,
  max_sim %~~% L,
  basal %~~% L.S,
  intermediate %~~% L.S,
  modularity %~~% L.S,
  mfcl %~~% L.S,
  gensd %~~% L.S,
  vulsd %~~% L.S,
  avg_ppmr %~~% L.S,
  max_sim %~~% L.S,
  basal %~~% mfcl,
  intermediate %~~% mfcl,
  modularity %~~% mfcl,
  gensd %~~% mfcl,
  vulsd %~~% mfcl,
  avg_ppmr %~~% mfcl,
  max_sim %~~% mfcl,
  modularity %~~% avg_ppmr,
  modularity %~~% basal,
  modularity %~~% intermediate,
  modularity %~~% max_sim,
  modularity %~~% gensd,
  modularity %~~% vulsd,
  gensd %~~% vulsd,
  gensd %~~% avg_ppmr,
  gensd %~~% intermediate,
  gensd %~~% basal,
  gensd %~~% max_sim,
  vulsd %~~% avg_ppmr,
  vulsd %~~% basal,
  vulsd %~~% intermediate,
  vulsd %~~% max_sim,
  log_mean_cv %~~% total_biomass,
  avg_ppmr %~~% intermediate,
  avg_ppmr %~~% max_sim,
  avg_ppmr %~~% basal,
  max_sim %~~% basal,
  intermediate %~~% basal,
  intermediate %~~% max_sim,
  data=cur_data_orig
)

summary_model_6 <- summary(model_6_spec, conserve=TRUE)

### because standardised coefficients are for some reason not obtained automatedly from psem for responses 'invaded', 'S', and 'L'
### we calculate these a posteriori
### To calculate standardised coefficients for binomial responses we use the following approach developed in :
### https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#standardized-coefficients
### section 2.2.4

glm.model <- glmer(invaded ~ temperature_c + S + L + L.S + mfcl + modularity + gensd + vulsd + avg_ppmr + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_orig, nAGQ = 0)

Beta.glm <- summary(glm.model)$coefficients[2:13, 1]
preds <- predict(glm.model, type = "link")

# Compute sd of error variance using theoretical variances
sd.y.LT <- sqrt(var(preds) + pi^2/3)
# Compute sd of x
sd.x <- apply(cur_data_orig[names(Beta.glm)], 2, sd)

summary_model_6$coefficients[1:12,8] <- round((Beta.glm * sd.x / sd.y.LT), 4)

### for S and L we follow the common standardisation of beta * (sdx/sdy)

sd.y <- apply(cur_data_orig['S'], 2, sd)
predictors <- summary_model_6$coefficients[which(summary_model_6$coefficients$Response == 'S'),]$Predictor
sd.x <- apply(cur_data_orig[as.character(predictors)], 2, sd)

summary_model_6$coefficients[which(summary_model_6$coefficients$Response == 'S'),]$Std.Estimate <- round(summary_model_6$coefficients[which(summary_model_6$coefficients$Response == 'S'),]$Estimate * (sd.x / sd.y), 4)

sd.y <- apply(cur_data_orig['L'], 2, sd)
predictors <- summary_model_6$coefficients[which(summary_model_6$coefficients$Response == 'L'),]$Predictor
sd.x <- apply(cur_data_orig[as.character(predictors)], 2, sd)

summary_model_6$coefficients[which(summary_model_6$coefficients$Response == 'L'),]$Std.Estimate <- round(summary_model_6$coefficients[which(summary_model_6$coefficients$Response == 'L'),]$Estimate * (sd.x / sd.y), 4)

save(summary_model_6, file='sem-model-6.RData')
write.csv(summary_model_6$coefficients, file='sem-model-6-coeffs.csv')

