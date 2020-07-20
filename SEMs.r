

require(lme4)
require(piecewiseSEM)
require(vegan)

load('data-for-sems.rda')

cur_data_norm <- decostand(cur_data_orig[c('temperature_c', 'S', 'L', 'C', 'L.S', 'mfcl', 'modularity', 'gensd', 'vulsd', 'log_com_cv', 'log_mean_cv', 'avg_ppmr', 'total_biomass', 'basal', 'intermediate', 'top', 'max_sim')], 'standardize')
cur_data_norm$id <- cur_data_orig$id
cur_data_norm$invaded <- cur_data_orig$invaded


#### This is the original figures for the SEMs as they appeared in the first draft of the paper
sem_invaded <- psem(
  glmer(invaded ~ temperature_c + S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_norm, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5))),
  lmer(S ~ temperature_c + (1|id), data=cur_data_norm),
  lmer(L ~ S + temperature_c + (1|id), data=cur_data_norm),
  lmer(L.S ~ S + temperature_c + (1|id), data=cur_data_norm),
  lmer(mfcl ~ S + temperature_c + (1|id), data=cur_data_norm),
  lmer(modularity ~ S + temperature_c + (1|id), data=cur_data_norm),
  lmer(gensd ~ S + temperature_c + (1|id), data=cur_data_norm),
  lmer(vulsd ~ S + temperature_c + (1|id), data=cur_data_norm),
  lmer(log_com_cv ~ S + temperature_c + (1|id), data=cur_data_norm),
  lmer(log_mean_cv ~ S + temperature_c + (1|id), data=cur_data_norm),
  lmer(avg_ppmr ~ S + temperature_c + (1|id), data=cur_data_norm),
  lmer(total_biomass ~ S + temperature_c + (1|id), data=cur_data_norm),
  lmer(basal ~ S + temperature_c + (1|id), data=cur_data_norm),
  lmer(intermediate ~ S + temperature_c + (1|id), data=cur_data_norm),
  lmer(top ~ S + temperature_c + (1|id), data=cur_data_norm),
  lmer(max_sim ~ S + temperature_c + (1|id), data=cur_data_norm),
  data=cur_data_norm
)

sem_invaded <- coefs(sem_invaded)
fisherC(sem_invaded)

#### This is what happens when variables are not normalised
#### since S and L become count data here we need to model them as poisson and negative binomial respectively

# This might be used for checking the model fit

cur_data_orig$L_boxcox <- T_box
my_model <- glmer(L ~ S + temperature_c + (1|id), data=cur_data_orig, nAGQ=0, family=MASS::negative.binomial(theta=1.75))
your_model <- glmer(L ~ S + temperature_c + (1|id), data=cur_data_orig, nAGQ=0, family='poisson')
your_model <- lmer(L_boxcox ~ S + temperature_c + (1|id), data=cur_data_orig)


your_model <- glmmTMB(L ~ S + temperature_c + (1|id), data=cur_data_orig, family=compois)

#### The model using compois seems to not suffer from under or overdispersion... 
#### let's see what happens when we split the range

require(blmeco)
dispersion_glmer(your_model)

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(your_model)

require(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = your_model, n=200, plot=F, refit = F)
testResiduals(simulationOutput)
plot(simulationOutput)
testDispersion(simulationOutput, alternative = 'less')


par(mfrow=c(2,2))
qqnorm(resid(your_model), main="normal qq-plot, residuals")
qqline(resid(your_model))

qqnorm(ranef(your_model)$id[,1])
qqline(ranef(your_model)$id[,1])

plot(fitted(your_model), resid(your_model)) #residuals vs fitted
abline(h=0)


fitted <- fitted(your_model)    #fitted vs observed
plot(fitted, jitter(cur_data_orig$L,0.1))
abline(0,1)


sem_invaded <- psem(
  glmer(invaded ~ temperature_c + S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=cur_data_orig, nAGQ = 0),
  glmer(S ~ temperature_c + (1|id), data=cur_data_orig, family='poisson', nAGQ = 0),
  glmer.nb(L ~ S + temperature_c + (1|id), data=cur_data_orig, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))),
  lmer(L.S ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(mfcl ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(modularity ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(gensd ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(vulsd ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(log_com_cv ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(log_mean_cv ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(avg_ppmr ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(total_biomass ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(basal ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(intermediate ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(top ~ S + temperature_c + (1|id), data=cur_data_orig),
  lmer(max_sim ~ S + temperature_c + (1|id), data=cur_data_orig),
  data=cur_data_orig
)


sem_invaded <- coefs(sem_invaded)
fisherC(sem_invaded)


#### Let's split the range between 5 to 25 and 26 to 40
#### and also include the other hypothesised relationships
data_first_half <- subset(cur_data_orig, (temperature_c >=5 & temperature_c <=25))
data_first_half$id <- as.factor(data_first_half$id)
sem_invaded_list <- psem(
  glmer(invaded ~ temperature_c + S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=data_first_half, nAGQ = 0),
  glmer(S ~ temperature_c + (1|id), data=data_first_half, nAGQ=0, family=poisson),
  glmer(L ~ S + temperature_c + (1|id), data=data_first_half, nAGQ=0, family=MASS::negative.binomial(theta=1e5)),
  lmer(L.S ~ S + L + temperature_c + (1|id), data=data_first_half),
  lmer(basal ~ S + L + temperature_c + (1|id), data=data_first_half),
  lmer(intermediate ~ S + L + temperature_c + (1|id), data=data_first_half),
  lmer(mfcl ~ S + L + temperature_c + (1|id), data=data_first_half),
  lmer(modularity ~ S + L + temperature_c + (1|id), data=data_first_half),
  lmer(gensd ~ S + L + temperature_c + (1|id), data=data_first_half),
  lmer(vulsd ~ S + L + temperature_c + (1|id), data=data_first_half),
  lmer(log_com_cv ~ S + L.S + modularity + mfcl + temperature_c + (1|id), data=data_first_half),
  lmer(log_mean_cv ~ S + L.S + modularity + mfcl + temperature_c + (1|id), data=data_first_half),
  lmer(avg_ppmr ~ S + L + temperature_c + (1|id), data=data_first_half),
  lmer(total_biomass ~ S + L + basal + intermediate + mfcl + temperature_c + (1|id), data=data_first_half),
  lmer(max_sim ~ S + L + temperature_c + (1|id), data=data_first_half),
  basal %~~% L.S,
  intermediate %~~% L.S,
  mfcl %~~% L.S,
  modularity %~~% L.S,
  gensd %~~% L.S,
  vulsd %~~% L.S,
  avg_ppmr %~~% L.S,
  total_biomass %~~% L.S,
  max_sim %~~% L.S,
  basal %~~% mfcl,
  intermediate %~~% mfcl,
  modularity %~~% mfcl,
  gensd %~~% mfcl,
  vulsd %~~% mfcl,
  avg_ppmr %~~% mfcl,
  max_sim %~~% mfcl,
  intermediate %~~% basal,
  gensd %~~% basal,
  vulsd %~~% basal,
  log_com_cv %~~% basal,
  log_mean_cv %~~% basal,
  avg_ppmr  %~~% basal,
  max_sim %~~% basal,
  modularity %~~% max_sim,
  modularity %~~% avg_ppmr,
  modularity %~~% vulsd,
  modularity %~~% total_biomass,
  modularity %~~% basal,
  gensd %~~% vulsd,
  gensd %~~% avg_ppmr,
  gensd %~~% intermediate,
  gensd %~~% max_sim,
  gensd %~~% log_mean_cv,
  gensd %~~% total_biomass,
  vulsd %~~% avg_ppmr,
  vulsd %~~% intermediate,
  vulsd %~~% max_sim,
  vulsd %~~% total_biomass,
  vulsd %~~% log_mean_cv,
  log_com_cv %~~% log_mean_cv,
  log_com_cv %~~% intermediate,
  log_com_cv %~~% gensd,
  log_com_cv %~~% vulsd,
  log_com_cv %~~% max_sim,
  log_com_cv %~~% avg_ppmr,
  log_mean_cv %~~% intermediate,
  log_mean_cv %~~% max_sim,
  log_mean_cv %~~% avg_ppmr,
  log_mean_cv %~~% total_biomass,
  avg_ppmr %~~% intermediate,
  avg_ppmr %~~% max_sim,
  data=data_first_half
)

sem_invaded <- summary(sem_invaded_list)

## the p-value for the Fisher C for this model is way above 0.05 (=0.478), indicating this model is a good fit
sem_invaded$Cstat

### because standardised coefficients are for some reason not obtained automatedly from psem for responses 'invaded', 'S', and 'L'
### we calculate these a posteriori
### To calculate standardised coefficients for binomial responses we use the following approach developed in :
### https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#standardized-coefficients
### section 2.2.4

glm.model <- glmer(invaded ~ temperature_c + S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=data_first_half, nAGQ = 0)
Beta.glm <- summary(glm.model)$coefficients[2:16, 1]
preds <- predict(glm.model, type = "link")

# Compute sd of error variance using theoretical variances
sd.y.LT <- sqrt(var(preds) + pi^2/3)
# Compute sd of x
sd.x <- apply(data_first_half[names(Beta.glm)], 2, sd)

sem_invaded$coefficients[1:15,8] <- round((Beta.glm * sd.x / sd.y.LT), 4)

### for S and L we follow the common standardisation of beta * (sdx/sdy)

sd.y <- apply(data_first_half['S'], 2, sd)
predictors <- sem_invaded$coefficients[which(sem_invaded$coefficients$Response == 'S'),]$Predictor
sd.x <- apply(data_first_half[as.character(predictors)], 2, sd)

sem_invaded$coefficients[which(sem_invaded$coefficients$Response == 'S'),]$Std.Estimate <- round(sem_invaded$coefficients[which(sem_invaded$coefficients$Response == 'S'),]$Estimate * (sd.x / sd.y), 4)

sd.y <- apply(data_first_half['L'], 2, sd)
predictors <- sem_invaded$coefficients[which(sem_invaded$coefficients$Response == 'L'),]$Predictor
sd.x <- apply(data_first_half[as.character(predictors)], 2, sd)

sem_invaded$coefficients[which(sem_invaded$coefficients$Response == 'L'),]$Std.Estimate <- round(sem_invaded$coefficients[which(sem_invaded$coefficients$Response == 'L'),]$Estimate * (sd.x / sd.y), 4)

write.csv(sem_invaded$coefficients, file = 'coefficients-first-half.csv')
write.csv(sem_invaded$R2, file = 'r-squared-first-half.csv')


#### We repeat the same analysis for the second have of the temperature range (from 26 to 40 degrees)
data_second_half <- subset(cur_data_orig, (temperature_c > 25))
data_second_half$id <- as.factor(data_second_half$id)
sem_invaded_list <- psem(
  glmer(invaded ~ temperature_c + S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=data_second_half, nAGQ = 0),
  glmer(S ~ temperature_c + (1|id), data=data_second_half, nAGQ=0, family=poisson),
  glmer(L ~ S + temperature_c + (1|id), data=data_second_half, nAGQ=0, family=MASS::negative.binomial(theta=1e5)),
  lmer(L.S ~ S + L + temperature_c + (1|id), data=data_second_half),
  lmer(basal ~ S + L + temperature_c + (1|id), data=data_second_half),
  lmer(intermediate ~ S + L + temperature_c + (1|id), data=data_second_half),
  lmer(mfcl ~ S + L + temperature_c + (1|id), data=data_second_half),
  lmer(modularity ~ S + L + temperature_c + (1|id), data=data_second_half),
  lmer(gensd ~ S + L + temperature_c + (1|id), data=data_second_half),
  lmer(vulsd ~ S + L + temperature_c + (1|id), data=data_second_half),
  lmer(log_com_cv ~ S + L.S + modularity + mfcl + temperature_c + (1|id), data=data_second_half),
  lmer(log_mean_cv ~ S + L.S + modularity + mfcl + temperature_c + (1|id), data=data_second_half),
  lmer(avg_ppmr ~ S + L + temperature_c + (1|id), data=data_second_half),
  lmer(total_biomass ~ S + L + basal + intermediate + mfcl + temperature_c + (1|id), data=data_second_half),
  lmer(max_sim ~ S + L + temperature_c + (1|id), data=data_second_half),
  basal %~~% L.S,
  intermediate %~~% L.S,
  mfcl %~~% L.S,
  modularity %~~% L.S,
  gensd %~~% L.S,
  vulsd %~~% L.S,
  avg_ppmr %~~% L.S,
  max_sim %~~% L.S,
  intermediate %~~% mfcl,
  modularity %~~% mfcl,
  gensd %~~% mfcl,
  vulsd %~~% mfcl,
  avg_ppmr %~~% mfcl,
  max_sim %~~% mfcl,
  intermediate %~~% basal,
  gensd %~~% basal,
  vulsd %~~% basal,
  log_com_cv %~~% gensd,
  log_mean_cv %~~% gensd,
  avg_ppmr  %~~% basal,
  max_sim %~~% basal,
  modularity %~~% max_sim,
  modularity %~~% avg_ppmr,
  modularity %~~% vulsd,
  modularity %~~% gensd,
  modularity %~~% basal,
  modularity %~~% intermediate,
  gensd %~~% vulsd,
  gensd %~~% intermediate,
  gensd %~~% max_sim,
  gensd %~~% total_biomass,
  vulsd %~~% avg_ppmr,
  vulsd %~~% intermediate,
  vulsd %~~% max_sim,
  vulsd %~~% log_mean_cv,
  log_com_cv %~~% log_mean_cv,
  avg_ppmr %~~% intermediate,
  avg_ppmr %~~% max_sim,
  avg_ppmr %~~% total_biomass,
  data=data_second_half
)

sem_invaded <- summary(sem_invaded_list)

## the p-value for the Fisher C for this model is way above 0.05 (=0.144), indicating this model is a good fit
sem_invaded$Cstat

### because standardised coefficients are for some reason not obtained automatedly from psem for responses 'invaded', 'S', and 'L'
### we calculate these a posteriori
### To calculate standardised coefficients for binomial responses we use the following approach developed in :
### https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#standardized-coefficients
### section 2.2.4

glm.model <- glmer(invaded ~ temperature_c + S + L + L.S + mfcl + modularity + gensd + vulsd + log_com_cv + log_mean_cv + avg_ppmr + total_biomass + basal + intermediate + max_sim + (1|id), family="binomial", data=data_second_half, nAGQ = 0)
Beta.glm <- summary(glm.model)$coefficients[2:16, 1]
preds <- predict(glm.model, type = "link")

# Compute sd of error variance using theoretical variances
sd.y.LT <- sqrt(var(preds) + pi^2/3)
# Compute sd of x
sd.x <- apply(data_second_half[names(Beta.glm)], 2, sd)

sem_invaded$coefficients[1:15,8] <- round((Beta.glm * sd.x / sd.y.LT), 4)

### for S and L we follow the common standardisation of beta * (sdx/sdy)

sd.y <- apply(data_second_half['S'], 2, sd)
predictors <- sem_invaded$coefficients[which(sem_invaded$coefficients$Response == 'S'),]$Predictor
sd.x <- apply(data_second_half[as.character(predictors)], 2, sd)

sem_invaded$coefficients[which(sem_invaded$coefficients$Response == 'S'),]$Std.Estimate <- round(sem_invaded$coefficients[which(sem_invaded$coefficients$Response == 'S'),]$Estimate * (sd.x / sd.y), 4)

sd.y <- apply(data_second_half['L'], 2, sd)
predictors <- sem_invaded$coefficients[which(sem_invaded$coefficients$Response == 'L'),]$Predictor
sd.x <- apply(data_second_half[as.character(predictors)], 2, sd)

sem_invaded$coefficients[which(sem_invaded$coefficients$Response == 'L'),]$Std.Estimate <- round(sem_invaded$coefficients[which(sem_invaded$coefficients$Response == 'L'),]$Estimate * (sd.x / sd.y), 4)

write.csv(sem_invaded$coefficients, file = 'coefficients-second-half.csv')
write.csv(sem_invaded$R2, file = 'r-squared-second-half.csv')




######## This code is to look at the residuals and dispersion of models to see if we are happy with these...

your_model <- glmer(L ~ S + temperature_c + (1|id), data=data_first_half, nAGQ=0, family=MASS::negative.binomial(theta=1e5))

simulationOutput <- simulateResiduals(fittedModel = your_model, n=2000, plot=F, refit = F)
testResiduals(simulationOutput)
plot(simulationOutput)
testDispersion(simulationOutput, alternative = 'less')

your_model <- glmer(S ~ temperature_c + (1|id), data=data_first_half, nAGQ=0, family=poisson)
dispersion_glmer(your_model)

your_model <- lmer(L.S ~S + temperature_c + (1|id), data=data_first_half)
dispersion_glmer(your_model)

par(mfrow=c(2,2))
qqnorm(resid(your_model), main="normal qq-plot, residuals")
qqline(resid(your_model))

qqnorm(ranef(your_model)$id[,1])
qqline(ranef(your_model)$id[,1])

plot(fitted(your_model), resid(your_model)) #residuals vs fitted
abline(h=0)

fitted <- fitted(your_model)    #fitted vs observed
plot(fitted, jitter(data_first_half$L.S,0.1))
abline(0,1)



