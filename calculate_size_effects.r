
## ---------------------------
##
## Script name: calculate_size_effects.r
##
## Purpose of script: This script takes the output tables from the temperature dependent simulation
## experiments and extracts and analyses the relevant information to investigate the size of the 
## effects of invasions on complex food webs. It additionally compares the effects of invasions
## on communities on which the invasive species was successfully established vs. those resulting
## from unsucessful invasion attempts, across a large temperature gradient (from 0 to 40 degrees)
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

#### Required libraries and functions:
require(reshape2)
require(RColorBrewer)
require(lme4)
require(piecewiseSEM)

##### The following function is defined to calculate some summary statistics over the data tables
##### It was written and provided by Arnaud Sentis
getPropertiesSummary <- function(x){
  require(plyr)
  str(x)
  x <- ddply(x, .(temperature_c,eutrophication), summarise,
             N = length(na.omit(S)),
             S.mean= mean(S, na.rm=T),
             S.sd=sd(S, na.rm=T),
             S.se=sd(S, na.rm=T)/sqrt(length(na.omit(S))),
             S.CI=S.se*1.96,
             L.mean= mean(L, na.rm=T),
             L.sd=sd(L, na.rm=T),
             L.se=sd(L, na.rm=T)/sqrt(length(na.omit(L))),
             L.CI=L.se*1.96,
             L.S.mean= mean(L.S, na.rm=T),
             L.S.sd=sd(L.S, na.rm=T),
             L.S.se=sd(L.S, na.rm=T)/sqrt(length(na.omit(L.S))),
             L.S.CI=L.S.se*1.96,
             C.mean= mean(C, na.rm=T),
             C.sd=sd(C, na.rm=T),
             C.se=sd(C, na.rm=T)/sqrt(length(na.omit(C))),
             C.CI=C.se*1.96,
             total_biomass.mean= mean((total_biomass), na.rm=T),
             total_biomass.se=sd((total_biomass), na.rm=T)/sqrt(length(na.omit((total_biomass)))),
             total_biomass.sd=sd((total_biomass), na.rm=T),
             total_biomass.CI=total_biomass.se*1.96,
             modularity.mean= mean(modularity, na.rm=T),
             modularity.se=sd(modularity, na.rm=T)/sqrt(length(na.omit(modularity))),
             modularity.sd=sd(modularity, na.rm=T),
             modularity.CI=modularity.se*1.96,
             com_cv.mean= mean(com_cv, na.rm=T),
             com_cv.se=sd(com_cv, na.rm=T)/sqrt(length(na.omit(com_cv))),
             com_cv.sd=sd(com_cv, na.rm=T),
             com_cv.CI=com_cv.se*1.96,
             com_cv_last_100.mean= mean(com_cv_last_100, na.rm=T),
             com_cv_last_100.se=sd(com_cv_last_100, na.rm=T)/sqrt(length(na.omit(com_cv_last_100))),
             com_cv_last_100.sd=sd(com_cv_last_100, na.rm=T),
             com_cv_last_100.CI=com_cv_last_100.se*1.96,
             mean_cv.mean= mean((mean_cv), na.rm=T),
             mean_cv.se=sd((mean_cv), na.rm=T)/sqrt(length(na.omit((mean_cv)))),
             mean_cv.sd=sd((mean_cv), na.rm=T),
             mean_cv.CI=mean_cv.se*1.96,
             mean_cv_last_100.mean= mean(mean_cv_last_100, na.rm=T),
             mean_cv_last_100.se=sd(mean_cv_last_100, na.rm=T)/sqrt(length(na.omit(com_cv_last_100))),
             mean_cv_last_100.sd=sd(mean_cv_last_100, na.rm=T),
             mean_cv_last_100.CI=mean_cv_last_100.se*1.96,
             log_com_cv.mean= mean(log_com_cv, na.rm=T),
             log_com_cv.se=sd(log_com_cv, na.rm=T)/sqrt(length(na.omit(log_com_cv))),
             log_com_cv.sd=sd(log_com_cv, na.rm=T),
             log_com_cv.CI=log_com_cv.se*1.96,
             log_mean_cv.mean= mean((log_mean_cv), na.rm=T),
             log_mean_cv.se=sd((log_mean_cv), na.rm=T)/sqrt(length(na.omit((log_mean_cv)))),
             log_mean_cv.sd=sd((log_mean_cv), na.rm=T),
             log_mean_cv.CI=log_mean_cv.se*1.96,
             mfcl.mean= mean(mfcl, na.rm=T),
             mfcl.se=sd(mfcl, na.rm=T)/sqrt(length(na.omit(mfcl))),
             mfcl.sd=sd(mfcl, na.rm=T),
             mfcl.CI=mfcl.se*1.96,
             max_sim.mean= mean(max_sim, na.rm=T),
             max_sim.se=sd(max_sim, na.rm=T)/sqrt(length(na.omit(max_sim))),
             max_sim.sd=sd(max_sim, na.rm=T),
             max_sim.CI=max_sim.se*1.96,
             inv_tp.mean= mean(inv_tp, na.rm=T),
             inv_tp.se=sd(inv_tp, na.rm=T)/sqrt(length(na.omit(inv_tp))),
             inv_tp.sd=sd(inv_tp, na.rm=T),
             inv_tp.CI=inv_tp.se*1.96,
             vulsd.mean=mean(vulsd, na.rm=T),
             vulsd.se=sd(vulsd, na.rm=T)/sqrt(length(na.omit(vulsd))),
             vulsd.sd=sd(vulsd, na.rm=T),
             vulsd.CI=vulsd.se*1.96,
             gensd.mean =mean(gensd, na.rm=T),
             gensd.se=sd(gensd, na.rm=T)/sqrt(length(na.omit(gensd))),
             gensd.sd=sd(gensd, na.rm=T),
             gensd.CI=gensd.se*1.96,
             avg_ppmr.mean=mean(avg_ppmr, na.rm=T),
             avg_ppmr.se=sd(avg_ppmr, na.rm=T)/sqrt(length(na.omit(avg_ppmr))),
             avg_ppmr.sd=sd(avg_ppmr, na.rm=T),
             avg_ppmr.CI=avg_ppmr.se*1.96,
             basal.mean=mean(basal, na.rm=T),
             basal.se=sd(basal, na.rm=T)/sqrt(length(na.omit(basal))),
             basal.sd=sd(basal, na.rm=T),
             basal.CI=basal.se*1.96,
             intermediate.mean=mean(intermediate, na.rm=T),
             intermediate.se=sd(intermediate, na.rm=T)/sqrt(length(na.omit(intermediate))),
             intermediate.sd=sd(intermediate, na.rm=T),
             intermediate.CI=intermediate.se*1.96,
             top.mean=mean(top, na.rm=T),
             top.se=sd(top, na.rm=T)/sqrt(length(na.omit(top))),
             top.sd=sd(top, na.rm=T),
             top.CI=top.se*1.96,
             mean_bs.mean=mean((mean_bs), na.rm=T),
             mean_bs.se=sd((mean_bs), na.rm=T)/sqrt(length(na.omit((mean_bs)))),
             mean_bs.sd=sd((mean_bs), na.rm=T),
             mean_bs.CI=mean_bs.se*1.96
  )
  
  head(x)
  return(x)
}

############### HERE STARTS THE CODE FOR THE ANALYSIS #####################

#### Here we load the output from the simulations
load('simulation-outputs.rda')

### Information on the data table:
##non-self-explanatory variables################
# id: network ID number
# type: type of network: 0 is the niche model before any simulation
#                        1 is the niche model after the transiant dynamics (xxxx days)
#                        3 is the networ after the invasion
# exp_type: type of the experiment: "biomass","body_size"
# exp_val: value of the experimental variable (e.g. body size of the invader)
# S: species richness
# C: connectance
# vulsd: vulnerability sd
# gensd: generality sd
# mfcl: mid food chain length
# avg_ppmr: average predator prey mass ratio
# basal,intermediate,top: proportions of basal, intemediate and top species
# max_sim: maximum similarity
# inv_biom: biomass of the invador at the end of the simulations        
# inv_tp: invador trophic position: 1 is basal, 2 intermediate, 3 top predator        
# inv_gen: invador generality
# inv_vul: invador generaility

#### This is where we obtain the output data for the invasion experiments
cur_data <- subset(networks, (exp_type=='neutral' & eutrophication == 0)) # & inv_extinct == FALSE))

#### Order the data frame to be more tidy
cur_data <- cur_data[with(cur_data, order(id, temperature, eutrophication, exp_val)),]

#### Some data type conversions to make sure we manipulate data correctly
cur_data$inv_extinct <- as.logical(cur_data$inv_extinct)
cur_data$inv_extinct <- as.numeric(cur_data$inv_extinct)

#### Since we want to look at the CV (i.e. variability, the inverse of stability) log10-transformed, here
#### we perform this transformation, both for the community and population stability measures
#### the 'abs' is the same as the log of the inverse, so we get stability
cur_data$log_mean_cv <- abs(log10(cur_data$mean_cv_last_100)) ### this is the same as the log of the inverse
cur_data$log_mean_cv[cur_data$log_mean_cv == Inf] <- 0

cur_data$log_com_cv <- abs(log10(cur_data$com_cv_last_100))
cur_data$log_com_cv[cur_data$log_com_cv == Inf] <- 0

##### And we obtain the corresponding networks before invasion so we can compare the properties before vs. after invasion
##### These are the original networks (i.e. before invasion) after transient dynamics (select type 1)
cur_data_orig <- subset(networks, (type == 1 & eutrophication == 0)) ##should contain 41(temp values)*3(eutrophication values)*140(unique networks) = 17220 rows

#### We order the data
cur_data_orig <- cur_data_orig[with(cur_data_orig, order(id, temperature, eutrophication)),]
#### Look at the data types...
str(cur_data_orig)
#### .... and based on that make some data type transformations
cur_data_orig$inv_extinct <- as.logical(cur_data_orig$inv_extinct)
cur_data_orig$inv_extinct <- as.numeric(cur_data_orig$inv_extinct)

#### Convert stability measures to log10-transformed as above
cur_data_orig$log_mean_cv <- abs(log10(cur_data_orig$mean_cv_last_100)) ### this is the same as the log of the inverse
cur_data_orig$log_mean_cv[cur_data_orig$log_mean_cv == Inf] <- 0

cur_data_orig$log_com_cv <- abs(log10(cur_data_orig$com_cv_last_100))
cur_data_orig$log_com_cv[cur_data_orig$log_com_cv == Inf] <- 0

########################### These make sure data are consistent  ###########################


#### To be able to compare like with like to obtain before vs. after ratios we need to identify the networks 
#### from the original dataset that match those in the invaded dataset. Here is how you do that:

###### we now know the invaded and non-invaded communities are very different.
###### This could have been because they were different in the first place.
###### Let's look at that and then at the size effects plots for each.

cur_data$unique_id <- interaction(cur_data$id, cur_data$temperature_c)
cur_data_orig$unique_id <- interaction(cur_data_orig$id, cur_data_orig$temperature_c)

#### We separate those where the ivnasion was successfull vs. where it wasn't
cur_data_invaded <- subset(cur_data, (inv_extinct == FALSE))
cur_data_not_invaded <- subset(cur_data, (inv_extinct == TRUE))

#### Once we know which networks were invaded and which ones weren't, we can use the unique ID to match them
#### hence identifying which ones to compare with which
cur_data_orig_not_invaded <- subset(cur_data_orig, unique_id %in% cur_data_not_invaded$unique_id)
cur_data_orig_invaded <- subset(cur_data_orig, unique_id %in% cur_data_invaded$unique_id)

#### If the above is correct, the line below whould yield 0
length(intersect(cur_data_orig_not_invaded$unique_id, cur_data_orig_invaded$unique_id))

#### Here we specify the network / community properties we are interested in
props <- c("S", "L", "L.S", "C", "mean_bs", "vulsd", "gensd", "mfcl", "ChnSD", 
           "ChnNo", "avg_ppmr", "basal", "intermediate", "top", "ca", "components",
           "modularity", "n_modules", "max_sim", "total_biomass", "log_mean_cv", "log_com_cv")

data_effects <- cur_data_invaded
for(p in props){
  eval(parse(text=paste0('data_effects$',p,' <- cur_data_invaded[match(cur_data_orig_invaded$unique_id, cur_data_invaded$unique_id),]$', p, '/cur_data_orig_invaded$',p)))
}

data_effects[(data_effects == Inf)] <- 0
data_effects[is.na(data_effects)] <- 0


data_effects_not_invaded <- cur_data_not_invaded
for(p in props){
    eval(parse(text=paste0('data_effects_not_invaded$',p,' <- cur_data_not_invaded[match(cur_data_orig_not_invaded$unique_id, cur_data_not_invaded$unique_id),]$', p, '/cur_data_orig_not_invaded$',p)))  
}

data_effects_not_invaded[(data_effects_not_invaded == Inf)] <- 0
data_effects_not_invaded[is.na(data_effects_not_invaded)] <- 0


##### Once we have had calculated the effect of successful and unsuccessful invaders on the recipient communities 
##### (quantified as the ratio of the value of each property before vs. after the invasion event), we can then check 
##### whether the effects of the invasion vs. non-invasion are different.
##### 
differences <- NULL
for(p in props){
  for(t in seq(0,40,by=5)){
    group_1 <- eval(parse(text=paste0('subset(data_effects, temperature_c == t)$',p)))
    group_2 <- eval(parse(text=paste0('subset(data_effects_not_invaded, temperature_c == t)$',p)))
    test_res <- wilcox.test(group_1, group_2)
    n1 <- length(group_1)
    n2 <- length(group_2)
    W <- test_res$statistic
    effect_size <- W/(n1*n2)
    
    cur_out <- data.frame(property = p, temperature = t, median_1=round(median(group_1),2), median_2=round(median(group_2),2), U=W, effect_size=round(effect_size,2), n1, n2, p.val=round(test_res$p.value,3))
    if(is.null(differences)){
      differences <- cur_out
    }else{
      differences <- rbind(differences, cur_out)
    }
    
  }
}

differences[which(differences$p.val < 0.05),]

write.csv(differences, file = 'diff-in-effects-ratios.csv')


## Here we calculate the mean and standard error / confident intervals of the effects
## and perform some transformations on the data format to facilitate plotting
data_effects_sum <- getPropertiesSummary(data_effects)
names(data_effects_sum)[12:15] <- c("LxS.mean", "LxS.sd", "LxS.se", "LxS.CI")

data_effects_sum_reformat <- melt(data_effects_sum, id.vars = 'temperature_c', measure.vars = names(data_effects_sum)[grep('\\.mean', names(data_effects_sum))])
data_effects_sum_reformat$variable <- unlist(lapply(as.character(data_effects_sum_reformat$variable), function(x) {strsplit(x, '\\.')[[1]][1]}))

names(data_effects_sum_reformat) <- c('temperature', 'property', 'mean')

data_effects_sum_reformat$CI <- melt(data_effects_sum, id.vars = 'temperature_c', measure.vars = names(data_effects_sum)[grep('.CI', names(data_effects_sum))])$value
data_effects_sum_reformat$se <- melt(data_effects_sum, id.vars = 'temperature_c', measure.vars = names(data_effects_sum)[grep('\\.se', names(data_effects_sum))])$value


#### Before deciding which data to plot on the effect size figures we can do a quick SEM to see
#### which of these properties are differentially affected by the invasion across the temperature range

#### we can at the same time determine whether the differences in effects due to temperature changes
#### are mediated by the number of species using SEMs as for the invasion success

data_effects$original_S <- cur_data_orig_invaded$S
data_effects_norm <- as.data.frame(scale(data_effects[c('temperature_c', 'original_S', "S", "L", "L.S", "C", "mean_bs", "vulsd", "gensd", "mfcl", "ChnSD", "ChnNo", "avg_ppmr", "basal", "intermediate", "top", "ca", "components",
                                                        "modularity", "n_modules", "max_sim", "total_biomass", "log_com_cv", "log_mean_cv", 'inv_tp', 'inv_bs', 'inv_biom', 'inv_gen_norm', 'inv_vul_norm')]))

data_effects_norm <- data_effects
data_effects_norm$id <- data_effects$id


modlist <- psem(
  lmer(original_S ~ temperature_c + (1|id), data=data_effects_norm),
  lmer(S ~ original_S + temperature_c + (1|id), data=data_effects_norm),
  lmer(L ~ S + original_S + temperature_c + (1|id), data=data_effects_norm),
  lmer(mean_bs ~ original_S + temperature_c + S + (1|id), data=data_effects_norm),
  lmer(C ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(L.S ~ original_S + S + L + temperature_c + (1|id), data=data_effects_norm),
  lmer(mfcl ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(modularity ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(gensd ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(vulsd ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(log_com_cv ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(log_mean_cv ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(avg_ppmr ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(total_biomass ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(basal ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(intermediate ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(top ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(max_sim ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  # L %~~% S,
  L %~~% mean_bs,
  gensd %~~% mean_bs,
  log_mean_cv %~~% mean_bs,
  avg_ppmr %~~% mean_bs,
  intermediate %~~% mean_bs,
  L.S %~~% mfcl,
  L.S %~~% modularity,
  L.S %~~% gensd,
  L.S %~~% vulsd,
  L.S %~~% log_mean_cv,
  L.S %~~% avg_ppmr,
  L.S %~~% total_biomass,
  L.S %~~% top,
  L.S %~~% basal,
  L.S %~~% intermediate,
  L.S %~~% max_sim,
  C %~~% L.S,
  C %~~% mfcl,
  C %~~% modularity,
  C %~~% gensd,
  C %~~% vulsd,
  C %~~% log_mean_cv,
  C %~~% avg_ppmr,
  C %~~% total_biomass,
  C %~~% top,
  C %~~% basal,
  C %~~% intermediate,
  C %~~% max_sim,
  mfcl %~~% modularity,
  mfcl %~~% gensd,
  mfcl %~~% vulsd,
  mfcl %~~% avg_ppmr,
  mfcl %~~% total_biomass,
  mfcl %~~% basal,
  mfcl %~~% intermediate,
  mfcl %~~% top,
  mfcl %~~% max_sim,
  modularity %~~% vulsd,
  modularity %~~% log_mean_cv,
  modularity %~~% avg_ppmr,
  modularity %~~% total_biomass,
  modularity %~~% basal,
  modularity %~~% max_sim,
  gensd %~~% vulsd,
  gensd %~~% log_com_cv,
  gensd %~~% avg_ppmr,
  gensd %~~% total_biomass,
  gensd %~~% basal,
  gensd %~~% intermediate,
  gensd %~~% top,
  gensd %~~% max_sim,
  vulsd %~~% log_com_cv,
  vulsd %~~% log_mean_cv,
  vulsd %~~% avg_ppmr,
  vulsd %~~% total_biomass,
  vulsd %~~% top,
  vulsd %~~% intermediate,
  vulsd %~~% max_sim,
  log_com_cv %~~% log_mean_cv,
  log_com_cv %~~% basal,
  log_com_cv %~~% intermediate,
  log_com_cv %~~% top,
  log_mean_cv %~~% basal,
  log_mean_cv %~~% top,
  log_mean_cv %~~% intermediate,
  log_mean_cv %~~% max_sim,
  avg_ppmr %~~% total_biomass,
  avg_ppmr %~~% basal,
  avg_ppmr %~~% intermediate,
  avg_ppmr %~~% top,
  avg_ppmr %~~% max_sim,
  total_biomass %~~% basal,
  total_biomass %~~% top,
  basal %~~% intermediate,
  basal %~~% top,
  basal %~~% max_sim,
  max_sim %~~% intermediate,
  top %~~% intermediate,
  max_sim %~~% top,
  data=data_effects_norm
)


sems_effects <- summary(modlist)
all_coeffs <- sems_effects$coefficients
write.csv(all_coeffs[order(all_coeffs$P.Value, -abs(all_coeffs$Estimate), decreasing = F),], file='sems-coeffs-effects.csv')

#### The following code is to create the effect size figures (Figures 2 and S1 in the manuscript)
#### For food web properties:
properties <- c("S", "L", "LxS", "C", "gensd", "vulsd", "max_sim", "basal", "intermediate", "top", "mfcl", "modularity")
properties_names <- c('S', 'L', 'L/S', 'C', 'GenSD', 'VulSD', 'MxSim', 'B', 'I', 'T', 'MFCL', 'Q')

#### These are the properties identified by the SEM above as being most significantly directly affected by temperature
properties <- c("S", "L", "LxS", "C", "vulsd", "basal", "intermediate", "top", "modularity")
properties_names <- c('S', 'L', 'L/S', 'C', 'VulSD', 'B', 'I', 'T', 'Q')

########################## This code generates Figure 2a of the manuscript ##########################
pdf('figure_2a.pdf', height = 12, width=10)
xdist=rev(seq(2,(length(properties))*4, by=4))
par(mar=c(7,9,0.5,0.5), mgp=c(0,2,0))
to_add <- FALSE
temps <- seq(0,40,by=5)
pal <- brewer.pal(length(temps), "RdBu")
pal[5] <- "#B3B3B3"
x_offset <- -1.3
for(t in rev(temps)){
  cur_xdist <- xdist + x_offset
  cur_col <- pal[match(t, rev(temps))]
  data_effects_sum_temp <- subset(data_effects_sum_reformat, (temperature == t & property %in% properties))
  data_effects_sum_temp <- data_effects_sum_temp[match(properties, data_effects_sum_temp$property),]
  pred <- data_effects_sum_temp$mean - 1     ### MEAN PREDICTED
  
  lwqt <- pred-data_effects_sum_temp$se
  hqt <- pred+data_effects_sum_temp$se
  
  matplot(rbind(lwqt, hqt), rbind(cur_xdist, cur_xdist), type="l", col=cur_col, xaxt = "n",lty=1, lwd=4, xlab="", ylab="", ylim=c(0.5, length(properties)*4), main="", axes=F, xlim=c(-.5,0.62), add=to_add)  
  
  if(!to_add){
    axis(1, cex.axis= 3, font=1, family='serif', las=1)
    axis(2, at = xdist, labels = paste(properties_names[which(properties == data_effects_sum_temp$property)]), cex.axis = 3, las=1, font=1, family='serif')
    text("a", x=-0.48, y=32, cex=3.5, font=1, family='serif')
    mtext("Effect size", side=1, line=5, cex=3.5, font=1, family='serif')
    
    abline(v=0)
    
    # abline(h=xdist+1.5)
    
    box(lwd=3)  
  }
  
  points(pred, cur_xdist, pch=15, cex=2, col= cur_col)
  to_add <- TRUE
  
  x_offset <- x_offset + .4
}

legend('topright', legend=(temps), cex=1.8, col = rev(pal[0:length(temps)]), pch=15, title=expression(paste("Temperature ("^"o","C)")), bty = 'n')

dev.off()
######################################################################################################

#### For community properties
properties <- c("total_biomass", "mean_bs", "avg_ppmr", "log_mean_cv", "log_com_cv")
properties_names <- c('Total\nbiomass', 'Avg BS', 'Avg PPMR', 'Population\nstability', 'Community\nstability')

########################## This code generates Figure 2b of the manuscript ##########################
pdf('figure_2b.pdf', height = 10, width=10)
xdist=rev(seq(2,(length(properties))*4, by=4))
par(mar=c(7,13,0.5,0.5), mgp=c(0,2,0))
to_add <- FALSE
temps <- seq(0,40,by=5)
pal <- brewer.pal(length(temps), "RdBu")
pal[5] <- "#B3B3B3"
x_offset <- -1.3
for(t in rev(temps)){
  cur_xdist <- xdist + x_offset
  cur_col <- pal[match(t, rev(temps))]
  data_effects_sum_temp <- subset(data_effects_sum_reformat, (temperature == t & property %in% properties))
  data_effects_sum_temp <- data_effects_sum_temp[match(properties, data_effects_sum_temp$property),]
  pred <- data_effects_sum_temp$mean - 1     ### MEAN PREDICTED
  
  lwqt <- pred-data_effects_sum_temp$se
  hqt <- pred+data_effects_sum_temp$se
  
  matplot(rbind(lwqt, hqt), rbind(cur_xdist, cur_xdist), type="l", col=cur_col, xaxt = "n",lty=1, lwd=4, xlab="", ylab="", ylim=c(0.5, length(properties)*4), main="", axes=F, xlim=c(-1,1.8), add=to_add)  
  
  if(!to_add){
    axis(1, cex.axis= 2.6, font=1, family='serif', las=1)
    axis(2, at = xdist, labels = paste(properties_names[which(properties == data_effects_sum_temp$property)]), cex.axis = 2.6, las=1, font=1, family='serif')
    text("b", x=-.9, y=20, cex=3.5, font=1, family='serif')
    mtext("Effect size", side=1, line=5, cex=3.5, font=1, family='serif')
    abline(v=0)
    box(lwd=3)  
  }
  
  points(pred, cur_xdist, pch=15, cex=2, col= cur_col)
  to_add <- TRUE
  
  x_offset <- x_offset + .4
}

legend('bottomright', legend=(temps), cex=1.8, col = rev(pal[0:length(temps)]), pch=15, title=expression(paste("Temperature ("^"o","C)")), bty = 'n')

dev.off()
######################################################################################################


### After comments by the reviewers from Proceedings B, we decided to change Figure 2 in the
### manuscript to show the log-ratio of the response instead of the ratio. 
### This requires the following implementation in the code

data_effects <- cur_data_invaded
for(p in props){
  eval(parse(text=paste0('data_effects$',p,' <- log(cur_data_invaded[match(cur_data_orig_invaded$unique_id, cur_data_invaded$unique_id),]$', p, '/cur_data_orig_invaded$',p,')')))
}

data_effects[(data_effects == Inf)] <- 0
data_effects[(data_effects == -Inf)] <- 0
data_effects[is.na(data_effects)] <- 0


data_effects_not_invaded <- cur_data_not_invaded
for(p in props){
  eval(parse(text=paste0('data_effects_not_invaded$',p,' <- log(cur_data_not_invaded[match(cur_data_orig_not_invaded$unique_id, cur_data_not_invaded$unique_id),]$', p, '/cur_data_orig_not_invaded$',p,')')))  
}

data_effects_not_invaded[(data_effects_not_invaded == -Inf)] <- 0
data_effects_not_invaded[(data_effects_not_invaded == Inf)] <- 0
data_effects_not_invaded[is.na(data_effects_not_invaded)] <- 0

##### Once we have had calculated the effect of successful and unsuccessful invaders on the recipient communities 
##### (quantified as the ratio of the value of each property before vs. after the invasion event), we can then check 
##### whether the effects of the invasion vs. non-invasion are different.
differences <- NULL
for(p in props){
  for(t in seq(0,40,by=5)){
    group_1 <- eval(parse(text=paste0('subset(data_effects, temperature_c == t)$',p)))
    group_2 <- eval(parse(text=paste0('subset(data_effects_not_invaded, temperature_c == t)$',p)))
    test_res <- wilcox.test(group_1, group_2)
    n1 <- length(group_1)
    n2 <- length(group_2)
    W <- test_res$statistic
    effect_size <- W/(n1*n2)
    
    cur_out <- data.frame(property = p, temperature = t, median_1=round(median(group_1),2), median_2=round(median(group_2),2), U=W, effect_size=round(effect_size,2), n1, n2, p.val=round(test_res$p.value,3))
    if(is.null(differences)){
      differences <- cur_out
    }else{
      differences <- rbind(differences, cur_out)
    }
    
  }
}

differences[which(differences$p.val < 0.05),]

write.csv(differences, file = 'diff-in-effects-log-ratios.csv')



data_effects$original_S <- cur_data_orig_invaded[match(cur_data_invaded$unique_id, cur_data_orig_invaded$unique_id),]$S
data_effects_norm <- data_effects
data_effects_norm$id <- data_effects$id

#### First we investigate the direct and indirect influence of temperature on the effects of invasions on communities using SEMs
#### This SEM produces Figure 3 in the manuscript
modlist <- psem(
  lmer(original_S ~ temperature_c + (1|id), data=data_effects_norm),
  lmer(S ~ original_S + temperature_c + (1|id), data=data_effects_norm),
  lmer(L ~ S + original_S + temperature_c + (1|id), data=data_effects_norm),
  lmer(mean_bs ~ original_S + temperature_c + S + (1|id), data=data_effects_norm),
  lmer(C ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(L.S ~ original_S + S + L + temperature_c + (1|id), data=data_effects_norm),
  lmer(mfcl ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(modularity ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(gensd ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(vulsd ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(log_com_cv ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(log_mean_cv ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(avg_ppmr ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(total_biomass ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(basal ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(intermediate ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(top ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  lmer(max_sim ~ original_S + temperature_c + S + L + (1|id), data=data_effects_norm),
  # L %~~% S,
  L %~~% mean_bs,
  gensd %~~% mean_bs,
  log_mean_cv %~~% mean_bs,
  avg_ppmr %~~% mean_bs,
  intermediate %~~% mean_bs,
  L.S %~~% mfcl,
  L.S %~~% modularity,
  L.S %~~% gensd,
  L.S %~~% vulsd,
  L.S %~~% log_mean_cv,
  L.S %~~% avg_ppmr,
  L.S %~~% total_biomass,
  L.S %~~% top,
  L.S %~~% basal,
  L.S %~~% intermediate,
  L.S %~~% max_sim,
  C %~~% L.S,
  C %~~% mfcl,
  C %~~% modularity,
  C %~~% gensd,
  C %~~% vulsd,
  C %~~% log_mean_cv,
  C %~~% avg_ppmr,
  C %~~% total_biomass,
  C %~~% top,
  C %~~% basal,
  C %~~% intermediate,
  C %~~% max_sim,
  mfcl %~~% modularity,
  mfcl %~~% gensd,
  mfcl %~~% vulsd,
  mfcl %~~% avg_ppmr,
  mfcl %~~% total_biomass,
  mfcl %~~% basal,
  mfcl %~~% intermediate,
  mfcl %~~% top,
  mfcl %~~% max_sim,
  modularity %~~% vulsd,
  modularity %~~% log_mean_cv,
  modularity %~~% avg_ppmr,
  modularity %~~% total_biomass,
  modularity %~~% basal,
  modularity %~~% max_sim,
  gensd %~~% vulsd,
  gensd %~~% log_com_cv,
  gensd %~~% avg_ppmr,
  gensd %~~% total_biomass,
  gensd %~~% basal,
  gensd %~~% intermediate,
  gensd %~~% top,
  gensd %~~% max_sim,
  vulsd %~~% log_com_cv,
  vulsd %~~% log_mean_cv,
  vulsd %~~% avg_ppmr,
  vulsd %~~% total_biomass,
  vulsd %~~% top,
  vulsd %~~% intermediate,
  vulsd %~~% max_sim,
  log_com_cv %~~% log_mean_cv,
  log_com_cv %~~% basal,
  log_com_cv %~~% intermediate,
  log_com_cv %~~% top,
  log_mean_cv %~~% basal,
  log_mean_cv %~~% top,
  log_mean_cv %~~% intermediate,
  log_mean_cv %~~% max_sim,
  avg_ppmr %~~% total_biomass,
  avg_ppmr %~~% basal,
  avg_ppmr %~~% intermediate,
  avg_ppmr %~~% top,
  avg_ppmr %~~% max_sim,
  total_biomass %~~% basal,
  total_biomass %~~% top,
  basal %~~% intermediate,
  basal %~~% top,
  basal %~~% max_sim,
  max_sim %~~% intermediate,
  top %~~% intermediate,
  max_sim %~~% top,
  data=data_effects_norm
)


sems_effects <- summary(modlist)
all_coeffs <- sems_effects$coefficients
write.csv(all_coeffs[order(all_coeffs$P.Value, -abs(all_coeffs$Estimate), decreasing = F),], file='sems-coeffs-effects.csv')


## Here we calculate the mean and standard error / confident intervals of the effects
## and perform some transformations on the data format to facilitate plotting
data_effects_sum <- getPropertiesSummary(data_effects)
names(data_effects_sum)[12:15] <- c("LxS.mean", "LxS.sd", "LxS.se", "LxS.CI")

data_effects_sum_reformat <- melt(data_effects_sum, id.vars = 'temperature_c', measure.vars = names(data_effects_sum)[grep('\\.mean', names(data_effects_sum))])
data_effects_sum_reformat$variable <- unlist(lapply(as.character(data_effects_sum_reformat$variable), function(x) {strsplit(x, '\\.')[[1]][1]}))

names(data_effects_sum_reformat) <- c('temperature', 'property', 'mean')

data_effects_sum_reformat$CI <- melt(data_effects_sum, id.vars = 'temperature_c', measure.vars = names(data_effects_sum)[grep('.CI', names(data_effects_sum))])$value
data_effects_sum_reformat$se <- melt(data_effects_sum, id.vars = 'temperature_c', measure.vars = names(data_effects_sum)[grep('\\.se', names(data_effects_sum))])$value

#### The following code is to create the effect size figures (Figures 2 and S1 in the manuscript)
#### For food web properties:
properties <- c("S", "L", "LxS", "C", "gensd", "vulsd", "max_sim", "basal", "intermediate", "top", "mfcl", "modularity")
properties_names <- c('S', 'L', 'L/S', 'C', 'GenSD', 'VulSD', 'MxSim', 'B', 'I', 'T', 'MFCL', 'Q')

#### These are the properties identified by the SEM above as being most significantly directly affected by temperature
properties <- c("S", "L", "LxS", "C", "vulsd", "intermediate", "top", "modularity")
properties_names <- c('S', 'L', 'L/S', 'C', 'VulSD', 'I', 'T', 'Q')

########################## This code generates new version of Figure 2a of the manuscript after changing to log-ratio responses ##########################
pdf('figure_2a-new.pdf', height = 12, width=10)
xdist=rev(seq(2,(length(properties))*4, by=4))
par(mar=c(7,9,0.5,0.5), mgp=c(0,2,0))
to_add <- FALSE
temps <- seq(0,40,by=5)
pal <- brewer.pal(length(temps), "RdBu")
pal[5] <- "#B3B3B3"
x_offset <- -1.3
for(t in rev(temps)){
  cur_xdist <- xdist + x_offset
  cur_col <- pal[match(t, rev(temps))]
  data_effects_sum_temp <- subset(data_effects_sum_reformat, (temperature == t & property %in% properties))
  data_effects_sum_temp <- data_effects_sum_temp[match(properties, data_effects_sum_temp$property),]
  pred <- data_effects_sum_temp$mean     ### MEAN PREDICTED
  
  lwqt <- pred-data_effects_sum_temp$se
  hqt <- pred+data_effects_sum_temp$se
  
  matplot(rbind(lwqt, hqt), rbind(cur_xdist, cur_xdist), type="l", col=cur_col, xaxt = "n",lty=1, lwd=4, xlab="", ylab="", ylim=c(0.5, length(properties)*4), main="", axes=F, xlim=c(-.77,0.5), add=to_add)  
  
  if(!to_add){
    axis(1, cex.axis= 2.7, font=1, family='serif', las=1)
    axis(2, at = xdist, labels = paste(properties_names[which(properties == data_effects_sum_temp$property)]), cex.axis = 3, las=1, font=1, family='serif')
    text("a", x=-0.75, y=32, cex=3.5, font=1, family='serif')
    mtext("Effect size", side=1, line=5, cex=3.5, font=1, family='serif')
    
    abline(v=0)
    
    # abline(h=xdist+1.5)
    
    box(lwd=3)  
  }
  
  points(pred, cur_xdist, pch=15, cex=2, col= cur_col)
  to_add <- TRUE
  
  x_offset <- x_offset + .4
}

legend('topright', legend=(temps), cex=1.8, col = rev(pal[0:length(temps)]), pch=15, title=expression(paste("Temperature ("^"o","C)")), bty = 'n')

dev.off()



#### For community properties
properties <- c("total_biomass", "mean_bs", "avg_ppmr", "log_mean_cv", "log_com_cv")
properties_names <- c('Total\nbiomass', 'Avg BS', 'Avg PPMR', 'Population\nstability', 'Community\nstability')

########################## This code generates Figure 2b of the manuscript ##########################
pdf('figure_2b-new.pdf', height = 10, width=10)
xdist=rev(seq(2,(length(properties))*4, by=4))
par(mar=c(7,13,0.5,0.5), mgp=c(0,2,0))
to_add <- FALSE
temps <- seq(0,40,by=5)
pal <- brewer.pal(length(temps), "RdBu")
pal[5] <- "#B3B3B3"
x_offset <- -1.3
for(t in rev(temps)){
  cur_xdist <- xdist + x_offset
  cur_col <- pal[match(t, rev(temps))]
  data_effects_sum_temp <- subset(data_effects_sum_reformat, (temperature == t & property %in% properties))
  data_effects_sum_temp <- data_effects_sum_temp[match(properties, data_effects_sum_temp$property),]
  pred <- data_effects_sum_temp$mean    ### MEAN PREDICTED
  
  lwqt <- pred-data_effects_sum_temp$se
  hqt <- pred+data_effects_sum_temp$se
  
  matplot(rbind(lwqt, hqt), rbind(cur_xdist, cur_xdist), type="l", col=cur_col, xaxt = "n",lty=1, lwd=4, xlab="", ylab="", ylim=c(0.5, length(properties)*4), main="", axes=F, xlim=c(-1.6,0.1), add=to_add)  
  
  if(!to_add){
    axis(1, cex.axis= 2.6, font=1, family='serif', las=1)
    axis(2, at = xdist, labels = paste(properties_names[which(properties == data_effects_sum_temp$property)]), cex.axis = 2.6, las=1, font=1, family='serif')
    text("b", x=-1.5, y=20, cex=3.5, font=1, family='serif')
    mtext("Effect size", side=1, line=5, cex=3.5, font=1, family='serif')
    abline(v=0)
    box(lwd=3)  
  }
  
  points(pred, cur_xdist, pch=15, cex=2, col= cur_col)
  to_add <- TRUE
  
  x_offset <- x_offset + .4
}

# legend('bottomright', legend=(temps), cex=1.8, col = rev(pal[0:length(temps)]), pch=15, title=expression(paste("Temperature ("^"o","C)")), bty = 'n')

dev.off()

######################################################################################################

### for non-invaded communities

## Here we calculate the mean and standard error / confident intervals of the effects
## and perform some transformations on the data format to facilitate plotting
data_effects_sum <- getPropertiesSummary(data_effects_not_invaded)
names(data_effects_sum)[12:15] <- c("LxS.mean", "LxS.sd", "LxS.se", "LxS.CI")

data_effects_sum_reformat <- melt(data_effects_sum, id.vars = 'temperature_c', measure.vars = names(data_effects_sum)[grep('\\.mean', names(data_effects_sum))])
data_effects_sum_reformat$variable <- unlist(lapply(as.character(data_effects_sum_reformat$variable), function(x) {strsplit(x, '\\.')[[1]][1]}))

names(data_effects_sum_reformat) <- c('temperature', 'property', 'mean')

data_effects_sum_reformat$CI <- melt(data_effects_sum, id.vars = 'temperature_c', measure.vars = names(data_effects_sum)[grep('.CI', names(data_effects_sum))])$value
data_effects_sum_reformat$se <- melt(data_effects_sum, id.vars = 'temperature_c', measure.vars = names(data_effects_sum)[grep('\\.se', names(data_effects_sum))])$value

properties <- c("S", "L", "LxS", "gensd", "mfcl", "basal", "intermediate")
properties_names <- c('S', 'L', 'L/S', 'GenSD', 'MFCL', 'B', 'I')

########################## This code generates Figure S4 of the manuscript ##########################
pdf('figure_s4.pdf', height = 14, width=10)
xdist=rev(seq(2,(length(properties))*4, by=4))
par(mar=c(7,9,0.5,0.5), mgp=c(0,2,0))
to_add <- FALSE
temps <- seq(0,40,by=5)
pal <- brewer.pal(length(temps), "RdBu")
pal[5] <- "#B3B3B3"
x_offset <- -1.3
for(t in rev(temps)){
  cur_xdist <- xdist + x_offset
  cur_col <- pal[match(t, rev(temps))]
  data_effects_sum_temp <- subset(data_effects_sum_reformat, (temperature == t & property %in% properties))
  data_effects_sum_temp <- data_effects_sum_temp[match(properties, data_effects_sum_temp$property),]
  pred <- data_effects_sum_temp$mean - 1     ### MEAN PREDICTED
  
  lwqt <- pred-data_effects_sum_temp$se
  hqt <- pred+data_effects_sum_temp$se
  
  matplot(rbind(lwqt, hqt), rbind(cur_xdist, cur_xdist), type="l", col=cur_col, xaxt = "n",lty=1, lwd=4, xlab="", ylab="", ylim=c(0.5, length(properties)*4), main="", axes=F, xlim=c(-.65,.3), add=to_add)  
  
  if(!to_add){
    axis(1, cex.axis= 2.6, font=1, family='serif', las=1)
    axis(2, at = xdist, labels = paste(properties_names[which(properties == data_effects_sum_temp$property)]), cex.axis = 2.6, las=1, font=1, family='serif')
    mtext("Effect size", side=1, line=5, cex=3.5, font=1, family='serif')
    abline(v=0)
    box(lwd=3)  
  }
  
  points(pred, cur_xdist, pch=15, cex=2, col= cur_col)
  to_add <- TRUE
  
  x_offset <- x_offset + .4
}

legend('topright', legend=(temps), cex=1.8, col = rev(pal[0:length(temps)]), pch=15, title=expression(paste("Temperature ("^"o","C)")), bty = 'n')

dev.off()
######################################################################################################



##### This is to generate Figure S4 using log ratios instead of ratios (as per suggestion from reviewers)
## Here we calculate the mean and standard error / confident intervals of the effects
## and perform some transformations on the data format to facilitate plotting
data_effects_sum <- getPropertiesSummary(data_effects_not_invaded)
names(data_effects_sum)[12:15] <- c("LxS.mean", "LxS.sd", "LxS.se", "LxS.CI")

data_effects_sum_reformat <- melt(data_effects_sum, id.vars = 'temperature_c', measure.vars = names(data_effects_sum)[grep('\\.mean', names(data_effects_sum))])
data_effects_sum_reformat$variable <- unlist(lapply(as.character(data_effects_sum_reformat$variable), function(x) {strsplit(x, '\\.')[[1]][1]}))

names(data_effects_sum_reformat) <- c('temperature', 'property', 'mean')

data_effects_sum_reformat$CI <- melt(data_effects_sum, id.vars = 'temperature_c', measure.vars = names(data_effects_sum)[grep('.CI', names(data_effects_sum))])$value
data_effects_sum_reformat$se <- melt(data_effects_sum, id.vars = 'temperature_c', measure.vars = names(data_effects_sum)[grep('\\.se', names(data_effects_sum))])$value

properties <- c("S", "L", "LxS", "gensd", "mfcl", "basal")
properties_names <- c('S', 'L', 'L/S', 'GenSD', 'MFCL', 'B')

########################## This code generates Figure S4 of the manuscript ##########################
pdf('figure_s4-new.pdf', height = 14, width=10)
xdist=rev(seq(2,(length(properties))*4, by=4))
par(mar=c(7,9,0.5,0.5), mgp=c(0,2,0))
to_add <- FALSE
temps <- seq(0,40,by=5)
pal <- brewer.pal(length(temps), "RdBu")
pal[5] <- "#B3B3B3"
x_offset <- -1.3
for(t in rev(temps)){
  cur_xdist <- xdist + x_offset
  cur_col <- pal[match(t, rev(temps))]
  data_effects_sum_temp <- subset(data_effects_sum_reformat, (temperature == t & property %in% properties))
  data_effects_sum_temp <- data_effects_sum_temp[match(properties, data_effects_sum_temp$property),]
  pred <- data_effects_sum_temp$mean     ### MEAN PREDICTED
  
  lwqt <- pred-data_effects_sum_temp$se
  hqt <- pred+data_effects_sum_temp$se
  
  matplot(rbind(lwqt, hqt), rbind(cur_xdist, cur_xdist), type="l", col=cur_col, xaxt = "n",lty=1, lwd=4, xlab="", ylab="", ylim=c(0.5, length(properties)*4), main="", axes=F, xlim=c(-1.2,.3), add=to_add)  
  
  if(!to_add){
    axis(1, cex.axis= 2.6, font=1, family='serif', las=1)
    axis(2, at = xdist, labels = paste(properties_names[which(properties == data_effects_sum_temp$property)]), cex.axis = 2.6, las=1, font=1, family='serif')
    mtext("Effect size", side=1, line=5, cex=3.5, font=1, family='serif')
    abline(v=0)
    box(lwd=3)  
  }
  
  points(pred, cur_xdist, pch=15, cex=2, col= cur_col)
  to_add <- TRUE
  
  x_offset <- x_offset + .4
}

legend('topright', legend=(temps), cex=1.8, col = rev(pal[0:length(temps)]), pch=15, title=expression(paste("Temperature ("^"o","C)")), bty = 'n')

dev.off()
######################################################################################################



