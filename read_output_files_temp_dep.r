

source('utils.R')

setwd('./output')

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
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
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



ObtainRow <- function(net, out_row){
  require(igraph)
  if(out_row[2] == 3){
    #print('I am in')
 
    #here we include the information about the type of experiment and the value
    if(net$experiment_type == 'neutral'){
      out_row <- append(out_row, c('neutral', 'NA'));
    }else if(net$experiment_type == 'efficiency'){
      out_row <- append(out_row, c('efficiency', net$eff_type));
    }else if(net$experiment_type == 'body_size'){
      out_row <- append(out_row, 'body_size');
      
      if(net$inv_type == 'better'){
        out_row <- append(out_row, net$bs_fraction);
      }else{
        out_row <- append(out_row, -(net$bs_fraction));
      }
      
    }else if(net$experiment_type == 'biomass'){
      out_row <- append(out_row, c('biomass', net$initial_biomass));
    }else if(net$experiment_type == 'hill_exponent'){
      out_row <- append(out_row, c('hill_exponent', net$hill_exp));
    }else if(net$experiment_type == 'generalism'){
      out_row <- append(out_row, 'generalism');
      
      if(net$generalism_type == 'add'){
        out_row <- append(out_row, net$generalism_frac);
      }else{
        out_row <- append(out_row, -(net$generalism_frac));
      }
    }else if(net$experiment_type == 'vulnerability'){
      out_row <- append(out_row, c('vulnerability', -(net$vulnerability_frac)));
    }
    
  }
  
  
  out_row <- append(out_row, as.integer(net$S));
  out_row <- append(out_row, net$C);
  out_row <- append(out_row, net$original_C);
  out_row <- append(out_row, mean(net$BS));
  
  #mean_ppmr=numeric(len), mfcl=numeric(len), 
  
  if(is.null(net$stats)){
    net <- CalculateFoodWebStats(net);
  }
  
  out_row <- append(out_row, net$stats$sd_vul);
  out_row <- append(out_row, net$stats$sd_gen);
  
  # if(net$C < 0.2 & net$S < 60){
  
  #print(out_row)
  fcstats <- FoodChainStats(net$M);
    
  out_row <- append(out_row, fcstats$overall_mean_length);
  out_row <- append(out_row, fcstats$path_length_sd);
    
  out_row <- append(out_row, (fcstats$overall_number_of_paths));
  # }else{
  #   out_row <- append(out_row, c(0.0,0.0,0.0));
  # }
  
  out_row <- append(out_row, AveragePPMR(net));
  
  out_row <- append(out_row, net$stats$B);
  out_row <- append(out_row, net$stats$I); 
  out_row <- append(out_row, net$stats$T);
  out_row <- append(out_row, net$stats$Ca);
  
  # out_row <- append(out_row, net$stats$M);
  
  g <-  graph.adjacency(net$M)
  
  out_row <- append(out_row, components(g)$no);
  
  mdlrt <- tryCatch({
    cluster_louvain(as.undirected(g))
  }, error = function(err){
    NA
  })
  
  if(is.na(mdlrt)){
    out_row <- append(out_row, c(NA, NA))
  }else{
    out_row <- append(out_row, c(modularity(mdlrt), length(mdlrt)))
  }
  
  
  
  out_row <- append(out_row, net$stats$MxSim);
  
  
  if(!is.null(net$invasive)){
    extinct <- FALSE;
    if( !(net$n_inv %in% net$niche) ){
      extinct <- TRUE;
    }
    out_row <- append(out_row, extinct);
    
    out_row <- append(out_row, net$orig_fw$invasive$BS);
    out_row <- append(out_row, net$dynamics[dim(net$dynamics)[1],-1][net$orig_fw$invasive$index]);
    out_row <- append(out_row, net$orig_fw$TP[net$orig_fw$invasive$index]);
    
    out_row <- append(out_row, net$orig_fw$stats$in_deg[net$orig_fw$invasive$index]);
    out_row <- append(out_row, net$orig_fw$stats$out_deg[net$orig_fw$invasive$index]);
    
    out_row <- append(out_row, net$orig_fw$stats$in_deg[net$orig_fw$invasive$index]/net$orig_fw$stats$LxS);
    out_row <- append(out_row, net$orig_fw$stats$out_deg[net$orig_fw$invasive$index]/net$orig_fw$stats$LxS);
    
    
    
  }else{    
    out_row <- append(out_row, c('NA', -1, -1, -1, -1, -1, -1, -1));
  }
  
  
  
  #print(net$dynamics[dim(net$dynamics)[1],-1])
  
  if(is.null(net$dynamics)){
    #print('dynamics is null')
    #print(net$dynamics[dim(net$dynamics)[1],-1])
    out_row <- append(out_row, c('NA','NA','NA','NA','NA'));
  }else{
    #print(net$dynamics[dim(net$dynamics)[1],-1])
    out_row <- append(out_row, sum(net$Final));
    
    cur_dynamics <- net$dynamics[,which(net$dynamics[21901,] != 0)]
    out_row <- append(out_row, sd(rowSums(cur_dynamics[,-1]))/mean(rowSums(cur_dynamics[,-1])));
    out_row <- append(out_row, mean(apply(cur_dynamics[,-1], 2, sd)/apply(cur_dynamics[,-1], 2, mean)));
    
    ##### variability for the last 100 years
    out_row <- append(out_row, sd(rowSums(cur_dynamics[18250:21900,-1]))/mean(rowSums(cur_dynamics[18250:21900,-1])));
    out_row <- append(out_row, mean(apply(cur_dynamics[18250:21900,-1], 2, sd)/apply(cur_dynamics[18250:21900,-1], 2, mean)));
  }
  
 
  #print(out_row)
  
  return(out_row);
  
}



filenames <- list.files(pattern='RFile-FW-');

networks <- data.frame('id'='', 'type'='',
                       'temperature'='', 'eutrophication'='',
                       'exp_type'='', 'exp_val'='', 
                       'S'='', 'C'='', 
                       'original_C'='', 
                       'mean_bs'='', 'vulsd'='', 
                       'gensd'='', 'mfcl'='',
                       'ChnSD'='', 'ChnNo'='',
                       'avg_ppmr'='', 'basal'='', 
                       'intermediate'='', 'top'='', 'ca'='', 
                       'components'='', 'modularity'='', 'n_modules'='', 
                       'max_sim'='', 'inv_extinct'='', 
                       'inv_bs' = '', 'inv_biom'='','inv_tp'='',
                       'inv_gen'='', 'inv_vul'='',
                       'inv_gen_norm'='', 'inv_vul_norm'='',
                       'total_biomass'='', 'com_cv'='', 
                       'mean_cv'='', 'com_cv_last_100'='', 
                       'mean_cv_last_100'='', 'sp_lost'='', stringsAsFactors=FALSE);



net_ids <- c();

to_include <- names(which(table(as.numeric(unlist(lapply(filenames, function(x) { strsplit(x, '-')[[1]][8]})))) == 984))

for(i in 1:length(filenames)){
  filename <- filenames[i];
  loaded <- TRUE;
  loaded <- tryCatch({
    load(filename);
  }, error = function(err){
    print(paste('Error loading file', filename));
    return(FALSE);
  })
  
  if(loaded == FALSE) next;
  
  
  filename_split <- strsplit(filename, '-')[[1]]
  
  if(! (filename_split[8] %in% to_include )){
    next
  }
  
  print(filename_split[8])
  
  tempr <- as.numeric(filename_split[which(filename_split == 'temp') + 1])
  eut <- strsplit(filename_split[which(filename_split == 'eutroph') + 1], '\\.')
  eut <- as.numeric(eut[[1]][1])
  if(grepl('before-inv', filename)){
    net <- FW;
    out_row <- c(as.integer(net$fwno), as.integer(1));
    
    ##### the temperature and eutrophication of the treatments
    out_row <- append(out_row, tempr)
    out_row <- append(out_row, eut)
    
    out_row <- append(out_row, c('NA', 'NA'));
    out_row <- ObtainRow(net, out_row);
    out_row <- append(out_row, as.integer(net$S_lost_b4inv$number));
    
    networks <- rbind(networks, out_row);
        
    net <- FW$orig_fw;
    out_row <- c(as.integer(FW$fwno), as.integer(0))
    
    out_row <- append(out_row, tempr)
    out_row <- append(out_row, eut)
    
    out_row <- append(out_row, c('NA', 'NA'));
    out_row <- ObtainRow(net, out_row);
    out_row <- append(out_row, as.integer(0));
    
    networks <- rbind(networks, out_row);
    
  }else{
    net <- invaded;
    current_id <- net$fwno;
    out_row <- c(current_id, as.integer(3));
    
    out_row <- append(out_row, tempr)
    out_row <- append(out_row, eut)
    
    out_row <- ObtainRow(net, out_row);
    out_row <- append(out_row, as.integer(net$S_lost_aftinv$number));
    
    networks <- rbind(networks, out_row);
    
    
    
    # if(!(current_id %in% net_ids) ){
    # net <- invaded$orig_fw;
    # out_row <- c(as.integer(net$fwno), as.integer(2));
    # 
    # out_row <- append(out_row, tempr)
    # out_row <- append(out_row, eut)
    # 
    # out_row <- append(out_row, c('NA', 'NA'));
    # out_row <- ObtainRow(net, out_row);
    # out_row <- append(out_row, as.integer(0));
    # networks <- rbind(networks, out_row);
    #   net_ids <- append(net_ids, current_id);
    # }

   
    
  }
  
}


networks <- networks[-1,]
networks$id <- as.numeric(networks$id)
networks$type <- as.numeric(networks$type)
networks$temperature <- as.numeric(networks$temperature)
networks$eutrophication <- as.numeric(networks$eutrophication)

networks$S <- as.numeric(networks$S)
networks$C <- as.numeric(networks$C)

networks$exp_val <- as.numeric(networks$exp_val)

networks <- networks[with(networks, order(id, temperature, eutrophication, type)),]

  
# write.csv(networks, file='output_new.csv');


require(ggplot2)

cur_data <- subset(networks, ( (exp_type=='body_size'))) # & (id %in% selected_ids)))

cur_data <- cur_data[with(cur_data, order(id, temperature, eutrophication, exp_val)),]

cur_data$mean_cv <- as.numeric(cur_data$mean_cv)
cur_data$inv_extinct <- as.logical(cur_data$inv_extinct)
cur_data$inv_extinct <- as.numeric(cur_data$inv_extinct)

cur_data$sp_lost <- as.numeric(cur_data$sp_lost)


ggplot(cur_data, aes(temperature, S, color=as.factor(eutrophication))) +
  facet_wrap(~exp_val) + geom_point()


cur_data_sum <- summarySE(cur_data, measurevar="inv_extinct", groupvars=c("temperature","eutrophication","exp_val"))

pd <- position_dodge(1)
ggplot(cur_data_sum, aes(temperature, (1-inv_extinct), colour=as.factor(eutrophication))) + 
  geom_errorbar(aes(ymin=(1-inv_extinct)-se, ymax=(1-inv_extinct)+se), width=.1, position=pd) +
  facet_wrap(~exp_val) +
  geom_line(position=pd) +
  geom_point(position=pd, size=2) +
  theme_bw()
  
ggplot(cur_data_sum, aes(temperature, 1-inv_extinct, color=as.factor(eutrophication))) +
  facet_wrap(~exp_val) + geom_point()

ggplot(cur_data_sum, aes(exp_val, (1-inv_extinct), color=as.numeric(temperature))) +
  facet_wrap(~eutrophication) + geom_point() + geom_path()


temp <- subset(cur_data_sum, (temperature == 273.15 & eutrophication == 10))

ggplot(temp, aes(exp_val, inv_extinct)) + geom_point()

cur_data_sum <- summarySE(cur_data, measurevar="inv_extinct", groupvars=c("temperature","eutrophication","exp_val"))





############### these are the original networks after transient dynamics

selected_ids <- as.vector(names(which(table(networks$id) == 1107)))
cur_data_orig <- subset(networks, ( (id %in% selected_ids) & (type == 1)))

cur_data_orig <- cur_data_orig[with(cur_data_orig, order(id, temperature, eutrophication)),]

cur_data_orig$mean_cv <- as.numeric(cur_data_orig$mean_cv)
cur_data_orig$com_cv <- as.numeric(cur_data_orig$com_cv)

cur_data_orig$inv_extinct <- as.logical(cur_data_orig$inv_extinct)
cur_data_orig$inv_extinct <- as.numeric(cur_data_orig$inv_extinct)

cur_data_orig$mean_bs <- as.numeric(cur_data_orig$mean_bs)
cur_data_orig$C <- as.numeric(cur_data_orig$C)
cur_data_orig$modularity <- as.numeric(cur_data_orig$modularity)
cur_data_orig$mfcl <- as.numeric(cur_data_orig$mfcl)
cur_data_orig$intermediate <- as.numeric(cur_data_orig$intermediate)

cur_data_orig$S <- as.numeric(cur_data_orig$S)
cur_data_orig$total_biomass <- as.numeric(cur_data_orig$total_biomass)

ggplot(cur_data_orig, aes(temperature, S, color=as.factor(eutrophication))) + 
  geom_point() + theme_bw()


cur_data_sum_orig <- summarySE(cur_data_orig, measurevar="S", groupvars=c("temperature","eutrophication"))

pd <- position_dodge(1)
ggplot(cur_data_sum_orig, aes(temperature, S, colour=as.factor(eutrophication))) + 
  geom_errorbar(aes(ymin=S-se, ymax=S+se), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=2) +
  theme_bw()

ggplot(cur_data_sum_orig, aes(temperature, S, color=as.factor(eutrophication))) + geom_point(size=3) + theme_bw() + geom_path()


merged_data <- merge(cur_data_sum, cur_data_sum_orig, by=c('temperature', 'eutrophication'))

ggplot(merged_data, aes(S, 1-inv_extinct, color=(temperature), shape=as.factor(eutrophication))) +
  facet_wrap(~exp_val) + geom_point(size=3) + theme_bw() + ylab('invasion success')



temp <- ddply(cur_data, .(temperature, eutrophication), summarise, N=length(temperature))

ggplot(temp, aes(temperature, N, colour=eutrophication)) + geom_point()

###### for the species lost

cur_data_sum_sp_lost <- summarySE(cur_data, measurevar='sp_lost', groupvars=c("temperature","eutrophication","exp_val"))


pd <- position_dodge(1)
ggplot(cur_data_sum_sp_lost, aes(temperature, sp_lost, colour=as.factor(eutrophication))) + 
  facet_wrap(~exp_val) +
  geom_errorbar(aes(ymin=sp_lost-se, ymax=sp_lost+se), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=2) +
  theme_bw()

ggplot(cur_data_sum_sp_lost, aes(temperature, fraction_lost, color=as.factor(eutrophication))) + 
  geom_point(size=3) + theme_bw() +facet_wrap(~exp_val) 




##################### NOW WITH PROPAGULE PRESSURE ###########################

cur_data <- subset(networks, exp_type=='biomass')

cur_data <- cur_data[with(cur_data, order(id, temperature, eutrophication, exp_val)),]

cur_data$mean_cv <- as.numeric(cur_data$mean_cv)
cur_data$inv_extinct <- as.logical(cur_data$inv_extinct)
cur_data$inv_extinct <- as.numeric(cur_data$inv_extinct)

cur_data$sp_lost <- as.numeric(cur_data$sp_lost)


ggplot(cur_data, aes(temperature, sp_lost, color=as.factor(eutrophication))) +
  facet_wrap(~exp_val) + geom_point()


cur_data_sum <- summarySE(cur_data, measurevar="inv_extinct", groupvars=c("temperature","eutrophication","exp_val"))

pd <- position_dodge(1)
ggplot(cur_data_sum, aes(log10(exp_val), (1-inv_extinct), colour=as.factor(temperature), shape=as.factor(eutrophication))) + 
  geom_errorbar(aes(ymin=inv_extinct-se, ymax=inv_extinct+se), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=2) +
  theme_bw()

ggplot(cur_data_sum, aes(temperature, (1-inv_extinct), color=as.factor(eutrophication))) +
  facet_wrap(~exp_val) + geom_point() + theme_bw() + ylab('invasion success')


temp <- subset(cur_data_sum, (temperature == 303.15 & eutrophication == 20))

ggplot(temp, aes(log10(exp_val), inv_extinct)) + geom_point()

cur_data_sum <- summarySE(cur_data, measurevar="inv_extinct", groupvars=c("temperature","eutrophication","exp_val"))



merged_data <- merge(cur_data_sum, cur_data_sum_orig, by=c('temperature', 'eutrophication'))

ggplot(merged_data, aes(S, 1-inv_extinct, color=(temperature), shape=as.factor(eutrophication))) +
  facet_wrap(~exp_val) + geom_point(size=3) + theme_bw() + ylab('invasion success')


#################### AND THE NEUTRAL EXPECTATION ###########################

cur_data <- subset(networks, (exp_type=='neutral')) # & inv_extinct == FALSE))

cur_data <- cur_data[with(cur_data, order(id, temperature, eutrophication, exp_val)),]

cur_data$mean_cv <- as.numeric(cur_data$mean_cv)
cur_data$inv_extinct <- as.logical(cur_data$inv_extinct)
cur_data$inv_extinct <- as.numeric(cur_data$inv_extinct)

cur_data$sp_lost <- as.numeric(cur_data$sp_lost)


ggplot(cur_data, aes(temperature, S, color=as.factor(eutrophication))) +
  geom_point()


cur_data$original_S <- cur_data$S + cur_data$sp_lost
cur_data$fraction_lost <- cur_data$sp_lost/cur_data$original_S

cur_data_sum <- summarySE(cur_data, measurevar="inv_extinct", groupvars=c("original_S","eutrophication"))

pd <- position_dodge(1)
ggplot(cur_data_sum, aes(temperature, S, colour=as.factor(eutrophication))) + 
  geom_errorbar(aes(ymin=(S)-se, ymax=(S)+se), width=.1, position=pd) +
  # geom_line(position=pd) +
  geom_point(position=pd, size=2) +
  theme_bw() + ylab('S')

ggplot(cur_data_sum, aes(original_S, (1-inv_extinct), color=as.factor(eutrophication))) +
  geom_point(size=3) + theme_bw() + geom_line() + ylab('invasion success')





glm1 <- glm(inv_extinct ~ temperature*original_S*eutrophication, data=cur_data, family='binomial')
glm2 <- glm(inv_extinct ~ temperature*S*eutrophication-temperature:S:eutrophication, data=cur_data, family='binomial')

require(car)
Anova(glm1, glm2)

summary(lm(inv_extinct ~ temperature+eutrophication, data=cur_data_sum))

merged_data_neutral <- merge(cur_data_sum, cur_data_sum_orig, by=c('temperature', 'eutrophication'))

ggplot(merged_data_neutral, aes(temperature ,(S.x-1)/S.y)) +
  geom_point(size=3) + theme_bw() + facet_wrap(~eutrophication)



