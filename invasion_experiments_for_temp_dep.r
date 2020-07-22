## ---------------------------
##
## Script name: invasion_experiments_for_temp_dep.r
##
## Purpose of script: This script implements the invasion experiments on complex food
## webs. Several different invasion types are implemented by different functions, according
## to different traits of the invader.
##
## Towards the end of the script, the experimentation protocol to run the invasion simulations
## across temperature regimes is implemented.
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
## Date Created: October-2012
##
## Copyright (c) Miguel Lurgi, 2018-2020
## Email: miguel.lurgi@swansea.ac.uk
##
## ---------------------------
##
## Notes:
## This script was first developed for the paper:
##
## Lurgi, Galiana, Lopez, Joppa, Montoya (2014) Network complexity and species traits mediate the effects of 
## biological invasions on dynamic food webs. Frontiers in Ecology and Evolution, 2:36, 1-11.
##
## Additions incoporating the temperature dependence effect of the food web on the invasions
## were developed for our current paper:
##
## Sentis, Montoya & Lurgi (2020) Warming indirectly incrases invasion success in food webs. Uploaded to BioRXiv. https://doi.org/10.1101/2020.07.20.211516
##
## for which this script is now presented as supplementary material
## ---------------------------

## sourcing scripts
source('temp-dep-food-web.r');
source('utils.r');

InvadeNetwork <- function(FW, C, fromNiche=TRUE, init_b=0.2){
  x <- 0;
  FW_inv <- FW;
  S_temp <- FW_inv$S + 1;
  y <- S_temp;
  dup_sp <- FALSE;
  
  inv_values <- FALSE;
  if(!is.null(FW_inv$n_inv) & !is.null(FW_inv$r_inv) & !is.null(FW_inv$c_inv)){
    inv_values <- TRUE; 
    print('InvadeNetwork :: I have got the values for the invader');
  }
  
  while(x < S_temp | y > 1 | dup_sp){
    dup_sp <- FALSE;
    
    if(inv_values){
      n <- FW_inv$n_inv;
    }else{
      n <- runif(1,0,1);      
      while(n %in% FW$orig_fw$niche)
        n <- runif(1,0,1);
      
      FW_inv$n_inv <- n;
    }
    
    n_temp <- FW_inv$niche;
    
    n_temp[length(n_temp)+1] <- n;
    n_temp <- sort(n_temp);
    
    index_inv <- which(n_temp == n);
    M_temp <- matrix(0, S_temp, S_temp);
    
    if(inv_values){
      r <- FW_inv$r_inv;
    }else{
      beta <- (1/(2*C)) - 1;
      r <- rbeta(1, 1, beta) * n;
      FW_inv$r_inv <- r;
    }
    
    r_temp <- append(FW_inv$radius, r, after=index_inv-1);
    
    #centre of the feeding range for the new species
    if(inv_values){
      c <- FW_inv$c_inv;
    }else{
      c <- runif(1, r/2, min(n, (1-r/2)));
      FW_inv$c_inv <- c;
    }
    c_temp <- append(FW_inv$centre, c, after=index_inv-1);
    
    for(i in 1:S_temp){
      offset <- r_temp[i]/2;     
      upper_bound <- c_temp[i] + offset;
      lower_bound <- c_temp[i] - offset;
      for(j in 1:S_temp){
        if(n_temp[j] > lower_bound & n_temp[j] < upper_bound)
          M_temp[j,i] = 1;
      } 
    }
    
    #########################################################################################    
    #after introducing the species we verify that the invasive species is not a duplicate of 
    #one of the species already present
    preys_i <- M_temp[,index_inv];
    predators_i <- M_temp[index_inv,];
    for(j in 1:S_temp){
      if(dup_sp){
        break;
      }
      if(index_inv == j) next;
      sim_prey <- preys_i == M_temp[,j];
      sim_preds <- predators_i == M_temp[j,];
      
      if(sum(sim_prey) == S_temp && sum(sim_preds) == S_temp ){
        dup_sp <- TRUE;
        print("InvadeNetwork :: the invasive is a duplicate of a species already present in the network");
        
        inv_values <- FALSE;
        
        break;
      }   
    }
    
    if(dup_sp) next;
    #########################################################################################
    
    #we verify that the network does not have cycles and that it is connected
    M <- M_temp;
    diag(M) <- 0;
    
    graf <- graph.adjacency(M);
    
    ##### When doing temperature dependent experiments I realised this constraint is too
    ##### stringent because after the initial dynamics for many values of T and K you might
    ##### end up with disconnected modules
    
    # y <- igraph::no.clusters(graf, mode='weak');
    # if(y > 1){
    #   print("InvadeNetwork :: the invaded network is not connected...");
    #   next;
    # }
    
    ##### so, I changed it to make sure that just the invasive species is connected
    
    if(degree(graf, index_inv, mode='all') > 0){
      y <- 1;
    }else{
      y <- 2;
      inv_values <- FALSE
      next;
    }
    
    
    clts <- igraph::clusters(graf, mode='strong');
    x <- clts$no;
    if(x < S_temp){
      clts_nos <- which(clts$csize > 1);  
      cycles_no_input <- FALSE;
      
      for(c in clts_nos){
        members <- which(clts$membership == c);
        cluster_ok <- FALSE;
        for(m in members){
          prey <- neighbors(graf, m, mode='in');
          if( length(intersect(prey, members))  > length(members) + 1 ){
            ## if this happens, this cycle/cluster has external energy input
            cluster_ok <- TRUE;
            break;
          }
        }
        
        if(!cluster_ok){
          print("InvadeNetwork :: the invaded network has cycles with no external energy input...");
          cycles_no_input <- TRUE;
          break;
        }
      }
      
      if(cycles_no_input){
        next;
      }else{
        x <- S_temp;
      }
      
    }
    
    #we check for whether the invasive is a basal species with cannibalistic links
    if(sum(preys_i) == 1 && M_temp[index_inv,index_inv] == 1){
      print("InvadeNetwork :: the invasive is a basal species with a cannibalistic link... removing it");
      M_temp[index_inv,index_inv] <- 0;
    } 
    
  }
  
  if(!inv_values & !is.null(FW_inv$n_inv)){
    #if this happens it means that we have had to obtain new values for the fundamental niche of
    #the species and hence we have to change the original value in the niche array for the new
    #one in order to be able to draw the network accurately further down
    
    FW_inv$original_niche[which(FW_inv$original_niche == FW$n_inv)] <- FW_inv$n_inv;
    
  }
  
  FW_inv$S <- S_temp;
  FW_inv$C <- sum(M_temp)/S_temp**2;
  
  FW_inv$M <- M_temp;
  FW_inv$niche <- n_temp;
  FW_inv$radius <- r_temp;
  FW_inv$centre <- c_temp;
  
  print('InvadeNetwork :: about to enter trophic positions');
  
  # print(FW$S);
  print(paste('InvadeNetwork :: invader index', index_inv));
  #print(FW_inv$M);
  
  FW_inv$invasive$index <- index_inv;
  FW_inv$invasive$original_index <- index_inv;
  FW_inv <- TrophicPositions(FW_inv);
  
  if(FW_inv$TP[index_inv] <= 1.0){
    print('InvadeNetwork :: the invasive is a basal species');
    FW_inv$invasive$producer <- TRUE;
  }else{
    print('InvadeNetwork :: the invasive is a consumer');
    FW_inv$invasive$producer <- FALSE;
  }
  
  #biomass of the invader (drawn randomly from an uniform distribution)
  #Xt <- runif(1, 0.05, 1);
  
  #init_b is the initial biomass of the invader species
  
  print(paste('the initial biomass of the invader is:', init_b));
  
  FW_inv$Initial <- append(FW_inv$Final, init_b, after=index_inv-1);
  FW_inv$Final <- append(FW_inv$Final, init_b, after=index_inv-1);
  
  # print(FW$Initial);
  # print(FW$Final);
  # 
  # print(FW_inv$Initial);
  # print(FW_inv$Final);
  
  return(FW_inv);

}


StructuralInvasion <- function(FW, gen_frac, vul_frac){
  x <- 0;
  FW_inv <- FW;
  S_temp <- FW_inv$S + 1;
  y <- S_temp;
  index_inv <- S_temp;
  dup_sp <- FALSE;
  
  while(x < S_temp | y > 1 | dup_sp){
    dup_sp <- FALSE;
    
    aprey <- 1.0-gen_frac;
    preys <- as.matrix(sample(0:1,FW_inv$S+1,TRUE,prob=c(aprey, gen_frac)));
    
    apred <- 1.0-vul_frac;
    predators <- t(as.matrix(sample(0:1,FW_inv$S,TRUE,prob=c(apred, vul_frac))));
    
    M_temp <- cbind(rbind(FW_inv$M, predators), preys);
    
    #########################################################################################    
    #after introducing the species we verify that the invasive species is not a duplicate of 
    #one of the species already present
    preys_i <- M_temp[,index_inv];
    predators_i <- M_temp[index_inv,];
    for(j in 1:S_temp){
      if(dup_sp){
        break;
      }
      if(index_inv == j) next;
      sim_prey <- preys_i == M_temp[,j];
      sim_preds <- predators_i == M_temp[j,];
      
      if(sum(sim_prey) == S_temp && sum(sim_preds) == S_temp ){
        dup_sp <- TRUE;
        print("StructuralInvasion :: the invasive is a duplicate of a species already present in the network");
        break;
      }   
    }
    
    if(dup_sp) next;
    #########################################################################################
    
    #we verify that the network does not have cycles and that it is connected
    M <- M_temp;
    diag(M) <- 0;
    
    graf <- graph.adjacency(M);
    y <- igraph::no.clusters(graf, mode='weak');
    if(y > 1){
      print("InvadeNetwork :: the invaded network is not connected...");
      next;
    }
    
    clts <- igraph::clusters(graf, mode='strong');
    x <- clts$no;
    if(x < S_temp){
      clts_nos <- which(clts$csize > 1);  
      cycles_no_input <- FALSE;
      
      for(c in clts_nos){
        members <- which(clts$membership == c);
        cluster_ok <- FALSE;
        for(m in members){
          prey <- neighbors(graf, m, mode='in');
          if( length(intersect(prey, members))  > length(members) + 1 ){
            ## if this happens, this cycle/cluster has external energy input
            cluster_ok <- TRUE;
            break;
          }
        }
        
        if(!cluster_ok){
          print("InvadeNetwork :: the invaded network has cycles with no external energy input ...");
          cycles_no_input <- TRUE;
          break;
        }
      }
      
      if(cycles_no_input){
        next;
      }else{
        x <- S_temp;
      }
      
    }
    
    #we check for whether the invasive is a basal species with cannibalistic links
    if(sum(preys_i) == 1 && M_temp[index_inv,index_inv] == 1){
      print("the invasive is a basal species with a cannibalistic link... removing it");
      M_temp[index_inv,index_inv] <- 0;
    } 
    
  }
   
  FW_inv$S <- S_temp;
  FW_inv$C <- sum(M_temp)/S_temp**2;  
  FW_inv$M <- M_temp;
  
  FW_inv$invasive$index <- index_inv;
  FW_inv <- TrophicPositions(FW_inv);
  
  if(FW_inv$TP[index_inv] <= 1.0){
    print('StructuralInvasion :: the invasive is a basal species');
    FW_inv$invasive$producer <- TRUE;
  }else{
    print('StructuralInvasion :: the invasive is a consumer');
    FW_inv$invasive$producer <- FALSE;
  }
  
  #biomass of the invader (drawn randomly from an uniform distribution)
  Xt <- runif(1, 0.05, 1);
  FW_inv$Initial <- append(FW_inv$Final, Xt, after=index_inv-1);
  FW_inv$Final <- append(FW_inv$Final, Xt, after=index_inv-1);
  
  return(FW_inv);
}

 
ObtainInvasiveTraits <- function(FW, FW_inv, invasion_type='neutral', efficiency_type='neutral', bs_factor=1.0, hill_exp=1.2){
  #we specify the body size of the invasive species
  #it could be specified as a parameter
  
  index_inv <- FW_inv$invasive$index;
  
  #recalculate these things...
  if(FW_inv$TP[index_inv] <= 1){  #body size of basal species == .01
    FW_inv$invasive$BS <- 0.01;
    FW_inv$BS <- append(FW_inv$BS, FW_inv$invasive$BS, after=index_inv-1);
    
    #we keep the same growth rate maximum...
    FW_inv$gr <- append(FW_inv$gr, mp['d','r']*(FW_inv$invasive$BS^mp['b','r'])*exp((mp['E','r']*temp_exp)), after=index_inv-1);
    
    #hill exponent
    # FW_inv$h <- append(FW_inv$h, hill_exp, after=index_inv-1);
    
    FW_inv$x <- append(FW_inv$x, 0, after=index_inv-1); 
    FW_inv$K <- append(FW_inv$K, mp['d','K']*(FW_inv$invasive$BS^mp['b','K'])*exp((mp['E','K']*temp_exp)), after=index_inv-1);
    
  }else{
    #average predator-prey mass ratios
    #first we obtain the current ('real') average of the predator-prey mass ratios
    ratios <- c();
    for(i in 1:FW$S){
      for(j in 1:FW$S){
        if(i==j)
          next;
        if(FW$M[i,j] == 1){
          ratios <- append(ratios, FW$BS[j]/FW$BS[i]);
        }
        
      }
    }
    
    FW$ratios <- ratios;
    ppmr <- mean(ratios);
    FW$real_ppmr <- ppmr;
    
    FW_inv$invasive$BS <- 0.01*(FW_inv$ppmr^( ((FW_inv$TP[i]-1))));
    
    if(invasion_type == 'better'){
      FW_inv$invasive$BS <- FW_inv$invasive$BS*bs_factor;
    }else if(invasion_type == 'worse'){
      FW_inv$invasive$BS <- FW_inv$invasive$BS/bs_factor;
    } 
    
    FW_inv$BS <- append(FW_inv$BS, FW_inv$invasive$BS , after=index_inv-1);
    
    
    #growth rate maximum
    FW_inv$gr <- append(FW_inv$gr, 0, after=index_inv-1);
    #hill exponent
    # FW_inv$h <- append(FW_inv$h, hill_exp, after=index_inv-1);
    
    
    FW_inv$x <- append(FW_inv$x, mp['d','x']*(FW_inv$invasive$BS^mp['b','x'])*exp((mp['E','x']*temp_exp)), after=index_inv-1); 
    FW_inv$K <- append(FW_inv$K, 0, after=index_inv-1);
  }
  
 
  FW_inv$e <- 0.85;
  
  S <- dim(FW_inv$M)[1]
  ##### here we calculate the interaction matrix of alphas and handling times
  FW_inv$alphas <- FW_inv$M
  FW_inv$Th <- FW_inv$M
  for(i in 1:S){
    for(j in 1:S){
      if(FW_inv$M[i,j] == 1){   #i is eaten by j
        FW_inv$alphas[i,j] <- mp['d','alpha']*FW_inv$BS[i]^mp['b','alpha']*FW_inv$BS[j]^mp['c','alpha']*exp(mp['E','alpha']*temp_exp);
        FW_inv$Th[i,j] <- mp['d','Th']*FW_inv$BS[i]^mp['b','Th']*FW_inv$BS[j]^mp['c','Th']*exp(mp['E','Th']*temp_exp);
      }
    }
  }
  
  return(FW_inv);
}


PerformInvasionExperiments <- function(filename=NULL, FW, type='niche'){
  print(paste('PerformInvasionExperiment :: file = ', filename));
  
  # Tmax <- 4000;
  # tlength <- 200;
  #obtain the network to be invaded
  if(! is.null(filename)){
    idx <- which(strsplit(filename, '')[[1]]=='C');
    
    original_C <- as.numeric(substr(filename, start=idx+2, stop=idx+6));
  
    print(paste("PerformInvasionExperiment :: connectance = ", original_C));
    load(filename);
    
    if(is.null(FW$dynamics)){
      FW <- ObtainDynamics(FW);
      FW <- RecalculateFW(FW);
    }
    
    FW <- RecalculateFW(FW);
    
    exec <- FW$dynamics;
    FW$total_biomass_b4inv <- sum(FW$Final);
    
    ##here we calculate the statistical properties of the original network
    FW$orig_fw <- CalculateFoodWebStats(FW$orig_fw);
    FW <- CalculateFoodWebStats(FW);
    
    if(is.null(FW$colours)){
      print('I am generating the colours')
      FW$colours <- rainbow(FW$orig_fw$S);
    }
    
  }else{
    exec <- FW$dynamics;
    FW$total_biomass_b4inv <- sum(FW$Final);
    
    ##here we calculate the statistical properties of the original network
    FW$orig_fw <- CalculateFoodWebStats(FW$orig_fw);
    FW <- CalculateFoodWebStats(FW);
    
    print(paste('species after recalculation using the same Final obtained in the previous step', FW$S));
    
    FW$colours <- rainbow(FW$orig_fw$S);
    
    #and then save the RData and plot the output.
    save(FW, file=paste('./output/RFile-FW-S-',FW$S,'-C-', round(FW$C, 3), '-fwno-', FW$fwno ,'-before-inv-temp-', FW$temp,'-eutroph-', FW$Eutroph, '.RData', sep=''));
    
    # print("PerformInvasionExperiment :: plot the dynamics");
    # # pdf(paste('Plot-dynamics-S',FW$S,'-C-', round(FW$C, 3), '-fwno-', FW$fwno, '-before-inv-stable-', FW$stability_a$State, '.pdf', sep=''), width=10, height=4);
    # par(mar=c(4,4,1,1));
    # plot(0,0,type='n', xlab='Time', ylab='Species biomass', xlim=range(exec[,1]), ylim=range(exec[,-1])); #, ylim=range(exec[,-1])
    # 
    # for(k in colnames(exec)){
    #   if(k != 'time'){
    #     index <- which(colnames(exec) == k);
    # 
    #     lines(exec[,1], exec[,index], col=FW$colours[as.integer(k)], lwd=1);
    #   }
    # }

    #legend("topleft", paste(c("Initial = ", "Extinctions = "), c(S, extinct)), bty='n');
    # dev.off();
    original_C <- FW$original_C;
  }
  
  
  #we proceed to obtain the invaded network....
  if(type == 'structural'){
    invaded <- StructuralInvasion(FW, 0.4, 0.2);
    while(invaded$invasive$producer)
      invaded <- StructuralInvasion(FW, 0.4, 0.2);
  }else{
    invaded <- InvadeNetwork(FW, original_C);
    # while(invaded$invasive$producer)
    #   invaded <- InvadeNetwork(FW, original_C);
  }
  
  invaded <- CalculateInvaderProperties(invaded);
  invaded <- CalculateFoodWebStats(invaded);
  
  ##neutral invasion
  
  invaded_temp <- ObtainInvasiveTraits(FW, invaded);
  invaded_temp$experiment_type <- 'neutral';
  
  after_inv <- GetNetworkOutput(invaded_temp, exec, files_suffix=paste0('neutral-temp-', FW$temp,'-eutroph-', FW$Eutroph));
  
  # DrawNetworks(FW, invaded_temp, after_inv, files_suffix='neutral');
  
  ###########################################################################################################
  ###efficiency experiments...
  # eff_types <- c('better', 'worse');
  # for(j in 1:length(eff_types)){
  #   invaded_temp <- ObtainInvasiveTraits(FW, invaded, efficiency_type=eff_types[j]);
  #   
  #   invaded_temp$experiment_type <- 'efficiency';
  #   invaded_temp$eff_type <- eff_types[j];
  #   
  #   after_inv <- GetNetworkOutput(invaded_temp, Tmax, tlength, exec, files_suffix=paste('efficiency-', eff_types[j], sep=''));
  #   
  #   DrawNetworks(FW, invaded_temp, after_inv, files_suffix=paste('efficiency-', eff_types[j], sep=''));
  #   
  # }
  ###########################################################################################################
  
  ###########################################################################################################
  #body size experiments
  #the neutral experiment...
  
  #invaded$experiment_type <- 'body_size';
  #invaded$inv_type <- 'neutral';
    
  #invaded_temp <- ObtainInvasiveTraits(FW, invaded, invasion_type='neutral');
  #after_inv <- GetNetworkOutput(invaded_temp, Tmax, tlength, exec);
  
  #DrawNetworks(FW, invaded_temp, after_inv);
  
  #and, for the other two...
  inv_types <- c('better', 'worse');
  bs_fractions <- c(2, 5, 10);
  for(i in 1:length(inv_types)){
    for(j in bs_fractions){
      invaded_temp <- ObtainInvasiveTraits(FW, invaded, invasion_type=inv_types[i], bs_factor=j);
      
      invaded_temp$experiment_type <- 'body_size';
      invaded_temp$inv_type <- inv_types[i];
      invaded_temp$bs_fraction <- j;
      
      after_inv <- GetNetworkOutput(invaded_temp, exec, files_suffix=paste('body_size-', inv_types[i],'-', j,'-temp-', FW$temp,'-eutroph-', FW$Eutroph, sep=''));
    
      # DrawNetworks(FW, invaded_temp, after_inv, files_suffix=paste('body_size-', inv_types[i],'-', j, sep=''));
      
    }
  }
  ###########################################################################################################
  
  ###########################################################################################################
  ##biomass experiments....
  init_bios <- c(0.001, 0.01, 0.1);
  for(b in init_bios){

    invaded_temp <- ObtainInvasiveTraits(FW, invaded);

    invaded_temp$experiment_type <- 'biomass';
    invaded_temp$initial_biomass <- b;

    ##applying the biomass
    invaded_temp$Initial[invaded_temp$invasive$index] <- b;
    invaded_temp$Final[invaded_temp$invasive$index] <- b;

    after_inv <- GetNetworkOutput(invaded_temp, exec, files_suffix=paste('biomass-', b,'-temp-', FW$temp,'-eutroph-', FW$Eutroph, sep=''));

    # DrawNetworks(FW, invaded_temp, after_inv, files_suffix=paste('biomass-', b, sep=''));
  }

  
  ###########################################################################################################
  
  
  ###########################################################################################################
  ##hill exponents experiments....
  # hill_exponents <- c(1.0, 2.0);
  # for(h in hill_exponents){
  #   invaded_temp <- ObtainInvasiveTraits(FW, invaded, hill_exp=h);
  #   
  #   invaded_temp$experiment_type <- 'hill_exponent';
  #   invaded_temp$hill_exp <- h;
  #   
  #   after_inv <- GetNetworkOutput(invaded_temp, Tmax, tlength, exec, files_suffix=paste('hill-', h, sep=''));
  #   
  #   DrawNetworks(FW, invaded_temp, after_inv, files_suffix=paste('hill-', h, sep=''));
  # }
  
  ###########################################################################################################
  
  
  ###########################################################################################################
  # ##generalism experiments....
  # index_inv <- invaded$invasive$index;
  # struc_fractions <- c(0.25, 0.5, 0.75);
  # struc_actions <- c('add', 'remove');
  # #we perform these set of experiments if the invader has at least 4 links in that we can remove
  # 
  # print(paste('invader index', index_inv));
  # 
  # print(paste('number of prey of the invader =', sum(invaded$M[,index_inv])));
  # 
  # if( sum(invaded$M[,index_inv]) >= 4 ){
  #   print('I am in the rewiring experiments... for gen');
  #   
  #   for(a in struc_actions){
  #     for(f in struc_fractions){
  #       invaded_temp <- ChangeInvaderLinks(invaded, type='gen', fraction=f, action=a);
  #       invaded_temp <- ObtainInvasiveTraits(FW, invaded_temp);
  #       
  #       invaded_temp$experiment_type <- 'generalism';
  #       invaded_temp$generalism_type <- a;
  #       invaded_temp$generalism_frac <- f;
  #       
  #       after_inv <- GetNetworkOutput(invaded_temp, Tmax, tlength, exec, files_suffix=paste('gen-',a, '-', f, sep=''));
  #       
  #       #DrawNetworks(FW, invaded_temp, after_inv, files_suffix=paste('gen-',a, '-', f, sep=''));
  #     }
  #   }
  #     
  # }
   
  ###########################################################################################################
  
  ###########################################################################################################
  ##vulnerability experiments....
  
  #for vulnerability experiments our only hypothesis is that predation release favours the success of invasive
  # species, and so we only look at removing predatory links from it
  
  #we make sure the invasive species has at least 4 predatory links since we do not want it to remain without
  #predators even for the largest value of the removal fraction (0.75) which for a species with 4 links will
  #leave only one predator, for a species with less than 4 links it will make it become a top predator
  
  # print(paste('number of predators of the invader =', sum(invaded$M[index_inv,])));
  # 
  # if( sum(invaded$M[index_inv,]) >= 4 ){
  #   print('I am in the rewiring experiments... for vul');
  #   
  #   for(f in struc_fractions){
  #     invaded_temp <- ChangeInvaderLinks(invaded, type='vul', fraction=f, action='remove');
  #     invaded_temp <- ObtainInvasiveTraits(FW, invaded_temp);
  #     
  #     invaded_temp$experiment_type <- 'vulnerability';
  #     invaded_temp$vulnerability_type <- 'remove';
  #     invaded_temp$vulnerability_frac <- f;
  #     
  #     after_inv <- GetNetworkOutput(invaded_temp, Tmax, tlength, exec, files_suffix=paste('vul-remove-', f, sep=''));
  #     
  #     #DrawNetworks(FW, invaded_temp, after_inv, files_suffix=paste('vul-remove-', f, sep=''));
  #   }  
  # }
  # 
  ###########################################################################################################
  
}
  

ChangeInvaderLinks <- function(invaded, type='gen', fraction=0.5, action='remove'){
  index_inv <- invaded$invasive$index;
  M_temp <- invaded$M;
  S <- dim(M_temp)[1];
  
  inv_tl <- invaded$TrLevels[index_inv];
  similar_sps <- setdiff(which(invaded$TrLevels == inv_tl), index_inv);
  dup_sp <- FALSE;
  
  if(type == 'gen'){
    preys_inv <- M_temp[,index_inv];
    
    n_links <- ceiling(sum(preys_inv)*fraction);
    if(action == 'remove'){
      if(n_links >= sum(preys_inv)){
        n_links <- sum(preys_inv) - 1;     
      }
      
      if(n_links == 0){
        print('ChangeInvaderLinks:: No links to resources removed');
        return(invaded);
      }
      
      l_removed <- 0;
      
      while(l_removed < n_links){
        for(sp in 1:S){
          
          if(preys_inv[sp] == 1 & M_temp[sp,index_inv] == 1 & runif(1,0,1) <= 0.5){
            M_temp[sp,index_inv] <- 0;
            l_removed <- l_removed + 1;
          }
          if(l_removed == n_links){break;}
        }
      }
      
      print(paste('ChangeInvaderLinks::', l_removed, 'links to resources were removed'));
      
    }
    
    if(action == 'add'){
      possible_prey <- rep(0,S);
      for(ssp in similar_sps){
        prey <- M_temp[,ssp];
        possible_prey <- possible_prey + prey;
      }
      
      possible_prey[which(possible_prey > 1)] <- 1;
      possible_prey <- possible_prey - preys_inv;
      possible_prey[which(possible_prey < 0)] <- 0;
      
      possible_prey[index_inv] <- 0;    #we don't want to add cannibalistic links
      
      if(n_links > sum(possible_prey)){
        for(sp in 1:S){
          if(possible_prey[sp] == 1){
            M_temp[sp,index_inv] <- 1;  
          }
        }
        print(paste('ChangeInvaderLinks::', sum(possible_prey), 'links to resources were added'));
      }else{
        l_added <- 0;
        while(l_added < n_links){
          for(sp in 1:S){
            if(sp == index_inv) next;   #we don't want to add cannibalistic links
            if(possible_prey[sp] == 1 &  M_temp[sp,index_inv] == 0 & runif(1,0,1) <= 0.5){
              M_temp[sp,index_inv] <- 1;
              l_added <- l_added + 1;
            }
            if(l_added == n_links){break;}
          }
        }
        
        print(paste('ChangeInvaderLinks::', l_added, 'links to resources were added'));
      }
    }
  }
  
  if(type == 'vul'){
    predators_inv <- M_temp[index_inv,];
    n_links <- ceiling(sum(predators_inv)*fraction);
    
    if(n_links == 0){
      print('ChangeInvaderLinks:: No links to consumers removed');
      return(invaded);
    }
    
    if(action == 'remove'){
      if(n_links >= sum(predators_inv)){
        M_temp[index_inv, ] <- 0;
        print('ChangeInvaderLinks:: All links to consumers were removed');
      }else{
        l_removed <- 0;
        
        while(l_removed < n_links){
          for(sp in 1:S){
            if(predators_inv[sp] == 1 & M_temp[index_inv, sp] == 1 & runif(1,0,1) <= 0.5){
              M_temp[index_inv, sp] <- 0;
              l_removed <- l_removed + 1;
            }
            if(l_removed == n_links){break;}
          }
        }
        print(paste('ChangeInvaderLinks::', l_removed, 'links to consumers were removed'));
      }
    }
    
    if(action == 'add'){
      possible_predators <- rep(0,S);
      for(ssp in similar_sps){
        predators <- M_temp[ssp,];
        possible_predators <- possible_predators + predators;
      }
      
      possible_predators[which(possible_predators > 1)] <- 1;
      possible_predators <- possible_predators - predators_inv;
      possible_predators[which(possible_predators < 0)] <- 0;
      
      possible_predators[index_inv] <- 0;   #make sure we don't add cannibalistic links...
      
      if(n_links > sum(possible_predators)){
        for(sp in 1:S){
          if(possible_predators[sp] == 1){
            M_temp[index_inv,sp] <- 1;  
          }
        }
        print(paste('ChangeInvaderLinks::', sum(possible_predators), 'links to consumers were added'));
      }else{
        
        l_added <- 0;
        while(l_added < n_links){
          for(sp in 1:S){
            if(possible_predators[sp] == 1 & M_temp[index_inv,sp] == 0 & runif(1,0,1) <= 0.5){
              M_temp[index_inv,sp] <- 1;
              l_added <- l_added + 1;
            }
            if(l_added == n_links){break;}
          }
        }
        print(paste('ChangeInvaderLinks::', l_added, 'links to consumers were added'));
      }
    }  
  }
  
  invaded$M <- M_temp;
  invaded <- TrophicPositions(invaded);
  return(invaded);
  
}

GetNetworkOutput <- function(invaded, original_dynamics, files_suffix=''){
  
  invaded$orig_fw <- invaded;
  
  #invaded$initS <- invaded$S;  
  invaded <- ExecuteModelFromFoodWeb(invaded, TRUE);
  
  invaded <- RecalculateFW(invaded);
  
  invaded$total_biomass_afterinv <- sum(invaded$Final);
  # invaded <- GetStabilityMeasures(invaded, Tmax, tlength);

  if(invaded$collapsed){
    print('GetNetworkOutput :: The network collapsed - 0 species left');
    print('GetNetworkOutput :: No output');
  }
  
  # if(invaded$stability_a$State == 0){
  #   print("GetNetworkOutput :: Invaded network not stable");
  # }else{
  #   print("GetNetworkOutput :: Invaded network is stable");
  # }
  
  if(sum(invaded$Final == 0)){
    print('GetNetworkOutput :: The network collapsed after the introduction of the invasive species')
    invaded$collapsed <- TRUE;
    #return;
  }
  
  
  print(paste('GetNetworkOutput :: invader index before =', invaded$invasive$index));
  invaded$invasive$index <- which(invaded$niche == invaded$n_inv);
  print(paste('GetNetworkOutput :: invader index after =', invaded$invasive$index));
  
  if(invaded$n_inv %in% invaded$niche){
    invaded$invasive$extinct <- FALSE;  
  }else{
    invaded$invasive$extinct <- TRUE;
    print('The invasive species went extinct');
  }
  
  invaded$S_lost_aftinv <- invaded$extinct_species;
  invaded$S_lost_aftinv$number <- invaded$orig_fw$S - invaded$S;
  
  #and we also calculate the statistical properties of the network...
  invaded <- CalculateFoodWebStats(invaded);
  
  #exec <- invaded$dynamics;
  
  print("GetNetworkOutput :: plot the dynamics for the invaded network...");
  print(files_suffix);
  # pdf(paste('Plot-dynamics-S',invaded$S,'-C-', round(invaded$C, 3), '-fwno-', invaded$fwno, '-after-inv-', files_suffix, '.pdf', sep=''), width=10, height=4);
  # par(mar=c(4,4,1,1));
  # plot(0,0,type='n', xlab='Time', ylab='Species biomass', xlim=range(invaded$dynamics[,1]*2), ylim=range(max(invaded$dynamics[,-1]), min(invaded$dynamics[,-1]), max(original_dynamics[,-1]), min(original_dynamics[,-1])  )); # ylim=range(invaded$dynamics[,-1]),
  # 
  # 
  # #colours <- rainbow((dim(original_dynamics)[2]-1));
  # #for(k in 1:(dim(original_dynamics)[2]-1)){
  # for(k in colnames(original_dynamics)){
  #   if(k != 'time'){
  #     index <- which(colnames(original_dynamics) == k);
  # 
  #     lines(original_dynamics[,1], original_dynamics[,index], col=invaded$colours[as.integer(k)], lwd=1);
  #   }
  # }
  # 
  # offset <- max(original_dynamics[,1]);
  # 
  # #we plot the invasive species dynamics in black
  # invaded$colours <- append(invaded$colours, '#000000FF');
  # for(k in colnames(invaded$dynamics)){
  #   if(k != 'time'){
  #     index <- which(colnames(invaded$dynamics) == k);
  # 
  #     if(k == ''){
  #       invader_dynamics <- invaded$dynamics[,index];
  #       next;
  #     }
  # 
  #     colour <- invaded$colours[as.integer(k)];
  #     lines((invaded$dynamics[,1])+offset, invaded$dynamics[,index], col=colour, lwd=1);
  #   }
  # 
  # }
  # 
  # if(exists('invader_dynamics')){
  #   colour <- invaded$colours[length(invaded$colours)];
  #   lines((invaded$dynamics[,1])+offset, invader_dynamics, col=colour, lwd=1);
  # }

#   #legend("topleft", paste(c("Initial = ", "Extinctions = "), c(S, extinct)), bty='n');
#   dev.off();
  
  #now that we have the new index of the invader, we can calculate its properties
  if(!invaded$invasive$extinct) invaded <- CalculateInvaderProperties(invaded);
  
  #and then save the RData and plot the output.
  save(invaded, file=paste('./output/RFile-FW-S-',invaded$S,'-C-', round(invaded$C, 3), '-fwno-', invaded$fwno, '-after-inv-', files_suffix, '.RData', sep=''));
  
  return(invaded);
  
}

CalculateInvaderProperties <- function(invaded){
  invaded$invasive$TP <- invaded$TP[invaded$invasive$index];
  return(invaded);
}

DrawNetworks <- function(FW, invaded, after_inv, files_suffix=''){
  
  # pdf(paste('Networks-S', after_inv$S,'-C-', round(after_inv$C, 3), '-fwno-', after_inv$fwno, '-after-inv-', files_suffix, '.pdf', sep=''), width=60, height=20);
  
  par(mai=c(.2,.2,.2,.2));
  layout(matrix(1:4, nr=1, byrow=TRUE));
  layout.show(4);
  
  ########################################################################
  #we first draw the original network as it came out of the niche model...
  vals <- FW$original_niche;
  vals[!(vals %in% FW$orig_fw$niche)] <- 0;
  v_labels <- (vals > 0);
  v_labels <- (1:length(vals))[v_labels];
  
  g <- graph.adjacency(FW$orig_fw$M);
  lay <- GetRankLayout(FW$orig_fw);
  plot(g, vertex.label=v_labels, vertex.size=6, edge.color="black",vertex.color="red", layout=lay, vertex.label.cex=2.5);
  
  ########################################################################
  #then we draw the network after the initial dynamics have been simulated
  #and the species that are lost before the invasion occurs are removed
  #from the network

  vals <- FW$original_niche;
  vals[!(vals %in% FW$niche)] <- 0;
  v_labels3 <- (vals > 0);
  v_labels3 <- (1:length(vals))[v_labels3];
  
  g1 <- graph.adjacency(FW$M);
  lay1 <- GetRankLayout(FW);
  plot(g1, vertex.label=v_labels3, vertex.size=6, edge.color="black",vertex.color="red", layout=lay1, vertex.label.cex=2.5);
  
  ########################################################################
  #the third network is the one with the new invasive species
  
  vals <- FW$original_niche;
  vals[!(vals %in% invaded$niche)] <- 0;
  v_labels1 <- (vals > 0);
  v_labels1 <- (1:length(vals))[v_labels1];
  
  g2 <- graph.adjacency(invaded$M);
  lay2 <- GetRankLayout(invaded);
  plot(g2, vertex.label=v_labels1, vertex.size=6, edge.color="black",vertex.color="red", layout=lay2, vertex.label.cex=2.5);
  
  ########################################################################
  #and finally... we draw the network that results after the dynamics
  #have been simulated after invasion
  
  vals <- FW$original_niche;
  vals[!(vals %in% after_inv$niche)] <- 0;
  v_labels2 <- (vals > 0);
  v_labels2 <- (1:length(vals))[v_labels2];
  
  g3 <- graph.adjacency(after_inv$M);
  lay3 <- GetRankLayout(after_inv);
  plot(g3, vertex.label=v_labels2, vertex.size=6, edge.color="black",vertex.color="red", layout=lay3, vertex.label.cex=2.5);
  
  # dev.off();
  
  # pdf(paste('DegreeDistributions-S', after_inv$S,'-C-', round(after_inv$C, 3), '-fwno-', after_inv$fwno, '-after-inv-', files_suffix, '.pdf', sep=''), width=12, height=3.5);
  
  layout(matrix(1:4, nr=1, byrow=TRUE));
  layout.show(4);
  
  plot(FW$orig_fw$stats$dd, log="xy", xlab="degree", ylab="cumulative frequency", col=1, main="degree distributions");
  points(FW$orig_fw$stats$in_dd, col=2, pch=2);
  points(FW$orig_fw$stats$out_dd, col=3, pch=3);
  legend(1, .015, c('sum','in','out'), col=1:3, pch=1:3, ncol=1, yjust=0, lty=0);
  
  plot(FW$stats$dd, log="xy", xlab="degree", ylab="cumulative frequency", col=1, main="degree distributions");
  points(FW$stats$in_dd, col=2, pch=2);
  points(FW$stats$out_dd, col=3, pch=3);
  legend(1, .015, c('sum','in','out'), col=1:3, pch=1:3, ncol=1, yjust=0, lty=0);
  
  plot(invaded$stats$dd, log="xy", xlab="degree", ylab="cumulative frequency", col=1, main="degree distributions");
  points(invaded$stats$in_dd, col=2, pch=2);
  points(invaded$stats$out_dd, col=3, pch=3);
  legend(1, .015, c('sum','in','out'), col=1:3, pch=1:3, ncol=1, yjust=0, lty=0);
  
  plot(after_inv$stats$dd, log="xy", xlab="degree", ylab="cumulative frequency", col=1, main="degree distributions");
  points(after_inv$stats$in_dd, col=2, pch=2);
  points(after_inv$stats$out_dd, col=3, pch=3);
  legend(1, 0.015, c('sum','in','out'), col=1:3, pch=1:3, ncol=1, yjust=0, lty=0);
  
  # dev.off();
  
}

ExecuteModelFromFoodWeb <- function(FW, invaded=FALSE){
  Tmax <- (365*24*60*60)*600;
  tlength <- (60*60*24*10)    ### 60*60*2.4 this timestep is one tenth of a day in seconds... after conversations with Arnaud
  
  if(!invaded){
    #obtain the body sizes for the species in the network
    FW <- ObtainBodySize(FW);
    #obtain the parameters for the allometric model
    FW <- ObtainATNParams(FW);
    
    FW$orig_fw <- FW;
    
    #### here we set initial abundances
    Xt <- FW$K
    Xt[which(FW$K == 0)] <- rep(mean(FW$K[which(FW$K != 0)])/8, length(which(FW$K == 0)))
    
    FW$Initial <- Xt;
    FW$Final <- Xt;
  }
  
  exec <- RecalculateX(RunDynamics(Tmax, tlength, FW=FW));
  
  FW$dynamics <- exec;
  FW$Final <- exec[dim(exec)[1], -1];
  
  FW <- RecalculateFW(FW);
  
  if(!invaded){
    FW$S_lost_b4inv <- FW$extinct_species;
    FW$S_lost_b4inv$number <- FW$orig_fw$S - FW$S;  
  }
  
  return(FW);
}


############### TEMPERATURE DEPENDENT INVASIONS START HERE #####################

replicates <- 100;
start_time <- proc.time()
# ExperimentsProtocol <- function(){
  
  #### First thing we do is to create the food webs we are going to use accross experiments.
  #### The same food webs will be used for each combination of temperature and eutrophication

  #### These ensures that results are comparable to Binzer et al. 2016 GCB in which they did it like this

  food_webs <- vector('list', replicates)
  for(r in 1:replicates){
    fwno <- sample(1:5000, 1);
    set.seed(fwno);
    FW <- NicheNetwork(30, 0.1);
    
    FW$fwno <- fwno
    FW$orig_fw <- FW;
    
    food_webs[[r]] <- FW
  }

  require(foreach)
  require(doParallel)

  #setup parallel backend to use many processors
  # #cores <- detectCores()
  cl <- makeCluster(35) #cores[1]-1) #not to overload your computer
  registerDoParallel(cl)

  
  
  #### defining the reference temperature
  To <- 293.15;
  #### then, we specify the values of temperature and eutrophication that we want to test
  #### Temperatures
  temps <- seq(-20,20,1)
  #### carrying capacities
  ks <- c(10, 20)
  
  for(t_idx in 1:length(temps)){
    
    Temp_new <- To + temps[t_idx]
    temp_exp <- (To-Temp_new)/(boltz*Temp_new*To)
    
    print(Temp_new)
    
    for(k_idx in 1:length(ks)){
      mp['d', 'K'] <- ks[k_idx]
      
      print(mp['d', 'K'])
      
      # for(r in 1:replicates){
      persistence <- foreach(r = 1:replicates, .combine=append) %dopar%{
        require(igraph)
        require(deSolve)
        
        FW <- food_webs[[r]]
        FW <- ExecuteModelFromFoodWeb(FW);
        
        FW$temp <- Temp_new
        FW$Eutroph <- mp['d', 'K']
        print('transient dynamics finished')
        
        if(FW$S > 0){
          PerformInvasionExperiments(FW=FW, type='niche');
        }
        
        FW$S
        
      }
      
      print(persistence)
    }
  }

  stop_time <- proc.time()
  print(stop_time - start_time)
  
  stopCluster(cl)
# }


# ExperimentsProtocol();




















