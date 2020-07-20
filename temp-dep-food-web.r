library(igraph)
library(deSolve)

########## model parameters and constants ##############
model_parameters <- data.frame(K=c(5, 0.28, NA, 0.71), 
                               r=c(exp(-15.68), -0.25, NA, -0.84), 
                               alpha=c(exp(-13.1), 0.25, -0.8, -0.38), 
                               Th=c(exp(9.66), -0.45, 0.47, 0.26), 
                               x=c(exp(-16.54), -0.31, NA, -0.69))
rownames(model_parameters) <- c('d', 'b', 'c', 'E')

##### just to write less
mp <- model_parameters

boltz <- 8.617e-5     # the boltzmann constant

##### model parameters end here ######


#build the niche model according to:
# Richard J. Williams and Neo D. Martinez (2000). "Simple rules yield complex food webs". Nature.

NicheNetwork <- function(S, C){
  S <- S+1;
  x <- 0;
  y <- S;
  dup_sp <- FALSE;
  
  while(x < S | y > 1 | dup_sp){
    dup_sp <- FALSE;
    M <- matrix(0, S, S);
    
    # first we obtain the n value for all the species
    n <- sort(runif(S, 0, 1));
    
    #we then obtain the feeding range for each of the species, drawn from a beta
    #distribution with beta value = (1/(2*C)) - 1
    beta <- (1/(2*C)) - 1;
    r <- rbeta(S, 1, beta) * n;    
    
    #we enforce the species with the lowest niche value to be a basal species
    r[1] <- 0;
    
    #the centre of the feeding range is drawn from an uniform distribution
    #between ri/2 and min(ni, 1-ri/2)
    c <- numeric(S);
    for(i in 1:S){
      c[i] <- runif(1, r[i]/2, min(n[i], (1-r[i]/2)));
      
      #once we have the r's and the c's we can obtain the boundaries of the feeding range
      offset <- r[i]/2;
      upper_bound <- c[i] + offset;
      lower_bound <- c[i] - offset;
      
      for(j in 1:S){
        if(n[j] > lower_bound & n[j] < upper_bound)
          M[j,i] = 1;
      } 
    }
    
    #we verify that the network (i) is connected and (2) does not have cycles with no external energy input
    M_temp <- M;
    diag(M_temp) <- 0;
    
    graf <- graph.adjacency(M);
    y <- igraph::no.clusters(graf, mode='weak');
    if(y > 1){
      #print("NicheNetwork :: the network is not connected...");
      next;
    }
    
    clts <- igraph::clusters(graf, mode='strong');
    x <- clts$no;
    if(x < S){
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
          #print("NicheNetwork :: the network has cycles with no external energy input...");
          cycles_no_input <- TRUE;
          break;
        }
      }
      
      if(cycles_no_input){
        next;
      }else{
        x <- S;
      }
      
    }
    
    #and we also verify that there are not duplicate species
    for(i in 1:S){
      if(dup_sp){
        break;
      }
      preys_i <- M[,i];
      predators_i <- M[i,];
      for(j in 1:S){
        if(i == j) next;
        sim_prey <- preys_i == M[,j];
        sim_preds <- predators_i == M[j,];
        
        if(sum(sim_prey) == S && sum(sim_preds) == S ){
          dup_sp <- TRUE;
          #print("NicheNetwork :: there is a duplicate species");
          break;
        }   
      }
      
      #as long as we are traversing the network we might as well check
      #for basal species with cannibalistic links...
      
      if(sum(preys_i) == 1 && M[i,i] == 1){
        #print("NicheNetwork :: a cannibalistic link was found in a basal species... removing it");
        M[i,i] <- 0;
      } 
    }
  }
  
  #we keep a reference of the original niche as a way of idetifying the species in the network:
  #each species is going to be identified by the index of its niche value within the 
  #original_niche array.
  original_niche <- n;
  
  ########################################################################################
  #select a random species to be a possible invader and remove it from the current network
  #keeping the values of n, r and c for the invasion event
  
  #we select the possible invader at random
  inv_idx <- sample(1:S, 1);
  
  #enforce the invasive species to be a consumer
  
  #here we can specify whether the invasive species should have at least x links
  #so we can perform the rewiring experiments afterwards....
  
  ### this was disabled for temperature dependent simulations because we didn't want to
  ### constraint the connectivity of the invader
  #while(sum(M[,inv_idx]) < 1 || sum(M[inv_idx,]) < 1) inv_idx <- sample(1:S, 1);
  
  
  #print(paste('the invader has', sum(M[,inv_idx]), 'prey and', sum(M[inv_idx,]), 'predators'));
  
  M <- M[-inv_idx, -inv_idx];
  n_inv <- n[inv_idx];
  n <- n[-inv_idx];
  c_inv <- c[inv_idx];
  c <- c[-inv_idx];
  r_inv <- r[inv_idx];
  r <- r[-inv_idx];
  
  ########################################################################################
  
  return(list(S=dim(M)[1], C=sum(M)/((dim(M)[1])**2), M=M, niche=n, centre=c, radius=r, original_C=C, n_inv=n_inv, r_inv=r_inv, c_inv=c_inv, original_niche=original_niche));
  
}

NormalizeMatrix <- function(M){
  colsum_M <- colSums(M);
  colsum_M[colsum_M==0] <- 1;
  return(t(t(M)/colsum_M));
  
}

##### Calculate the trophic positions using the prey-averaged trophic level metric from:
##### Williams and Martinez (2004). American Naturalist.
##### Limits to trophic levels and omnivory in complex food webs: theory and data.

TrophicPositions <- function(FW){
  S <- dim(FW$M)[1];
  if(S < 3) return(FW);
  M <- NormalizeMatrix(FW$M);
  if(det(diag(S) - t(M)) != 0){
    FW$TP <- solve(diag(S)-t(M), rep(1,S));
  }else{
    tmp <- diag(S);
    for(i in 1:9){
      tmp <- tmp %*% t(M) + diag(S);
    }
    FW$TP <- tmp %*% rep(1,S);
  }
  FW$W <- M;
  FW <- TrophicLevels(FW);
  return(FW);
}


##### We can also obtain the trophic level of each species
##### e.g. basal, herbivore, omnivore, top predator
TrophicLevels <- function(FW){
  S <- dim(FW$M)[1];
  
  if(S < 2) return(FW)
  
  TrLevels <- rep(-1, S)
  
  #we find the neighbours of each node in the network
  M_temp <- FW$M;
  diag(M_temp) <- 0;
  
  #all the species with trophic position 1 are basal species, and therefore
  #belong to the trophic level 0
  TrLevels[FW$TP == 1.0] <- 0;
  producers <- which(TrLevels == 0);
  for(i in 1:S){
    if(TrLevels[i] == -1){
      herb <- TRUE;
      top <- TRUE;
      
      #for each species we verify two things:
      #first: if any of its prey does not belong to trophic level 0,
      #then it is not a herbivore
      if(sum(TrLevels[which(M_temp[,i] != 0)]) != 0) herb <- FALSE;
      
      #second: if it has any predators, then it is not a top predator
      if(sum(M_temp[i,]) > 0) top <- FALSE;
      
      #after we've found out whether it is a herbivore or a top
      #predator, or none of those, assign the value accordingly
      if(herb){
        TrLevels[i] = 1;
      }else if(top){
        TrLevels[i] = 3;
      }else{
        TrLevels[i] = 2;
      }
    }
  }
  FW$TrLevels <- TrLevels;
  return(FW);
}

##### Body sizes of species in the network following the method used in:

##### The susceptibility of species to extinctions in model communities. Binzer et al. (2011)
##### Basic and Applied Ecology

##### Simple prediction of interaction strengths in complex food webs. Berlow et al. (2009). PNAS.

ObtainBodySize <- function(FW){
  S <- dim(FW$M)[1];
  FW <- TrophicPositions(FW);
  FW$original_TP <- FW$TP;
  BS <- numeric(S);
  
  #the body sizes are assigned starting with the species with the lowest trophic position
  tp_order <- order(FW$TP)
  #average predator-prey mass ratios
  FW$ppmr <- 100; #rlnorm(1, meanlog=log(10), sdlog=log(10));
  
  #### we changed here the scaling of the body mass based on Binzer et al. 2015 Global Change Biology
  
  for(i in tp_order){
    if(FW$TP[i] == 1){  #body size of basal species == 1
      BS[i] <- 0.01; #*(FW$ppmr^(rnorm(1)));
    }else{
      BS[i] <- 0.01*(FW$ppmr^( ((FW$TP[i]-1)))); # + (rnorm(1, sd=0.1))  ))); #+ (rnorm(1)) ) ) );
    }
  }
  
  FW$BS <- BS;
  return(FW);
  
}

#This function obtains the parameters for the allometric predator-prey model as explained in:
# The susceptibility of species to extinctions in model communities. (2011). Binzer et al.
#   Basic and applied ecology

ObtainATNParams <- function(FW){
  S <- dim(FW$M)[1]
  ##### here we calculate the interaction matrix of alphas and handling times
  FW$alphas <- FW$M
  FW$Th <- FW$M
  for(i in 1:S){
    for(j in 1:S){
      if(FW$M[i,j] == 1){   #i is eaten by j
        FW$alphas[i,j] <- mp['d','alpha']*FW$BS[i]^mp['b','alpha']*FW$BS[j]^mp['c','alpha']*exp(mp['E','alpha']*temp_exp);
        FW$Th[i,j] <- mp['d','Th']*FW$BS[i]^mp['b','Th']*FW$BS[j]^mp['c','Th']*exp(mp['E','Th']*temp_exp);
      }
    }
  }
  
  #### this yields the weighted interaction matrix by the handling times
  FW$W <- FW$alphas * FW$Th
  
  ## Hill exponent
  FW$h <- 1.2
  
  FW$gr <- rep(0,S); # maximum growth rate (only for basal/resource species)
  FW$gr[FW$TP == 1] <- mp['d','r']*(FW$BS[FW$TP == 1]^mp['b','r'])*exp((mp['E','r']*temp_exp));
  
  FW$x <- rep(0,S); # death or metabolic or respiration rate (only for consumer species)
  FW$x[FW$TP > 1] <- mp['d','x']*(FW$BS[FW$TP > 1]^mp['b','x'])*exp((mp['E','x']*temp_exp));
  
  #assimilation efficiencies: e_ij = 0.85 for carnivores and e_ij = 0.45 for herbivores
  #### in the new version of the model (after Binzer et al. 2015) all e's = .85
  FW$e <- 0.85;  #this is a carnivorous link
  
  # we assign the carrying capacities for the basal species
  K <- numeric(S)
  K[FW$TP == 1] <- mp['d','K']*(FW$BS[FW$TP == 1]^mp['b','K'])*exp((mp['E','K']*temp_exp))
  FW$K <- K;
  return(FW);
  
}

BuildATN <- function(S, C){
  #create the niche model web
  network <- NicheNetwork(S,C);
  #obtain the body sizes for the species in the network
  network <- ObtainBodySize(network);
  #obtain the parameters for the allometric model
  network <- ObtainATNParams(network);
}



#### This is the function we give to the solver to solve the dynamics through time
#### It implements the equations that govern the biomass changes of species through time

SolveEcoDynamics <- function(t, Xt, FW){
  S <- dim(FW$M)[1]
  Xt[Xt <= 1e-9] <- 0;
  B <- Xt;
  
  #the functional response F_ji specifies the realised fraction
  #of consumption of species i on species j (i - predator, j - prey)
  F <- matrix(0, S, S);
  for(i in 1:S){
    if(FW$TP[i] > 1){
      fraction <- 1 + FW$W[,i] %*% (B^FW$h);
      for(j in 1:S){
        #### we do this to void the solver to fail. It nonetheless makes sense because the only way for fractions
        #### to be 0 is that all the other prey are 0 and then shouldn't affect its interaction with the current prey
        if(fraction == 0){
          print('in fraction')
          F[j,i] <- 0;  
        }else{
          interaction_j <- FW$alphas[j,i] * ((B[j])^(FW$h));
          F[j,i] <- interaction_j / fraction;
        }
      }
    }
  }
  
  
  # once we calculate the denominator of the functional response and the alphas
  # we can then calculate the changes of biomass for each species
  dB <- numeric(S);
  for(i in 1:S){
    tmp <- F[i,] * B;
    if(FW$original_TP[i] == 1.0 & FW$K[i] > 0.0){   # resources
      dB[i] <- FW$gr[i] * (1 - (B[i]/FW$K[i])) * B[i] - sum(tmp[tmp!=0]);
    }else{  # consumers
      dB[i] <- - FW$x[i]*B[i] + sum((F[,i]*B[i]*FW$e)) - sum(tmp[tmp!=0]);
    }
  }
  return(list(dB));
}


##### This function recalculates different aspects of the food web
##### Usually called after some changes have occurred like for example
##### the execution of ecological dynamics
##### If there have been extinctions the corresponding pointers
##### to the extinct species are removed
RecalculateFW <- function(FW, ext_threshold=1e-9){
  FW$collapsed <- FALSE;
  X <- FW$Final;
  X[X <= ext_threshold] <- 0;
  S <- length(X);
  
  #we first look at the species that went extinct
  index <- (X == 0.0);
  index <- (1:S)[index];
  FW$extinct_species$indexes <- index;
  FW$extinct_species$BSs <- FW$BS[index];
  FW$extinct_species$TPs <- FW$TP[index];
  FW$extinct_species$TLs <- FW$TrLevels[index];
  FW$extinct_species$niche <- FW$niche[index];
  
  if(sum(X) == 0){
    print('RecalculateFW :: the network collapsed - 0 species left');
    FW$S <- 0;
    FW$collapsed <- TRUE;
    return(FW);
  }
  
  index <- (X>0);
  index <- (1:S)[index];
  
  if(length(index) > 1){
    
    M_temp <- FW$M[index, index];
    
    #we need to check if nodes become isolated after extinctions have occurred
    #here we remove the nodes that are not connected anymore
    diag(M_temp) <- 0;
    remove <- c();
    for(i in 1:length(index)){
      if(sum(M_temp[i,]) == 0 & sum(M_temp[,i]) == 0 ){
        remove <- append(remove, index[i]);
        #this species has become disconnected and we have decided here to treat it as an extinction
        FW$extinct_species$indexes <- append(FW$extinct_species$indexes, index[i]);
        FW$extinct_species$BSs <- append(FW$extinct_species$BSs, FW$BS[index[i]]);
        FW$extinct_species$TPs <- append(FW$extinct_species$TPs, FW$TP[index[i]]);
        FW$extinct_species$TLs <- append(FW$extinct_species$TLs, FW$TrLevels[index[i]]);
        FW$extinct_species$niche <- append(FW$extinct_species$niche, FW$niche[index[i]]);
      }
    }
    
    index <- setdiff(index, remove);
  }
  
  
  SFinal <- length(index);
  
  if(SFinal == 0){
    print('RecalculateFW :: the network collapsed - 0 species left');
    FW$S <- 0;
    FW$collapsed <- TRUE;
    return(FW);
  }
  
  FW$S <- SFinal;
  FW$M <- FW$M[index, index];
  FW$C <- sum(FW$M)/FW$S**2;
  FW$radius <- FW$radius[index];
  FW$niche <- FW$niche[index];
  FW$centre <- FW$centre[index];
  
  #FW <- TrophicPositions(FW);   #within this function normalised matrix and trophic levels are recalculated as well
  FW$BS <- FW$BS[index];
  FW$gr <- FW$gr[index];
  FW$K <- FW$K[index];
  FW$Final <- FW$Final[index];
  
  return(FW);
}

RecalculateX <- function(out, ext_threshold=1e-9){
  m <- dim(out)[2];
  for(i in 1:(m-1)){
    out[out[,(i+1)]<=ext_threshold, (i+1)] <- 0;
  }
  return(out);
}

RunDynamics <- function(MaxTime=100, length, FW){
  times <- seq(0, MaxTime, by=length);
  out <- ode(y=FW$Initial, times=times, func=SolveEcoDynamics, parms=FW);
  return(out);
}


GetRankLayout <- function(FW){
  nodes <- nrow(as.matrix(FW$TrLevels));
  layout <- matrix(0, nodes, 2);
  n_per_level <- table(FW$TrLevels);
  gap <- matrix(0,4,1);
  
  middle <- 6;
  
  next_right <- matrix(middle,5,1);
  next_left <- matrix(middle,5,1);
  for(i in 1:4){  #the number of levels in the ranking algorithm
    gap[i] <- (middle*2)/n_per_level[names(n_per_level)==i-1];
  }
  
  gap[3] <- gap[3]*2;
  
  count <- as.vector(table(FW$TrLevels));
  
  even <- FALSE;
  for(i in 1:nodes){
    if(FW$TrLevels[i] == 0){
      layout[i,2] <- -6;
      
      if(abs(middle - next_left[1]) > abs(middle - next_right[1])){
        
        next_right[1] <- next_right[1] + gap[1];
        layout[i,1] <- next_right[1];
        
      }else{
        
        layout[i,1] <- next_left[1];
        next_left[1] <- next_left[1] - gap[1];
      }
      
      if(count[1]%%2==0) layout[i,1] <- layout[i,1] - gap[1]*0.5;
      
    }else if(FW$TrLevels[i] == 1){
      layout[i,2] <- -3.5;
      if(abs(middle - next_left[2]) > abs(middle - next_right[2])){  
        next_right[2] <- next_right[2] + gap[2];
        layout[i,1] <- next_right[2];
        
      }else{
        
        layout[i,1] <- next_left[2];
        next_left[2] <- next_left[2] - gap[2];
      }
      
      if(count[2]%%2==0) layout[i,1] <- layout[i,1] - gap[2]*0.5;
    }else if(FW$TrLevels[i] == 2){
      if(even){
        layout[i,2] <- 1;
        if(abs(middle - next_left[3]) > abs(middle - next_right[3])){  
          next_right[3] <- next_right[3] + gap[3];
          layout[i,1] <- next_right[3];
          
        }else{          
          
          layout[i,1] <- next_left[3];
          next_left[3] <- next_left[3] - gap[3];
        }
      }else{
        layout[i,2] <- -0.5;
        if(abs(middle - next_left[4]) > abs(middle - next_right[4])){  
          next_right[4] <- next_right[4] + gap[3];
          layout[i,1] <- next_right[4];
          
        }else{          
          
          layout[i,1] <- next_left[4];
          next_left[4] <- next_left[4] - gap[3];
        }
      }
      even <- !even;
      
      if(count[3]%%2==0) layout[i,1] <- layout[i,1] - gap[3]*0.5;
      
    }else{
      layout[i,2] <- 5;
      if(abs(middle - next_left[5]) > abs(middle - next_right[5])){   
        next_right[5] <- next_right[5] + gap[4];
        layout[i,1] <- next_right[5];
      }else{
        layout[i,1] <- next_left[5];
        next_left[5] <- next_left[5] - gap[4];
      }
      if(count[4]%%2==0) layout[i,1] <- layout[i,1] - gap[4]*0.5;
      
    }
  }
  return(layout);
}
