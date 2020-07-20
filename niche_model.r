

library(igraph)
library(deSolve)

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
      print("NicheNetwork :: the network is not connected...");
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
          print("NicheNetwork :: the network has cycles with no external energy input...");
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
          print("NicheNetwork :: there is a duplicate species");
          break;
        }   
      }
      
      #as long as we are traversing the network we might as well check
      #for basal species with cannibalistic links...
      
      if(sum(preys_i) == 1 && M[i,i] == 1){
        print("NicheNetwork :: a cannibalistic link was found in a basal species... removing it");
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
  
  #here we can specify whether the invasive species should have at least 4 links
  #so we can perform the rewiring experiments afterwards....
  
  while(sum(M[,inv_idx]) < 1 || sum(M[inv_idx,]) < 1) inv_idx <- sample(1:S, 1);
  
  
  print(paste('the invader has', sum(M[,inv_idx]), 'prey and', sum(M[inv_idx,]), 'predators'));
  
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


#calculate the trophic positions using the prey-averaged trophic level metric from:
# Richar J. Williams and Neo D. Martinez (2004). American Naturalist.
# Limits to trophic levels and omnivory in complex food webs: theory and data.

TrophicPositions <- function(FW){
  
  S <- dim(FW$M)[1];
  
  M <- NormalizeMatrix(FW$M);
 
  #print('printing matrix');
  #print(FW$M);
  #print(M);
  #print(diag(S));
  
  #TP[j] = 1 + sum_over_i ( M[i,j] TL[i])
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


#we can also obtain the trophic level of each species
#e.g. basal, herbivore, omnivore, top predator
TrophicLevels <- function(FW){
  S <- dim(FW$M)[1];
  TrLevels <- rep(-1, S)
  
  #we find the neighbours of each node in the network
  M_temp <- FW$M;
  diag(M_temp) <- 0;
  
  #all the species with trophic position 1 are basal species, and therefore
  #belong to the trophic level 0
  TrLevels[FW$TP<=1.0] <- 0;
  
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

#we estimate the body sizes of species in the network following the method used in:

#The susceptibility of species to extinctions in model communities. Binzer et al. (2011)
#   Basic and Applied Ecology

#and Simple prediction of interaction strengths in complex food webs. Berlow et al. (2009). PNAS.

ObtainBodySize <- function(FW){
  S <- dim(FW$M)[1];
  FW <- TrophicPositions(FW);
  
  FW$original_TP <- FW$TP;
  
  BS <- numeric(S);
  
  #the body sizes are assigned starting with the species with the lowest trophic position
  tp_order <- order(FW$TP)
  
  #average predator-prey mass ratios
  FW$ppmr <- 10^2; #rlnorm(1, meanlog=log(10), sdlog=log(10));
  
  for(i in tp_order){
    if(FW$TP[i] <= 1){  #body size of basal species == 1
      BS[i] <- 1;
    }else{
      BS[i] <- FW$ppmr^(FW$TP[i]-1);
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
  nprods <- length(FW$TP[FW$TP<=1]);
  FW$B0 <- 0.5;
  
  FW$Wless <- FW$M/2;
  
  ######################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###################
  
  ###CHANGE THESE THINGS TO BE ONE PER SPECIES TO SEE HOW MUCH THIS AFFECTS THE
  ###BEHAVIOUR OF THE SYSTEM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #Hill exponent and predator interference term
  #pred_int <- -1;
  #while(pred_int < 0 | pred_int > 2){
  #  pred_int <- rnorm(1, mean=0.1, sd=0.05);
  #}
  #h <- 0;
  #while(h < 1 | h > 2){
  #  h <- rnorm(1, mean=1.1, sd=0.05);
  #}
  
  #for now, we are going to keep the parameters constant, so...
  pred_int <- 0.1;
  #h <- 1.1;
  
  #######################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###################
  
  FW$pred_int <- pred_int;
  #FW$h <- h;
  
  ##we give a particular hill exponent to each one of the species
  ##to tell how efficient they are capturing their prey
  ##the default value is 1.1
  FW$h <- rep(1.5, S); #maximum growth rate
  
  FW$y <- 8; #maximum consumption rate
  FW$gr <- rep(0,S); #maximum growth rate
  FW$gr[FW$TP<=1] <- 1;
  
  FW$Mp <- mean(FW$BS[FW$TP<=1.0])  #mean body mass producers
  #FW$Mc <- mean(FW$BS[FW$TP>1.0])   #mean body mass of consumers
  
  FW$metabolic <- numeric(S)  #metabolic rate
  #alpha_x / alpha_r = 0.314 for all normalized metabolic rates
  #this approach was wrong....
  #FW$metabolic[FW$TP<=1.0] <- 0.314 * (FW$Mp / FW$Mp) ^ -0.25;  #producers
  #FW$metabolic[FW$TP>1.0] <- 0.314 * (FW$Mc / FW$Mp) ^ -0.25;   #consumers
  
  FW$metabolic <- 0.314 * (FW$BS / FW$Mp) ^ -0.25;
  
  #assimilation efficiencies: e_ij = 0.85 for carnivores and e_ij = 0.45 for herbivores
  FW$e <- matrix(0, S, S);
  for(i in 1:S){
    for(j in 1:S){
      if(FW$M[i,j] == 1){   #i is eaten by j
        if(FW$TP[i] <= 1.0){  #this link is herbivorous
          FW$e[i,j] <- 0.45;
        }else{
          FW$e[i,j] <- 0.85;  #this is a carnivorous link
        }    
      }
    }
  }
  
  #we assign the carrying capacities for the basal species
  #K_sys = S
  K <- numeric(S)
  K[FW$TP <= 1] <- 5/nprods;
  FW$K <- K;
  return(FW);
  
}

#this function allows us to obtain a food-web like layout
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

#build the allometric trophic network based on the parameters above

BuildATN <- function(S, C){
  #create the niche model web
  network <- NicheNetwork(S,C);
  #obtain the body sizes for the species in the network
  network <- ObtainBodySize(network);
  #obtain the parameters for the allometric model
  network <- ObtainATNParams(network);
}


#implement the dynamics of the model

SolveEcoDynamics <- function(t, Xt, FW){
  
  #print(paste("solving dynamics...", t));
  
  S <- dim(FW$M)[1]
  Xt[Xt <= 1e-9] <- 0;
  Xt[Xt >= 1e9] <- 0;
  
 # Xt[is.nan(Xt)] <- 0;
  
  
  B <- Xt;
  
  #print(t);
  #print(B);
  
  #the functional response F_ji specifies the realised fraction
  #of consumption of species i on species j (i - predator, j - prey)
  
  F <- matrix(0, S, S);
  
  #B[(is.nan(B))] <- 0.0;
  #B[(is.na(B))] <- 0.0;
  
  #for(i in 1:length(B)){
  #  if(is.nan(B[i]) || is.na(B[i])){
  #    B[i] <- 0;
  #  }
  #}
  
  for(i in 1:S){
    if(FW$TP[i] > 1){
      
      #print(B)
      #print(FW$h)
      
      fraction <- (FW$B0)^(FW$h[i]) + FW$pred_int*B[i] + FW$W[,i] %*% (B^FW$h[i]);
      
      #print(paste("B0", FW$B0, "h", FW$h, "pred int", FW$pred_int));
      #print(B^FW$h);
      #print(FW$W[,i]);
      
      
      #print(paste("fraction" ,fraction));
      
      for(j in 1:S){
        interaction_j <- FW$W[j,i] * ((B[j])^(FW$h[i]));
        if(is.nan(fraction) || is.na(fraction)){
          
          
          
          print(paste("fraction" ,fraction));
          print(paste("B0", FW$B0, "h", FW$h, "pred int", FW$pred_int));
          print(t);
          
          print(FW$Final);
          
          print(Xt);
          
          print(B);
          print(B^FW$h);
          print(FW$W[,i]);
          
          print(FW$metabolic);
          #print(FW$e);
          
          print(FW$Mc);
          print(FW$Mp);
          print(FW$ppmr);
          print(FW$ratios);
          
          if(fraction == 0){
            F[j,i] <- 0;  
          }
          
        }else{
          F[j,i] <- interaction_j / fraction;
        }
        
        #print(F[j,i]);
        
      }
      
    }
  }
  
  #then, we calculate the changes of biomass for each species
  dB <- numeric(S);
  for(i in 1:S){
    tmp <- F[i,] * FW$metabolic * B * FW$y;
    if(FW$original_TP[i] <= 1.0 & FW$K[i] > 0.0){   #producers...
      dB[i] <- FW$gr[i] * (1 - B[i]/FW$K[i]) * B[i] - sum(tmp[tmp!=0]/FW$e[i,tmp!=0]);
      
      #dB[i] <- FW$gr[i] * (1 - B[i]/FW$K[i]) * B[i] - sum(tmp[tmp!=0]);
      
      #if(is.nan(dB[i]) || is.na(dB[i])){
        #print(paste(i,"dB_i is nan or na..."));
        #print(tmp);
        #print(paste("tmp", tmp, "gr", FW$gr[i], "B_i", B[i], "K", FW$K[i], sum(tmp[tmp!=0]/FW$e[i,tmp!=0])));
        #print(FW);
        
        #if this happens an animal is now on the basal level and is in the process of dying
        #hence its change of rate has to be negative
      #  dB[i] <- -1;
     # }
      
    }else{  #consumers
      dB[i] <- -FW$metabolic[i]*B[i] + sum(F[,i]*B[i]*FW$metabolic[i]*FW$y) - sum(tmp[tmp!=0]/FW$e[i,tmp!=0]);
      
      #dB[i] <- -FW$metabolic[i]*B[i] + sum(F[,i]*B[i]*FW$metabolic[i]*FW$y*FW$e[,i]) - sum(tmp[tmp!=0]);
      
      
     # if(is.nan(dB[i]) || is.na(dB[i])){
        #print("dB_i is nan or na... in predator");
        #print(tmp);
        #print(paste("tmp", tmp, "F_i", F[,i], "B_i", B[i], "metabolic ", FW$metabolic[i], "y", FW$y, sum(tmp[tmp!=0]/FW$e[i,tmp!=0])));
         
     #   dB[i] <- -1;
     # }
      
      
    }
  }

  #print(paste("dB ", list(dB)));
  
  return(list(dB));
}


#this is to recalculate the biomasses of the species and the network features
#once the updates have been done
RecalculateFW <- function(FW){
  FW$collapsed <- FALSE;
  
  X <- FW$Final;
  X[X <= 1e-9] <- 0;
  X[X >= 1e9] <- 0;
  
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
  
  M_temp <- FW$M[index, index];
  #we need to check if nodes become isolated after extinctions have occurred
  #here we remove the nodes that are not connected anymore
  diag(M_temp) <- 0;
  
  remove <- c();
  for(i in 1:length(index)){
    if(sum(M_temp[i,]) == 0 & sum(M_temp[,i]) == 0 ){
      remove <- append(remove, index[i]);
      #this species has become disconnected and we have decided here to treat it as an
      #extinct species.
      FW$extinct_species$indexes <- append(FW$extinct_species$indexes, index[i]);
      FW$extinct_species$BSs <- append(FW$extinct_species$BSs, FW$BS[index[i]]);
      FW$extinct_species$TPs <- append(FW$extinct_species$TPs, FW$TP[index[i]]);
      FW$extinct_species$TLs <- append(FW$extinct_species$TLs, FW$TrLevels[index[i]]);
      FW$extinct_species$niche <- append(FW$extinct_species$niche, FW$niche[index[i]]);
      
    }
  }
  
  index <- setdiff(index, remove);
  
  SFinal <- length(index);
  FW$S <- SFinal;
  FW$M <- FW$M[index, index];
  FW$C <- sum(FW$M)/FW$S**2;
  FW$radius <- FW$radius[index];
  FW$niche <- FW$niche[index];
  FW$centre <- FW$centre[index];
  
  
  #print('this is the call within Recalculate FW');
  
  #print(paste('S',FW$S, 'C', FW$C, 'n', FW$niche, 'r', FW$radius, 'index', index));
  FW <- TrophicPositions(FW);   #within this function normalised matrix and trophic levels are recalculated as well
  FW$BS <- FW$BS[index];
  
  FW$gr <- FW$gr[index];
  FW$h <- FW$h[index];
  
  #FW$Mp <- mean(FW$BS[FW$TP <= 1.0]);
  #FW$Mc <- mean(FW$BS[FW$TP >1.0]);
  
  #FW$metabolic <- numeric(SFinal)  #metabolic rate
  #alpha_x / alpha_r = 0.314 for all normalized metabolic rates
  #FW$metabolic[FW$TP<=1.0] <- 0.314 * (FW$Mp / FW$Mp) ^ -0.25;  #producers
  #FW$metabolic[FW$TP>1.0] <- 0.314 * (FW$Mc / FW$Mp) ^ -0.25;   #consumers
  
  #when we update the network after possible extinctions we shouldn't recalculate
  #the metabolic rates, since these should be maintained by the species,
  #specially if they are not changing body size
  
  FW$metabolic <- FW$metabolic[index];
  
  FW$e <- FW$e[index, index];
  FW$K <- FW$K[index];
  
  #FW$Initial <- FW$Initial[index];
  FW$Final <- FW$Final[index];
  
  return(FW);
  
}

RecalculateX <- function(out){
  m <- dim(out)[2];
  for(i in 1:(m-1)){
    out[out[,(i+1)]<=1e-9, (i+1)] <- 0;
    out[out[,(i+1)]>=1e9, (i+1)] <- 0;
    
    #out[is.nan(out[,(i+1)]), (i+1)] <- 0;
    
  }
  return(out);
}

CheckStability <- function(FW, Tmax, threshold=1e-2){
  tlength <- 100;
  t <- 0;
  delta <- 1;
  
  while(delta > threshold && t < Tmax){
    times <- seq(0, 1000, length=tlength);
    
    #print(paste('inside check stability', 'S', FW$S, 'C', FW$C, 'final', FW$Final));
    
    #print('FW$Final = ');
    #print(FW$Final);
    
    #print(times);
    
    out1 <- ode(y=FW$Final, times=times, func=SolveEcoDynamics, parms=FW);
    out2 <- ode(y=as.numeric(out1[tlength, -1]), times=times, func=SolveEcoDynamics, parms=FW);
    
    out1 <- RecalculateX(out1);
    out2 <- RecalculateX(out2);
    
    avg1 <- apply(out1[,-1], 2, mean);
    avg2 <- apply(out2[,-1], 2, mean);
    
    #print(out1[tlength,]);
    #print(out2[1,]);
    
    delta <- max(abs(avg1-avg2)/apply(rbind(avg1, avg2), 2, max))
    
    #print(delta);
    
    FW$Final <- out2[tlength, -1];
    FW <- RecalculateFW(FW);
    
    if(FW$S == 0){
      return(FW);
    }
    
    if(delta < threshold){
      FW$State <- 1;
      break;
    }else{    #not at equilibrium or limit cycle
      FW$State <- 0;
    }
    
    t <- t+1000;
    
  }
  
  FW$StabilityThreshold <- threshold;
  FW$TimeToStability <- t;
  
  return(FW);
}


GetStabilityMeasures <- function(FW, TMax=4000, tlength=200){
  
  if(! is.null(FW$dynamics)){
    
    #steps <- dim(FW$dynamics)[1];
    #init <- as.integer(steps/2);
    
    #final_dynamics <- FW$dynamics[init:steps,];
    
    final_dynamics <- FW$dynamics;
    tlength <- dim(final_dynamics)[1];
    
    
      
  }else{
    #each_t <- as.integer(TMax/2);
    
    times <- seq(0, TMax, length=tlength);
    #initial_dynamics <- RecalculateX(ode(y=FW$Final, times=times, func=SolveEcoDynamics, parms=FW));
    
    #FW$Final <- initial_dynamics[tlength, -1];
    #FW <- RecalculateFW(FW);
    
    final_dynamics <- RecalculateX(ode(y=FW$Final, times= times, func=SolveEcoDynamics, parms=FW));
    
    FW$Final <- final_dynamics[tlength, -1];
  }
  
  
  FW$overall_stability <- FALSE;
  final <- tlength;
  
  tlength <- 200;
  start <- 1;
  end <- start + tlength;
  
  half_interval <- as.integer(tlength/2);
  
  FW$stability_a$State <- 0;
  FW$stability_b$State <- 0;
  FW$stability_c$State <- 0;
  
  
  while(!FW$overall_stability & start < final){
    
    current_dynamics <- final_dynamics[start:end,];
    
    if(FW$stability_a$State == 0){
      sum_abundances <- apply(current_dynamics[,-1], 1, sum);
      sd_sum <- sd(sum_abundances);
      
      FW$stability_a$sums <- sum_abundances;
      FW$stability_a$sds <- sd_sum;
      if(sd_sum < 0.1){
        FW$stability_a$State <- 1;
        FW$stability_a$time_to_stable <- end;
      }#else{
        #FW$stability_a$State <- 0;
      #}
    }
    
    if(FW$stability_b$State == 0){
      cv_abs_time <- apply(current_dynamics[,-1], 2, sd)/apply(current_dynamics[,-1], 2, mean);
      FW$stability_b$cvs <- cv_abs_time;
      
      FW$stability_b$cv_mean <- mean(cv_abs_time);
      if(FW$stability_b$cv_mean < 0.1){
        FW$stability_b$State <- 1;
        FW$stability_b$time_to_stable <- end;
      }#else{
      #  FW$stability_b$State <- 0;
      #}
    
    }
    
    
    if(FW$stability_c$State == 0){
      #half_interval <- as.integer(tlength/2);
      first_half <- current_dynamics[1:half_interval,-1];
      second_half <- current_dynamics[(half_interval+1):(dim(current_dynamics)[1]),-1];
      
      avg1 <- apply(first_half, 2, mean);
      avg2 <- apply(second_half, 2, mean);
      
      delta <- max(abs(avg1-avg2)/apply(rbind(avg1, avg2), 2, max));
      
      FW$stability_c$delta <- delta;
      if(FW$stability_c$delta < 0.1){
        FW$stability_c$State <- 1;
        FW$stability_c$time_to_stable <- end;
      }#else{
      #  FW$stability_c$State <- 0;
      #}
    }
    
    if(FW$stability_a$State == 1 & FW$stability_b$State == 1 & FW$stability_c$State == 1){
      FW$overall_stability <- TRUE;
    }
    
    start <- end;
    end <- start + tlength;
    
  }
  return(FW);
  
}

GetStabilityDeltas <- function(FW, TMax){
  tlength <- 100;
  t <- 0;
  
  results <- list(times=c(), deltas=c());
  
  while(t < Tmax){
    times <- seq(0, 1000, length=tlength);
    out1 <- RecalculateX(ode(y=FW$Final, times=times, func=SolveEcoDynamics, parms=FW));
    out2 <- RecalculateX(ode(y=as.numeric(out1[tlength, -1]), times=times, func=SolveEcoDynamics, parms=FW));
    
    avg1 <- apply(out1[,-1], 2, mean);
    avg2 <- apply(out2[,-1], 2, mean);
    
    delta <- max(abs(avg1-avg2)/apply(rbind(avg1, avg2), 2, max))
    
    print(delta);
    
    FW$Final <- out2[tlength, -1];
    FW <- RecalculateFW(FW);
    
    t <- t+1000;
    
    results$times <- append(results$times, t);
    results$deltas <- append(results$deltas, delta);
      
    
  }
  
  return(results);
  
}

GetStableFoodWeb <- function(S=100, C=0.08, fraction=1.0, stability=TRUE){
  fwno <- sample(1:1000, 1);
  Tmax <- 1000;
  tlength <- 200;
  
  set.seed(fwno);
  finished <- FALSE;
  
  trials_sps <- 0;
  trials_stab <- 0;
  trials_not_c <- 0;
  
  while(!finished){
    FW <- BuildATN(S,C); 
    FW$orig_fw <- FW;
    
    Xt <- runif(S, 0.05, 1);
    FW$Initial <- Xt;
    FW$Final <- Xt;
    
    FW <- ObtainDynamics(FW);
    
    #exec <- RecalculateX(RunDynamics(2000, tlength, FW=FW));
    
    #FW$Final <- exec[tlength, -1];
    FW <- RecalculateFW(FW);
    
    M_temp <- FW$M;
    diag(M_temp) <- 0;
    graf <- graph.adjacency(M_temp);
    comps <- igraph::no.clusters(graf, mode='weak');
    if(comps > 1){
      trials_not_c <- trials_not_c + 1;
      print("GetStableFoodWeb :: the network is not connected...");
      next;
    }
    
    if(FW$S < as.integer(S*fraction)){
      trials_sps <- trials_sps + 1;
      print("GetStableFoodWeb :: too few species, starting over...");
      next;
    }
    
    #look at the dynamics to find out whether the network is stable
    FW <- GetStabilityMeasures(FW, Tmax, tlength);
    
    if(stability & FW$stability_a$State == 0){
      trials_stab <- trials_stab + 1;
      print('GetStableFoodWeb :: Model not stable, looking for a stable one...');
      next;
    }
    
    if(!stability & FW$stability_a$State == 1){
      trials_stab <- trials_stab + 1;
      print('GetStableFoodWeb :: Model stable, but looking for an unstable one...');
      next;
    }
    
    finished <- TRUE;
  }
  
  
  FW$S_lost_b4inv <- FW$extinct_species;
  FW$S_lost_b4inv$number <- FW$orig_fw$S - FW$S;
  
  #here we include inside the food web values a variable keeping record of
  #how many iterations where needed to find a network with the appropriate
  #number of species
  FW$unsucc_trials_sps <- trials_sps;
  
  #and within those that had the right number of species, how many of those
  #were unstable and hence needed further trials
  FW$unsucc_trials_stab <- trials_stab;
  
  FW$unsucc_trials_nc <- trials_not_c;
  
  return(FW); 
}

RunDynamics <- function(MaxTime=100, length, FW){
  times <- seq(0, MaxTime, length=length);
  out <- ode(y=FW$Initial, times=times, func=SolveEcoDynamics, parms=FW);
  return(out);
  
}

ObtainDynamics <- function(FW){
  Tmax <- 2000;
  tlength <- 2001;
  times <- seq(0, Tmax, length=tlength);
  out <- ode(y=FW$Initial, times=times, func=SolveEcoDynamics, parms=FW);
  
  FW$dynamics <- out;
  
  FW$Final <- FW$dynamics[tlength, -1];

  return(FW);
  
  
}

#write.table(a$M, file = "network.csv", append=FALSE, sep=",", col.names=FALSE, row.names=FALSE)
#analyse.single(filename="network.csv")
#plotweb(col=1:9, radii=c(7,8,9,10,11,12,13,14,15))

RunModel <- function(S=100, C=0.08){
  fwno <- sample(1:1000, 1);
  
  #S <- 100;
  #C <- 0.08;
  
  Tmax <- 4000;
  tlength <- 200;
  
  set.seed(fwno);
  
  
  finished <- FALSE;
  
  while(!finished){
    FW <- BuildATN(S,C);
    FW$orig_fw <- FW;
    
    Xt <- runif(S, 0.05, 1);
    FW$Initial <- Xt;
    FW$Final <- Xt;
    
    exec <- RecalculateX(RunDynamics(Tmax, tlength, FW=FW));
    
#    print(sum(exec[tlength, -1]));
    
    FW$Final <- exec[tlength, -1];
    FW <- RecalculateFW(FW);
    
    #print("The execution was ok");
    #print(paste("number of species = ", FW$S));
    
    #if(FW$S < S/2){
    #  print("too few species, starting over...");
    #  next;
    #}
    
    #print(FW$Final);
    #print(FW$W);
    
    #integrate the dynamics
#    FW <- CheckStability(FW, Tmax);
#    FW <- RecalculateFW(FW);
    
    #print(FW$S);
    
#    if(FW$State == 0){
#      print("Model not stable");
#      next;
#    }
    
    finished = TRUE;
    
    extinct <- S - FW$S; 
    
  }
  
  # save(FW, file=paste('RFile-FW-S-',FW$S,'-C-', round(FW$C, 3), '-no-', fwno, '.RData', sep=''))
  
  print("RunModel :: plot the dynamics");
  # pdf(paste('Plot-dynamics-S',FW$S,'-C-', round(FW$C, 3), '-no-', fwno, '.pdf', sep=''), width=10, height=4);
  par(mar=c(4,4,1,1));
  plot(0,0,type='n', xlab='Time', ylab='Species biomass', xlim=range(exec[,1]), ylim=range(exec[,-1]))
  
  colours <- rainbow(FW$S);
  for(k in 1:FW$S){
    lines(exec[,1], exec[,1+k], col=colours[k], lwd=2);
  }
  
  legend("topleft", paste(c("Initial = ", "Extinctions = "), c(S, extinct)), bty='n');
  # dev.off();
}


#a <- BuildATN(100, 0.08);
#a <- TrophicPositions(a)
#sna::gplot(a$M, diag=TRUE, coord=GetRankLayout(a))

#for(i in 1:10000){
# print(paste("running", i));
 
# S <- as.integer(msm::rtnorm(1, mean=100, sd=20, lower=80, upper=120));
# C <- msm::rtnorm(1, mean=0.08, sd=0.01, lower=0.06, upper=0.1);
 
# RunModel(S,C);
#}


################ This was a test for the dingo food web!!!! #######################

BuildATNModule <- function(){
  #create the niche model web
  
  M <- matrix(c(0,1,1,1,0,0,0, 0,0,0,0,1,1,1, 0,0,0,0,1,1,0, 0,0,0,0,0,0,1, 0,0,0,0,0,0,1, 0,0,0,0,0,0,1, 0,0,0,0,0,0,0), ncol=7, nrow = 7, byrow=T)
  
  network <- list(S=dim(M)[1], C=sum(M)/((dim(M)[1])**2), M=M, original_C=sum(M)/((dim(M)[1])**2))
  #obtain the body sizes for the species in the network
  network <- ObtainBodySize(network);
  #obtain the parameters for the allometric model
  network <- ObtainATNParams(network);
}


RunModule <- function(){
  fwno <- sample(1:1000, 1);
  set.seed(fwno);
  Tmax <- 1000;
  tlength <- 200;
  finished <- FALSE;
  
  while(!finished){
    FW <- BuildATNModule();
    FW$orig_fw <- FW;
    
    Xt <- runif(7, 0.05, 1);
    FW$Initial <- Xt;
    FW$Final <- Xt;
    
    exec <- RecalculateX(RunDynamics(Tmax, tlength, FW=FW));
    FW$Final <- exec[tlength, -1];
    FW <- RecalculateFW(FW);
    finished = TRUE;
    extinct <- 7 - FW$S;
  }
  
  FW <- TrophicPositions(FW)
  sna::gplot(FW$M, diag=TRUE, coord=GetRankLayout(FW))
  
  # save(FW, file=paste('RFile-FW-S-',FW$S,'-C-', round(FW$C, 3), '-no-', fwno, '.RData', sep=''))
  
  print("RunModel :: plot the dynamics");
  # pdf(paste('Plot-dynamics-S',FW$S,'-C-', round(FW$C, 3), '-no-', fwno, '.pdf', sep=''), width=10, height=4);
  par(mar=c(4,4,1,1));
  plot(0,0,type='n', xlab='Time', ylab='Species biomass', xlim=range(exec[,1]), ylim=range(exec[,-1]))
  
  colours <- c('black', 'red', 'blue', 'green', 'yellow', 'magenta', 'cyan') #rainbow(FW$S);
  for(k in 1:FW$S){
    lines(exec[,1], exec[,1+k], col=colours[k], lwd=2);
  }
  
  legend("topleft", paste(c("Initial = ", "Extinctions = "), c(7, extinct)), bty='n');
  # dev.off();
}

RunModule()
                           
                           