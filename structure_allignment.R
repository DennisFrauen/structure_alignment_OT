library(Rpdb)
library(wordspace)
library(rgl)
library(pracma)
library(viridis)
library(hrbrthemes)
library(plotly)
source("distances.R")
source("monotone_rearrangement.R")
source("structure_align_plotting.R")

#Calculates structure alignment using monotone rearrangement
#Input: p1, p2 Protein 3d coordinates (Nx3 matrix, Mx3 matrix)
#       kappa: Number of aligned residuals (if 0, kappa is choosen by 
#                                     optimising an alignement measure)
#       dist_cutoff: distance cutoff value (0 means no cutoff)
#       cutoff_by_dist: If TRUE: distances larger than dist_cutoff will be cut off when 
#                       calculating the inter-residue distributions
#                       If FALSE: for each residue, only take the 1...dist_cutoff smallest
#                                  distances when calculating the inter-res distances
#       seq_cutoff: for each residue, only consider distances to residues which are within
#                    a neighborhood of size seq_cutoff 
#       iter: Number of iterations
#       alpha, p: parameters of partial Wasserstein distance
#       dist: Distance of distributions chosen to compare inter-residue distributions
#       wstein_matr: Matrix of Wasserstein distances, may speed up computation if provided
#       init_align: Initial alignment. By default monotone rearrangement is applied
#       remove_res: Indicates whether inter-residue distances of unaligned distances
#                   are removed from inter-residue distribtions after every interation
#       return_meas: Indicates whether alignment measures over all iterations are returned
#       print_info: Indicates whether progress information is printed
align_proteins <- function(p1,p2,kappa = 0,dist_cutoff = 0,seq_cutoff=0,
                      iter=1,lambda=0,alpha=1,p=1,dist = "wasserstein", wstein_matr = c(),
                      init_align = c(),remove_res=FALSE,return_meas=FALSE,print_info=TRUE){
  
  #First use all inter-residue distances
  remove_res1 <- c()
  remove_res2 <- c()
  rmsd <- 0
  #Calculate Loss matrix (e.g. partial Wasserstein matrix) if not provided
  if(length(wstein_matr)==0){
    if(print_info==TRUE){
      print("Compute IRD")
    }
    #Get inter-residue distances
    ird1 <- get_inter_res_dist(p1,dist_cutoff,seq_cutoff,c())
    ird2 <- get_inter_res_dist(p2,dist_cutoff,seq_cutoff,c())
    if(print_info==TRUE){
      print("Compute cost matrix")
    }
    #Create (Wasserstein) loss matrix
    wstein_matr <- wstein_dist_matr(ird1,ird2,p = p,alpha = alpha,dist=dist)
  }
  
  #Empty vectors to save alignment measures over iterations
  nali_iter <- c()
  rmsd_iter <- c()
  so_iter <- c()
  
  #Check for initial alignment
  if(length(init_align)>0){
    print("Take initial alignment")
    alignment <- init_align
    start_iter <- 2
    if(remove_res==TRUE){
      remove_res1 <- setdiff(1:dim(p1)[1],alignment$I)
      remove_res2 <- setdiff(1:dim(p2)[1],alignment$J)
    }
    #Save alignment measures
    align_meas <- get_align_measures(p1,p2,alignment)
    nali_iter <- c(nali_iter,align_meas$nali)
    rmsd_iter <- c(rmsd_iter,align_meas$rmsd)
    so_iter <- c(so_iter,align_meas$so)
  } else{
    start_iter <- 1
  }
  
  #Start iterative procedure
  for(i in start_iter:iter){
    #From the second iteration step on update the Loss matrix by adding the lambda penalty
    if(i > 1){
      #Compute pairwise eucledian dist after superposition
      eucl_dist_matr <- eukl_dist_matrix(p1,p2,list(I=alignment$I,J=alignment$J))
      #lambda = 1000 means that only euclidean distances are considered
      if(lambda < 1000){
        #In case residues are removed, recalculate the loss matrix
        if(remove_res==TRUE){
          ird1 <- get_inter_res_dist(p1,dist_cutoff,seq_cutoff,remove_res1)
          ird2 <- get_inter_res_dist(p2,dist_cutoff,seq_cutoff,remove_res2)
          print("Compute Wasserstein")
          wstein_matr <- wstein_dist_matr(ird1,ird2,p = p,alpha = alpha,dist=dist)
        }
        #Add euclidean distances weighted with lambda
        D <- wstein_matr + lambda * eucl_dist_matr
      } else{
        D <- eucl_dist_matr
      }
    } else{
      #In the first iteration, only consider the Loss matrix
      D <- wstein_matr
    }
    #Calculate H tensor (dynamic programming/ monotone rearrangement)
    if(print_info==TRUE){
      print("Compute H")
    }
    H <- monotone_rearr(D,kappa)
    #Recover alignment (traceback procedure)
    if(print_info==TRUE){
      print("Recover alignment")
    }
    #Check whether number of aligned residues is provided
    if(kappa == 0){
      #Get the alignment for every kappa
      max_kappa <- min(dim(p1)[1],dim(p2)[1])
      so_vec <- vector(length = max_kappa-1)
      for(k in 3:(max_kappa+1)){
        algn <- get_optimal_align(H[,,1:k])
        so_vec[k-2] <- get_align_measures(p1,p2,algn,"so")
      }
      #get optimal kappa (in terms of so) and alignment, in case of ties take maximal kappa
      kappa_opt <- max(which(so_vec == max(so_vec)))+1
      alignment <- get_optimal_align(H[,,1:(kappa_opt+1)]) 
    } else{
      #Get the alignment for fixed kappa
      alignment <- get_optimal_align(H)
    }
    
    #For next iterations set the unmatched residuals to be removed
    if(remove_res==TRUE){
      remove_res1 <- setdiff(1:dim(p1)[1],alignment$I)
      remove_res2 <- setdiff(1:dim(p2)[1],alignment$J)
    }
    
    #Save alignment measures for the current iteration step
    align_meas <- get_align_measures(p1,p2,alignment)
    nali_iter <- c(nali_iter,align_meas$nali)
    rmsd_iter <- c(rmsd_iter,align_meas$rmsd)
    so_iter <- c(so_iter,align_meas$so)
    #Print alignment measures of current iteration step
    if(print_info==TRUE){
      print(paste("Iteration",i))
      print(align_meas$so)
      print(align_meas$rmsd)
      print(align_meas$nali)
    }
    #Check whether SO was reduced -> stop iterations
    so_reduced <- FALSE
    if(i > 1){
      if(align_meas$so < so_iter[length(so_iter)-1]){
        print("SO reduced")
        so_reduced <- TRUE
        nali_iter <- nali_iter[1:length(nali_iter)-1]
        rmsd_iter <- rmsd_iter[1:length(rmsd_iter)-1]
        so_iter <- so_iter[1:length(so_iter)-1]
      }
    }
    #If the rmsd does not change, stop iterations
    if(align_meas$rmsd==rmsd | so_reduced == TRUE){
      if(return_meas==FALSE){
        return(alignment)
      } else{
        return(list(alignment=alignment,measures=list(so = so_iter,rmsd = rmsd_iter,nali=nali_iter)))
      }
    }
    rmsd <- align_meas$rmsd
  }
  #return alignment
  if(return_meas==FALSE){
    return(alignment)
  } else{
    return(list(alignment=alignment,measures=list(so = so_iter,rmsd = rmsd_iter,nali=nali_iter)))
  }
}


#Calculation of alignment measures (superimposition + nali, rmsd, so, sp-score)
#Input: pr1, pr2 Protein 3d coordinates (Nx3 matrix, Mx3 matrix)
#       alignment: Structural alignment of the proteins (list of two index sets I and J) 
#       type: Type of alignment measure to be returned ("nali","rmsd","so","sp_score","all")
#       superpos: boolean to indicate whether protein 1 is translated + rotated by kabsch
#                 algorithm before calculating the alignment measures
get_align_measures <- function(pr1,pr2,alignment,type = "all",superpos = TRUE){
  algn_measures <- list()
  #Distances of the aligned residues
  dist_vec <- get_pairwise_dist_align(pr1,pr2,alignment)
  #Alignment measures calculation
  #Number of alignmed residues (NaLi)
  if(type == "nali" | type == "all"){
    nali <- length(alignment$I) 
    if(type == "nali"){
      return(nali)
    }
    algn_measures$nali <- nali
  }
  #Root mean square deviation (RMSD)
  if(type == "rmsd" | type == "all"){
    rmsd <- sqrt((1/(length(dist_vec)))*sum(dist_vec^2))
    if(type == "rmsd"){
      return(rmsd)
    }
    algn_measures$rmsd <- rmsd
  }
  
  #Structure overlap (so)
  if(type == "so" | type == "all"){
    m <- min(dim(pr1)[1],dim(pr2)[1])
    n <- length(dist_vec)
    ones_count <- ifelse(dist_vec <= 3.5,rep(1,n),rep(0,n))
    so <- (100/m)*sum(ones_count)
    if(type == "so"){
      return(so)
    }
    algn_measures$so <- so
  }
 
  #SP Score
  if(type == "sp_score" | type == "all"){
    n <- length(dist_vec)
    summands <- ifelse(dist_vec < 2*3.5, (1/(1+(dist_vec / 3.5)^2)) - 0.2  ,rep(0,n))
    #Compute the L
    #First get core aligned residuals
    core_algn_res <- ifelse(dist_vec < 2*3.5, rep(1,n)  ,rep(0,n))
    core_res_1 <- alignment$I[as.logical(core_algn_res)]
    core_res_2 <- alignment$J[as.logical(core_algn_res)]
    #Get unaligned residuals from both proteins
    unalgn_1 <- setdiff(1:(dim(pr1)[1]),alignment$I)
    unalgn_2 <- setdiff(1:(dim(pr2)[1]),alignment$J)
    
    #For each core aligned residual, get number of "close" unaligned residuals
    #First protein
    core_res_neighbors_1 <- rep(0,length = length(core_res_1))
    if(length(core_res_1)>0){
      for(i in 1:length(core_res_1)){
        for(j in unalgn_1){
          if(sqrt(sum((pr1[core_res_1[i],] - pr1[j,])^2))<3*3.5){
            core_res_neighbors_1[i] <- core_res_neighbors_1[i]+1
          }
        }
        for(j in unalgn_2){
          if(sqrt(sum((pr1[core_res_1[i],] - pr2[j,])^2))<3*3.5){
            core_res_neighbors_1[i] <- core_res_neighbors_1[i]+1
          }
        }
      }
    }
    #Second protein
    core_res_neighbors_2 <- rep(0,length = length(core_res_2))
    if(length(core_res_1)>0){
      for(i in 1:length(core_res_2)){
        for(j in unalgn_1){
          if(sqrt(sum((pr2[core_res_2[i],] - pr1[j,])^2))<3*3.5){
            core_res_neighbors_2[i] <- core_res_neighbors_2[i]+1
          }
        }
        for(j in unalgn_2){
          if(sqrt(sum((pr2[core_res_2[i],] - pr2[j,])^2))<3*3.5){
            core_res_neighbors_2[i] <- core_res_neighbors_2[i]+1
          }
        }
      }
    }
    core_res_neighbors <- c(core_res_neighbors_1,core_res_neighbors_2)
    
    #L= number of core algn res + average unalgn residuals within 3d_0 from any core res
    if(length(core_res_neighbors)==0){
      mean_core_res_n <- 0
    } else{
      mean_core_res_n <- mean(core_res_neighbors)
    }
    L <- sum(core_algn_res) + mean_core_res_n
    #Calculate sp-score
    if(L > 0){
      sp_score <- (1/L^(1-0.3))*sum(summands)
    } else{
      sp_score <- 0
    }
    if(type == "sp_score"){
      return(sp_score)
    }
    algn_measures$sp_score <- sp_score
  }
  return(algn_measures)
}








