#Code to compute several distances, e.g. partial Wasserstein, median difference, 
#euclidean distances after superimposition,...

library(transport)
library(tidyverse)
library(Rsymphony)
library(Barycenter)


#Euclidean distance
eucl_dist <- function(x,y,p=2){
  return((sum((x-y)^p))^(1/p))
}

#Partial Wasserstein distance
partial_Wasserstein1d <- function(x,y,p=1,alpha=1){
  #dimensions
  n <- length(x)
  m <- length(y)
  #Create cost matrix c (nxm)
  C <- sapply(x,"-",y) %>% abs() %>% t()
  C <- C^p
  #Extent C (to reduce partial transport to full transport)
  eta <- 0
  A <- max(C)
  C <- rbind(C,rep(eta,m))
  C <- cbind(C, c(rep(eta,n),2*eta + A))
  #Distribution (with auxiliary points added)
  a = c(rep(1/n,n),1-alpha)
  b = c(rep(1/m,m),1-alpha)
  #Compute coupling
  plan <- transport(a = a,b = b,costm = C,fullreturn = TRUE)
  #return transport cost
  return(plan$cost^(1/p))
}

#Median difference
median_diff1d <- function(x,y){
  return(abs(eucl_dist(median(x),median(y),1)))
}

#Euclidean distance matrix (Inter-protein-residue distances) after superimposition
eukl_dist_matrix <- function(pr1,pr2,alignment){
  N <- dim(pr1)[1]
  M <- dim(pr2)[1]
  I <- alignment$I
  J <- alignment$J
  #Superimposition
  kabsch <- pracma::kabsch(A = t(pr1[I,]),B = t(pr2[J,]))
  Trans_mat <- matrix(kabsch$R, nrow=dim(pr1)[1], ncol=3, byrow=TRUE)
  pr1 <- pr1 %*% t(kabsch$U) +Trans_mat
  #Create empty matrix
  eucl_dist_matr <- matrix(nrow = N,ncol = M)
  #Compute inter-protein-residue distances
  for (i in 1:M){
    eucl_dist_matr[,i] <- apply(FUN = eucl_dist,X = pr1,MARGIN = 1,y = pr2[i,])
  }
  return(eucl_dist_matr)
}

#Creates Cost matrix for monotone rearrangement depending on cost functional dist
#Possible values for dist: wasserstein, partial_wasserstein, median_diff
wstein_dist_matr <- function(ird1,ird2,p=1,alpha=1, dist = "wasserstein"){
  #Create empty matrix
  wstein_matrix <- matrix(NA,nrow = length(ird1),ncol = length(ird2))
  #Fill column wise, depending on the cost functional
  #Wasserstein distance (alpha = 1)
  if(dist == "wasserstein" | alpha == 1){
    for (i in 1:length(ird2)){
      wstein_matrix[,i] <- sapply(FUN = wasserstein1d,X = ird1,b = ird2[[i]],p=p)
    }
  } 
  #Partial Wassersetin distance (alpha < 1)
  if(dist == "partial_wasserstein" & alpha < 1){
    for (i in 1:length(ird2)){
      wstein_matrix[,i] <- sapply(FUN = partial_Wasserstein1d,X = ird1,y = ird2[[i]],
                                  p=p,alpha=alpha)
    }
  } 
  #Median difference
  if(dist == "median_diff"){
    for (i in 1:length(ird2)){
      wstein_matrix[,i] <- sapply(FUN = median_diff1d,X = ird1,y = ird2[[i]])
    }
  }
  return(wstein_matrix)
}

#Computes inter-residue distributions for given protein structure p
#Input: dist_cutoff: larger inter-residue distances will be cut off, 0 indicates no cutoff
#       seq_cutoff: inter-residue distances between residues larger apart in the sequence
#                   will be cut off, 0 indicates no cutoff
#       remove_res: distances to indicated residues will be cut off (for adjusted 
#                   monotone rearrangement)
#Output: List containing all inter-residue distributions (as vectors)
get_inter_res_dist <- function(p,dist_cutoff = 0,seq_cutoff=0,remove_res = c()){
  len <- dim(p)[1]
  #Function to compute inter-residue distance bewteen two residues
  inter_res_dist <- function(p,res1,res2){
    return(dist(rbind(p[res1,],p[res2,]))%>% as.numeric())
  }
  #empty inter-residue distance matrix
  ird_matr <- matrix(0,nrow = len,ncol = len)
  #Fill ird matrix column wise
  for (i in 1:len){
    ird_matr[i:len,i] <- sapply(inter_res_dist,X = i:len,res1 = i,p = p)
  }
  ird_matr[upper.tri(ird_matr)] <- t(ird_matr)[upper.tri(t(ird_matr))]
  #Change rownames to residue numbers
  rownames(ird_matr) <- paste(1:len)
  #Remove specified columns
  if(length(remove_res)!=0){
    ird_matr <- ird_matr[-remove_res,]
  }
  select_res <- as.integer(rownames(ird_matr))
  #Convert ird matrix to a list of distributions
  ird <- split(ird_matr,rep(1:ncol(ird_matr),each = nrow(ird_matr)))
  #Sequence cutoff
  #Only consider inter residue distances withing neighborhood in the sequence
  if(seq_cutoff > 0){
    #Go through all residues
    for(i in 1:length(ird)){
      left_res <- max(1,i-seq_cutoff)
      right_res <- min(length(ird),i + seq_cutoff)
      left_ind <- min(which(left_res <= select_res))
      right_ind <- max(which(right_res >= select_res))
      ird[[i]] <- ird[[i]][left_ind:right_ind]
    }
  }
  #Distance cutoff
  if(dist_cutoff>0){
    #Cut off by distance
    #Go through all ird distributions
    for(i in 1:length(ird)){
      #Cut off too large ird distances
      if(length(which(ird[[i]] <= dist_cutoff))>1){
        ird[[i]] <- ird[[i]][ird[[i]] <= dist_cutoff]
      }
    }
  }
  return(ird)
}

#get inter-protein-residue distances of aligned residues (before/ after superposition)
get_pairwise_dist_align <- function(pr1,pr2,alignment,superpos=TRUE){
  #First apply kabsch algorithm if necissary, to perform superimposition
  if(superpos==TRUE){
    #apply katsch on aligned residues to get rotation matrix and translation vector
    kabsch <- pracma::kabsch(A = t(pr1[alignment$I,]),B = t(pr2[alignment$J,]))
    #translate/ rotate first protein accordingly
    Trans_mat <- matrix(kabsch$R, nrow=dim(pr1)[1], ncol=3, byrow=TRUE)
    pr1 <- pr1 %*% t(kabsch$U) +Trans_mat
  }
  
  #Differences of the aligned coordinate vectors
  coord_diff <- as.matrix(pr1[alignment$I,] - pr2[alignment$J,])
  if(dim(coord_diff)[2] == 1){
    coord_diff <- t(coord_diff)
  }
  #Distances of the aligned residues
  dist_vec <- rowNorms(coord_diff,method = "euclidean", p = 2)
  return(dist_vec)
}
