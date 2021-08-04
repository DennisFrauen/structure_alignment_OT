#Code for Gromov-Wasserstein alignment method
library(matrixcalc)
library(FNN)
library(numDeriv)
library(tensor)

#Gromov-Wasserstein objective function
gromov_wstein_cost <- function(pi,Loss_tensor){
  return(frobenius.prod(tensor(A=Loss_tensor,B=pi,alongA = c(3,4),alongB = c(1,2)),pi))
}

#Greedy algorithm, converts a coupling pi to a (possibly non-sequentail) alignment (I,f)
coupling_to_alignNS <- function(pi){
  n <- dim(pi)[1]
  m <- dim(pi)[2]
  
  #Add artificial residue to second protein structure
  pi_added <- cbind(pi,rep(1/n,n)-pi%*%rep(1,m))
  colnames(pi_added) <- 1:(m+1)
  rownames(pi_added) <- 1:n
  
  #Compute alignment
  I <- c()
  J <- c()
  while(dim(pi_added)[1]>0 & dim(pi_added)[2]>0){
    #choose maximum coupling element location
    max_location <- as.matrix(which(pi_added == max(pi_added), arr.ind = TRUE))[1,]
    dims <- dim(pi_added)
    #Check whether residue i is matched to artificial residue
    if(max_location[2] < dims[2]){
      #add to alignment and remove corresponding row + column
      I <- c(I,as.integer(rownames(pi_added)[max_location[1]]))
      J <- c(J,as.integer(colnames(pi_added)[max_location[2]]))
      pi_added <- pi_added[-max_location[1],-max_location[2],drop=FALSE]
    } else{
      #remove row
      pi_added <- pi_added[-max_location[1],,drop=FALSE]
    }
  }
  return(list(I=I,J=J))
}

#Partial gromov wasserstein computation using Bregman gradient descent
#Input:   pr1,pr2: protein structures
#         mu_s,mu_t: probability measures (default to uniform)
#         alpha: fraction of mass transported
#         epsilon: entropic regularization parameter
#         tau: stepsize for Bregman gradient descent
#         iter_n: number of gradient step iterations
#         iter_k: number of iterations of Dykstra's algorithm
#         tol: indicates tolerance when to stop the algorithm
#         Loss_tensor: Gromov-Wasserstein loss tensor
#Output: Coupling pi
partial_gromov_entropic <-function(pr1,pr2,mu_s = c(),mu_t = c(),alpha=1,epsilon=0,tau=1,
                                   iter_n=100,iter_k=100,tol=10^-5,Loss_tensor = c()){
  n <- dim(pr1)[1]
  m <- dim(pr2)[1]
  #If measures are not specified, put uniform
  if(length(mu_s)==0){
    mu_s <- rep(1/n,n)
  }
  if(length(mu_t)==0){
    mu_t <- rep(1/m,m)
  }
  #Random initialization coupling
  Pi_n <- matrix(runif(n*m,min = 0,max = min(1/n,1/m)),nrow = n,ncol = m)
  print("Start Iteration")
  for(i in 1:iter_n){
    old <- Pi_n
    #Compute gradient of objective function
    if(length(Loss_tensor)!= 0){
      grad_xi <- tensor(A=Loss_tensor,B=Pi_n,alongA = c(3,4),alongB = c(1,2))
    } 
    #Input for Dykstras algorithm
    proj_input <- Pi_n * exp(-tau*(grad_xi + epsilon * log(Pi_n)))
    #Compute KL projection using Dykastras algorithm
    Pi_n <- dykstras(eta = proj_input,iter = iter_k,measure1 = mu_s,measure2 = mu_t,alpha = alpha)
    #Check whether computation was succesful
    if(anyNA(Pi_n)==TRUE){
      print("Error, step size too large")
      return(Pi_n)
    }
    #Check Tolerance
    if(base::norm(Pi_n - old,"I")<tol){
      return(Pi_n)
    }
  }
  return(Pi_n)
}

#Compute Gromov-Wasserstein loss tensor for protein structures pr1, pr2 (for p = 1)
get_loss_tensor <- function(pr1,pr2){
  #Loss function
  loss_gromov <- function(x,y){
    return(abs(x-y))
  }
  n <- dim(pr1)[1]
  m <- dim(pr2)[1]
  #Get inter-residue distance matrices
  C_s <- matrix(unlist(get_inter_res_dist(pr1)),nrow = n,ncol = n)
  C_t <- matrix(unlist(get_inter_res_dist(pr2)),nrow = m,ncol = m)
  #Empty tensor
  Loss_tensor <- array(dim = c(n,m,n,m))
  #Fill loss tensor
  for(i in 1:n){
    for(i_bar in 1:m){
      for(j in 1:n){
        for(j_bar in 1:m){
          Loss_tensor[i,i_bar,j,j_bar] <- loss_gromov(C_s[i,j],C_t[i_bar,j_bar])
        }
      }
    }
  }
  return(Loss_tensor)
}


#KL projections-------------------------------------------------------------------------------------
proj_C1 <- function(gamma,mu_s){
  ones_n <- rep(1,dim(gamma)[1])
  ones_m <- rep(1,dim(gamma)[2])
  proj <- diag(as.vector(pmin(mu_s/ (gamma %*% ones_m),ones_n)))
  return(proj%*%gamma)
}

proj_C2 <- function(gamma,mu_t){
  ones_n <- rep(1,dim(gamma)[1])
  ones_m <- rep(1,dim(gamma)[2])
  proj <- diag(as.vector(pmin(mu_t/ (t(gamma) %*% ones_n),ones_m)))
  return(gamma%*%proj)
}

proj_C3 <- function(gamma,alpha){
  ones_n <- rep(1,dim(gamma)[1])
  ones_m <- rep(1,dim(gamma)[2])
  proj <- (alpha / (t(ones_n)%*%gamma%*%ones_m))[1,1]
  return(gamma*proj)
}

#Dykstras algorithm
dykstras <- function(eta,iter,measure1,measure2,alpha,tol=10^-6){
  n <- dim(eta)[1]
  m <- dim(eta)[2]
  q_list <- list(a = ones(n,m),b = ones(n,m),c = ones(n,m))
  gam <- eta
  for(i in 1:iter){
    proj_nr <- mod(i-1,3) +1
    gam_old <- gam
    if(proj_nr==1){
      #First projection
      gam <- proj_C1(gam*q_list$a,measure1)
      q_list$a <- (gam_old/gam) * q_list$a
    } else{
      if(proj_nr == 2){
        #Second projection
        gam <- proj_C2(gam*q_list$b,measure2)
        q_list$b <- (gam_old/gam) * q_list$b
      } else{
        if(proj_nr==3){
          #Third projection
          gam <- proj_C3(gam*q_list$c,alpha)
          q_list$c <- (gam_old/gam) * q_list$c
        }
      }
    }
  }
  return(gam)
}

#Iterative Bregman iterations
bregman <- function(eta,iter,measure1,measure2,alpha){
  n <- dim(eta)[1]
  m <- dim(eta)[2]
  gam <- eta
  for(i in 1:iter){
    proj_nr <- mod(i-1,3) +1
    gam_old <- gam
    if(proj_nr==1){
      #First projection
      gam <- proj_C1(gam,measure1)
    } else{
      if(proj_nr == 2){
        #Second projection
        gam <- proj_C2(gam,measure2)
      } else{
        if(proj_nr==3){
          #Third projection
          gam <- proj_C3(gam,alpha)
        }
      }
    }
  }
  return(gam)
}

#Quantized partial gromov wasserstein --------------------------------------------------------------
#Main function: Computes quantized partial Gromov-Wasserstein coupling
#Input:     pr1,pr2: protein structures
#           runtimes: times bregman descent is applied to fin the global coupling (avoiding local minima)
#           m_repr: number of representatives
#           alpha: fraction of transported mass
#           tau: learning rate for Bregman descent
#           iter_n: number of iterations for Bregman descent
quantized_partial_gromov <- function(pr1,pr2,runtimes=1,m_repr,alpha=1,tau=0.03,
                                     iter_n=200,quantize=TRUE){
  n <- dim(pr1)[1]
  m <- dim(pr2)[1]
  #Determine representatives (most uniformly over both proteins)
  d_1 <- ceil(n/m_repr)
  d_2 <- ceil(m/m_repr)
  repr1 <- vector(length = m_repr)
  repr2 <- vector(length = m_repr)
  for(i in 1:m_repr){
    repr1[i] <- 1 + (i-1)*d_1
    repr2[i] <- 1 + (i-1)*d_2
  }
  remove_1 <- setdiff(1:n,repr1)
  remove_2 <- setdiff(1:m,repr2)
  x1 <- length(remove_1) - (n - m_repr)
  x2 <- length(remove_2) - (m - m_repr)
  remove_1 <- setdiff(remove_1,sample(remove_1,size = x1,replace = FALSE))
  remove_2 <- setdiff(remove_2,sample(remove_2,size = x2,replace = FALSE))
  repr1 <- setdiff(1:n,remove_1)
  repr2 <- setdiff(1:m,remove_2)
  
  #Create pointed partitions
  partition1 <- get_partition(pr1,m_repr,repr1)
  partition2 <- get_partition(pr2,m_repr,repr2)
  #Get quantization coupling
  Pi <- get_quantized_coupling(pr1 = pr1,pr2 = pr2,partition1 = partition1,
                                    partition2 = partition2,runtimes = runtimes,
                                      alpha = alpha,tau = tau,iter_n = iter_n)
  #return coupling
  return(Pi)
}

#Creates quantization couplings from given pointed partitions
get_quantized_coupling <- function(pr1,pr2,partition1,partition2,alpha=1,tau=0.03,
                                   iter_n=200,runtimes){
  #Global coupling
  print("Global coupling")
  Pi_global <- get_global_coupling(pr1 = pr1,pr2 = pr2,partition1 = partition1,
                                   partition2 = partition2, alpha = sqrt(alpha),
                                   tau = tau,iter_n = iter_n,runtimes = runtimes)
  #Local couplings
  print("Local couplings")
  Pi_loc <- get_local_couplings(pr1,pr2,partition1,partition2,sqrt(alpha))
  #Create quantized coupling by combining global with local couplings
  Pi <- matrix(0,nrow = dim(pr1)[1],ncol = dim(pr2)[1])
  for(p in 1:length(partition1$repr)){
    for (q in 1:length(partition1$repr)) {
      Pi <- Pi + as.matrix(Pi_loc[[p]][[q]]) * Pi_global[p,q]
    }
  }
  return(as.matrix(Pi))
}

#Creates a partition given a protein structure and representative residues
get_partition <- function(protein,m,repr = c()){
  point_ind <- c(1:dim(protein)[1])
  #If no representatives given, sample representatives randomly
  if(length(repr)==0){
    repr <- sort(sample(point_ind,m))
  }
  #indices of residues which are no representatives
  point_ind <- setdiff(point_ind,repr)
  #Create partition classes
  partition <- list()
  for(i in 1:m){
    partition[[i]] <- c(repr[i])
  }
  #Asign all residues to the closest partition class w.r.t. the euclidean metric
  if(length(point_ind)>0){
    for (i in 1:length(point_ind)) {
      dist_to_repr <- apply(X=protein[repr,],FUN = eucl_dist,y = protein[point_ind[i],],MARGIN = 1)
      closest_repr <- min(which(dist_to_repr==min(dist_to_repr)))
      partition[[closest_repr]] <- sort(c(partition[[closest_repr]],point_ind[i]))
    }
  }
  #return partition classes and representatives
  return(list(repr = repr, partition = partition))
}

#Computes the global coupling given pointed partitions
get_global_coupling <- function(pr1,pr2,partition1,partition2,alpha=1,tau=0.03,
                                iter_n=200,runtimes=1){
  #Extract reprisentative residues
  repr1 <- pr1[partition1$repr,]
  repr2 <- pr2[partition2$repr,]
  #Push forward measures
  mu_s <- unlist(lapply(X = partition1$partition,FUN = length)) / dim(pr1)[1]
  mu_t <- unlist(lapply(X = partition2$partition,FUN = length)) / dim(pr2)[1]
  
  #Get loss tensor for representatives
  Loss_tensor <- get_loss_tensor(repr1,repr2)
  #Get gw coupling, try different starting points of gradient descent (avoid local minima)
  gw_couplings <- list()
  gw_cost <- vector(length = runtimes)
  k <- 1
  while(k <= runtimes){
    gw_couplings[[k]] <- partial_gromov_entropic(pr1 =repr1,pr2 =repr2,
                alpha = alpha,iter_n = 200,Loss_tensor = Loss_tensor,tau = tau,
                mu_s = mu_s,mu_t = mu_t)
    if(anyNA(gw_couplings[[k]])==FALSE){
      gw_cost[k] <- gromov_wstein_cost(gw_couplings[[k]], Loss_tensor)
      k <- k + 1
    } else{
      learn_rate <- learn_rate - 0.02
      print(learn_rate)
    }
  }
  #Choose coupling minimizing the Gromov-Wasserstein cost
  best_coupl <- gw_couplings[[min(which(gw_cost==min(gw_cost)))]]
  return(best_coupl)
}

#Computes an partial optimal transport plan given cost matrix C, measures mu_s and mu_t and alpha
partial_OT <- function(C,mu_s,mu_t,alpha=1){
  n <- dim(C)[1]
  m <- dim(C)[2]
  if(n > 1 | m > 1){
    #check whether partial or full transport needs to be applied
    if(alpha < 1){
      #Extent Cost matrix to use optimal transport
      A <- max(C) + 1
      C <- rbind(C,rep(eta,m))
      C <- cbind(C, c(rep(eta,n),A))
      #Extent measures
      mu_s = c(mu_s,1-alpha)
      mu_t = c(mu_t,1-alpha)
      #compute extended plan
      plan <- transport(a = mu_s,b = mu_t,costm = C,fullreturn = TRUE)$primal
      #restrict
      plan <- plan[-(n+1),-(m+1),drop=FALSE]
      return(plan[-(n+1),-(m+1)])
    } else{
      if(n > 1 & m > 1){
        plan <- transport(a = mu_s,b = mu_t,costm = C,fullreturn = TRUE)$primal
        return(plan)
      } else{
        if(n == 1){
          plan <- matrix(nrow = 1,ncol = m)
          plan[1,] <- mu_t
          return(plan)
        } else{
          plan <- matrix(nrow = n,ncol = 1)
          plan[,1] <- mu_s
          return(plan)
        }
      }
    }
  } else{
    plan <- matrix(alpha,1,1)
    return(plan)
  }

}

#Computes the local couplings given pointed partitions
get_local_couplings <- function(pr1,pr2,partition1,partition2,alpha=1){
  #Dimensions
  n <- dim(pr1)[1]
  m <- dim(pr2)[1]
  m_rep <- length(partition1$repr)
  #Local couplings, two dimensional list (mrep x mrep) of sparse matrices (n x m)
  loc_couplings <- list()
  for(i in 1:m_rep){
    loc_couplings[[i]] <- list()
  }

  #Compute couplings
  for (p in 1:m_rep) {
    for(q in 1:m_rep){
      #Partition classes
      U_p <- partition1$partition[[p]]
      V_q <- partition2$partition[[q]]
      #Measures
      mu_Up <- rep(1/length(U_p), length(U_p))
      mu_Vq <- rep(1/length(V_q), length(V_q))
      #Define Cost matrix
      C <- matrix(nrow=length(U_p),ncol=length(V_q))
      for(x in 1:length(U_p)){
        for(y in 1:length(V_q)){
          C[x,y] <- (eucl_dist(pr1[U_p[x],], pr1[partition1$repr[p],]) - 
                       eucl_dist(pr2[V_q[y],], pr2[partition2$repr[q],]))^2
        }
      }
      #Solve partial optimal transport problem
      Pi <- partial_OT(C=C,mu_s = mu_Up,mu_t = mu_Vq,alpha = alpha)
      #Extent to full coupling (sparse matrix)
      Pi_tilde <- matrix(0,nrow = n,ncol = m)
      Pi_tilde[U_p,V_q] <- Pi
      Pi_tilde <- Matrix(Pi_tilde,sparse = TRUE)
      #Add to local coupling list
      loc_couplings[[p]][[q]] <- Pi_tilde
    }
  }
  return(loc_couplings)
  
}

