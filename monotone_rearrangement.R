source("distances.R")

#Monotone rearrangement-------------------------------------------------------------------
#Input: Loss matrix D (either median differences or partial Wasserstein distances)
#       kappa (optionl) number of aligned residues,
#       if kappa = 0, the full tensor H is returned
#Output: Monotone rearrangement tensor H
monotone_rearr <- function(D,kappa=0){
  #Sequence lenghts
  N <- dim(D)[1]
  M <- dim(D)[2]
  
  #kappa = 0 <-> no kappa provided -> calculate for all kappas and full H matrix
  if(kappa == 0){
    k_lim <- min(N,M)+1
    i_lim <- N+1
    j_lim <- M+1
  } else{
    k_lim <- kappa+1
  }
  
  #Dynamic Programming matrix (H matrix)
  H <- array(dim = c(N+1,M+1,k_lim))
  H[,,1] <- matrix(0,N+1,M+1)
  
  for(k in 2:k_lim){
    #If kappa provided -> don't calculate full H matrix
    if(kappa != 0){
      i_lim <- N - kappa + k
      j_lim <- M - kappa + k
    }
    #Calculate elements of H by recursion formula
    for(i in k:i_lim){
      for(j in k:j_lim){
        #Create list to minimise over
        minlist <- list(a = H[i,j-1,k],b = H[i-1,j,k],
                     c = H[i-1,j-1,k-1] + D[i-1,j-1])
        #Check that k <= min(i,j)
        if(j == k){
          minlist <- minlist[names(minlist) != "a"]
        }
        if(i == k){
          minlist <- minlist[names(minlist) != "b"]
        }
        minvec <- unlist(minlist)
        H[i,j,k] <- min(minvec)
      }
    }
  }
  return(H)
}

#Monotone rearrangement traceback procedure-------------------------------------------------
#Input: Monotone rearrangement tensor H
#Output: Alignment corresponding to H
get_optimal_align <- function(H){
  i <- dim(H)[1]
  j <- dim(H)[2]
  k <- dim(H)[3]
  #Alignment domains
  A1 <- vector(length = k-1)
  A2 <- vector(length = k-1)
  #Traceback
  while(k > 1){
    #Check if a is infinity
    if(j == k){
      #if b is infinity
      if(i == k){
        #take c
        A1[k-1] <- i-1
        A2[k-1] <- j-1
        i <- i-1
        j <- j-1
        k <- k-1
      } else{
        #check if b is optimal
        if(H[i,j,k] == H[i-1,j,k]){
          #take b
          i <- i-1
        } else{
          #take c
          A1[k-1] <- i-1
          A2[k-1] <- j-1
          i <- i-1
          j <- j-1
          k <- k-1
        }
      }
    } else{
      #check if a is optimal
      if(H[i,j,k] == H[i,j-1,k]){
        #take a
        j <- j-1
      } else{
        #check if b is infinity
        if(i == k){
          #take c
          A1[k-1] <- i-1
          A2[k-1] <- j-1
          i <- i-1
          j <- j-1
          k <- k-1
        } else{
          #if b is optimal
          if(H[i,j,k] == H[i-1,j,k]){
            #take b
            i <- i-1
          } else{
            #take c
            A1[k-1] <- i-1
            A2[k-1] <- j-1
            i <- i-1
            j <- j-1
            k <- k-1
          }
        }
      }
    }
  }
  return(list(I = A1, J = A2))
}




