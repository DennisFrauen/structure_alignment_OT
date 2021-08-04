#Functions for artificial protein generation/ modification
library(uniformly)

#Randomly samples a 3-d point given an old coordinate
#new coordinate is sampled uniformly on a sphere around the old coordinate
#radius of the sphere is log-normal distributed
sample_coordinate <- function(old_coord,mean,sd){
  #Random distance to old coordinate
  radius <- exp(rnorm(1,mean = mean,sd = sd))
  new_coord <- runif_on_sphere(1,3,radius) + old_coord
  return(new_coord)
}

#Simulates an artificial protein structure by sampling coordinates using the function above
simulate_protein <- function(N,mean,sd){
  protein <- matrix(nrow = N,ncol = 3)
  protein[1,] <- c(0,0,0)
  for(i in 2:N){
    protein[i,] <- sample_coordinate(protein[i-1,],mean,sd)
  }
  return(protein)
}

#Simulates an artificial protein structure of size N on a sphere with radius r
simulate_protein_on_sphere <- function(N,r){
  protein <- matrix(nrow = N,ncol = 3)
  for(i in 1:N){
    protein[i,] <- runif_on_sphere(n = 1,d = 3,r = r)
  }
  return(protein)
}

#Perturbs a given protein structure by adding a vectors drawn from a multivariate normal
#distribution with mean 0 and standard deviation sd
perturb_protein <- function(protein, sd){
  for(i in 1:dim(protein)[1]){
    for(j in 1:dim(protein)[2]){
      protein[i,j] <- protein[i,j] + rnorm(1,0,sd)
    }
  }
  return(protein)
}


#Writes a protein structure to .pdb file
write_protein_to_pdb <- function(protein,name){
  bio3d::write.pdb(file = name,xyz = c(t(protein)), eleno = 1:dim(protein)[1],
                   resno = 1:dim(protein)[1],resid = rep("MET",dim(protein)[1]),
                   chain = rep("A",dim(protein)[1]),type = rep("ATOM",dim(protein)[1]),
                   o = rep(1,dim(protein)[1]),b = rep(50,dim(protein)[1]),
                   elesy = rep("C",dim(protein)[1]))
}
