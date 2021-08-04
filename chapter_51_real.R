#Code for all results based on real protein data (section 5.1)
#For running, the source folders of the protein structures might need to be adjusted
library(bio3d)
source("structure_allignment.R")
source("gromov_wstein.R")

#Load protein structures-------------------------------------------------------------------
#Function for extracting coordinate information of PDB files
extract_coordinates <- function(pdb_list){
  #Extract coordinates
  protein_list <- list()
  
  for(i in 1:length(pdb_list)){
    protein_list[[i]] <- pdb_list[[i]]$atom %>% select(x,y,z) %>% 
      slice(which(pdb_list[[i]]$atom$elety=="CA")) %>% as.matrix()
  }
  return(protein_list)
}
#Import Phycocyanin family
load_phycocyanin <- function(){
  #PDB files
  phycocyanin_pdb <- list()
  phycocyanin_pdb[[1]] <-  bio3d::read.pdb(file = "./protein_structures/Phycocyanin/1all.pdb",ATOM.only = TRUE)
  phycocyanin_pdb[[2]] <-  bio3d::read.pdb(file = "./protein_structures/Phycocyanin/1b33.pdb",ATOM.only = TRUE)
  phycocyanin_pdb[[3]] <-  bio3d::read.pdb(file = "./protein_structures/Phycocyanin/1b8d.pdb",ATOM.only = TRUE)
  phycocyanin_pdb[[4]] <-  bio3d::read.pdb(file = "./protein_structures/Phycocyanin/1cpc.pdb",ATOM.only = TRUE)
  phycocyanin_pdb[[5]] <-  bio3d::read.pdb(file = "./protein_structures/Phycocyanin/1lia.pdb",ATOM.only = TRUE)
  phycocyanin_pdb[[6]] <-  bio3d::read.pdb(file = "./protein_structures/Phycocyanin/1phn.pdb",ATOM.only = TRUE)
  phycocyanin_pdb[[7]] <-  bio3d::read.pdb(file = "./protein_structures/Phycocyanin/1qgw.pdb",ATOM.only = TRUE)
  return(phycocyanin_pdb)
}
phycocyanin <- extract_coordinates(load_phycocyanin())
#Import ferredoxin family 
load_ferredoxin <- function(){
  #PDB files
  ferredoxin_pdb <- list()
  ferredoxin_pdb[[1]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/1blu.pdb",ATOM.only = TRUE)
  ferredoxin_pdb[[2]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/1bwe.pdb",ATOM.only = TRUE)
  ferredoxin_pdb[[3]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/1clf.pdb",ATOM.only = TRUE)
  ferredoxin_pdb[[4]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/1dur.pdb",ATOM.only = TRUE)
  ferredoxin_pdb[[5]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/1dwl.pdb",ATOM.only = TRUE)
  ferredoxin_pdb[[6]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/1fxd.pdb",ATOM.only = TRUE)
  ferredoxin_pdb[[7]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/1fxr.pdb",ATOM.only = TRUE)
  ferredoxin_pdb[[8]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/1h98.pdb",ATOM.only = TRUE)
  ferredoxin_pdb[[9]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/1hfe.pdb",ATOM.only = TRUE)
  ferredoxin_pdb[[10]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/1k0t.pdb",ATOM.only = TRUE)
  ferredoxin_pdb[[11]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/1vjw.pdb",ATOM.only = TRUE)
  ferredoxin_pdb[[12]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/1xer.pdb",ATOM.only = TRUE)
  ferredoxin_pdb[[13]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/2fdn.pdb",ATOM.only = TRUE)
  ferredoxin_pdb[[14]] <-  bio3d::read.pdb(file = "./protein_structures/ferredoxin/7fd1.pdb",ATOM.only = TRUE)
  return(ferredoxin_pdb)
}
ferrodoxin <- extract_coordinates(load_ferredoxin())

#Import GNAT family
load_gnat <- function(){
  #PDB files
  phycocyanin_pdb <- list()
  phycocyanin_pdb[[1]] <-  bio3d::read.pdb(file = "./protein_structures/Gnat/1b87.pdb",ATOM.only = TRUE)
  phycocyanin_pdb[[2]] <-  bio3d::read.pdb(file = "./protein_structures/Gnat/1bo4.pdb",ATOM.only = TRUE)
  phycocyanin_pdb[[3]] <-  bio3d::read.pdb(file = "./protein_structures/Gnat/1cjw.pdb",ATOM.only = TRUE)
  phycocyanin_pdb[[4]] <-  bio3d::read.pdb(file = "./protein_structures/Gnat/1cm0.pdb",ATOM.only = TRUE)
  phycocyanin_pdb[[5]] <-  bio3d::read.pdb(file = "./protein_structures/Gnat/1qst.pdb",ATOM.only = TRUE)
  phycocyanin_pdb[[6]] <-  bio3d::read.pdb(file = "./protein_structures/Gnat/1ygh.pdb",ATOM.only = TRUE)
  return(phycocyanin_pdb)
}
gnat <- extract_coordinates(load_gnat())

#Import Lipocalin family
load_lipocalin <- function(){
  #PDB files
  lipocalin_pdb <- list()
  lipocalin_pdb[[1]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1aqb.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[2]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1bbp.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[3]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1beb.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[4]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1bj7.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[5]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1dzk.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[6]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1e5p.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[7]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1epa.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[8]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1ew3.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[9]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1exs.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[10]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1hn2.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[11]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1i4u.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[12]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1iiu.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[13]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1jv4.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[14]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/1qqs.pdb",ATOM.only = TRUE)
  lipocalin_pdb[[15]] <-  bio3d::read.pdb(file = "./protein_structures/Lipocalin/2a2u.pdb",ATOM.only = TRUE)

  return(lipocalin_pdb)
}
lipocalin <- extract_coordinates(load_lipocalin())

#Import legume_lectin family
load_legume <- function(){
  #PDB files
  legume_pdb <- list()
  legume_pdb[[1]] <-  bio3d::read.pdb(file = "./protein_structures/legume_lectin/1avb.pdb",ATOM.only = TRUE)
  legume_pdb[[2]] <-  bio3d::read.pdb(file = "./protein_structures/legume_lectin/1lec.pdb",ATOM.only = TRUE)
  legume_pdb[[3]] <-  bio3d::read.pdb(file = "./protein_structures/legume_lectin/1lte.pdb",ATOM.only = TRUE)
  legume_pdb[[4]] <-  bio3d::read.pdb(file = "./protein_structures/legume_lectin/1lu1.pdb",ATOM.only = TRUE)
  legume_pdb[[5]] <-  bio3d::read.pdb(file = "./protein_structures/legume_lectin/1sbf.pdb",ATOM.only = TRUE)
  legume_pdb[[6]] <-  bio3d::read.pdb(file = "./protein_structures/legume_lectin/1wbl.pdb",ATOM.only = TRUE)
  return(legume_pdb)
}
legume_lectin <- extract_coordinates(load_legume())

#Import RIPC proteins
load_ripcs <- function(){
  #PDB files
  ripc_pdb <- list()
  ripc_pdb[[1]] <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d1biaa1.pdb",ATOM.only = TRUE)
  ripc_pdb[[2]] <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d1qdma1.pdb",ATOM.only = TRUE)
  ripc_pdb[[3]] <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d1an9a1.pdb",ATOM.only = TRUE)
  ripc_pdb[[4]] <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d148le_.pdb",ATOM.only = TRUE)
  ripc_pdb[[5]] <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d1cpoa2.pdb",ATOM.only = TRUE)
  ripc_pdb[[6]] <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d3chya_.pdb",ATOM.only = TRUE)
  ripc_pdb[[7]] <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d3coxa2.pdb",ATOM.only = TRUE)
  ripc_pdb[[8]] <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d4clna_.pdb",ATOM.only = TRUE)
  ripc_pdb[[9]] <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d1b09a_.pdb",ATOM.only = TRUE)
  ripc_pdb[[10]] <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d1puja_.pdb",ATOM.only = TRUE)
  ripc_pdb[[11]] <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d1cpoa2.pdb",ATOM.only = TRUE)
  ripc_pdb[[12]] <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d4clna_.pdb",ATOM.only = TRUE)
  return(ripc_pdb)
}
ripcs <- extract_coordinates(load_ripcs())

#Defining alignment pairs
match1 <- legume_lectin[c(1,3,5)]
match2 <- legume_lectin[c(2,4,6)]
match1[4:10] <- lipocalin[c(1,3,5,7,9,11,13)]
match2[4:10] <- lipocalin[c(2,4,6,8,10,12,14)]
match1[11:13] <- gnat[c(1,3,5)]
match2[11:13] <- gnat[c(2,4,6)]
match1[14:20] <- ferrodoxin[c(1,3,5,7,9,11,13)]
match2[14:20] <- ferrodoxin[c(2,4,6,8,10,12,14)]
match1[21:23] <- phycocyanin[c(1,3,5)]
match2[21:23] <- phycocyanin[c(2,4,6)]
match1[24:30] <- ripcs[c(1,3,5,7,9,11,1)]
match2[24:30] <- ripcs[c(2,4,6,8,10,12,6)]

#Experiments--------------------------------------------------------------------------------
#Structure alignment methods are applied on protein data for different hyperparameters
#
#Monotone rearrangement with median difference
lambda_seq <- exp(seq(log(0.0001), log(1000), length.out = 20))
#Matrix to save alignment measures
meas_matr <- matrix(nrow = length(match1),ncol = 5)
#Run over all protein pairs
for (i in 1:30) {
  print(i)
  print(Sys.time())
  so_vec <- vector(length = length(lambda_seq))
  #Run for different lambda
  for(l in 1:length(lambda_seq)){
    alignment <- align_proteins(match1[[i]],match2[[i]],return_meas = TRUE,lambda = lambda_seq[l],iter = 10,
                                dist = "median_diff",print_info = FALSE)
    so_vec[l] <- alignment$measures$so[length(alignment$measures$so)]
  }
  #Choose optimal lambda
  lambda_opt <- max(lambda_seq[which(so_vec==max(so_vec))])
  print(lambda_opt)
  start_time <- Sys.time()
  alignment <- align_proteins(match1[[i]],match2[[i]],return_meas = TRUE,lambda = lambda_opt,iter = 10,
                              dist = "median_diff",print_info = FALSE)
  end_time <- Sys.time()
  #Save alignment measures
  meas_matr[i,1] <- alignment$measures$so[length(alignment$measures$so)]
  meas_matr[i,2] <- alignment$measures$rmsd[length(alignment$measures$rmsd)]
  meas_matr[i,3] <- alignment$measures$nali[length(alignment$measures$nali)]
  meas_matr[i,4] <- end_time - start_time
  meas_matr[i,5] <- lambda_opt
}

#Monotone rearrangement with wasserstein
lambda_seq <- exp(seq(log(0.0001), log(1000), length.out = 20))
alpha_seq <- seq(0.5,1,0.1)
#Matrix to save alignment measures
meas_matr_wstein <- matrix(nrow = length(match1),ncol = 6)
#Run on all protein pairs
for (i in 30:30) {
  print(i)
  print(Sys.time())
  so_matr <- matrix(nrow = length(alpha_seq),ncol = length(lambda_seq))
  nali_matr <- matrix(nrow = length(alpha_seq),ncol = length(lambda_seq))
  #Run for different alpha
  for(k in 1:length(alpha_seq)){
    print(paste("Alpha =",alpha_seq[k]))
    #Compute partial Wasserstein matrix for alpha
    ird1 <- get_inter_res_dist(match1[[i]])
    ird2 <- get_inter_res_dist(match2[[i]])
    wstein_matr <- wstein_dist_matr(ird1,ird2,p = 1,alpha = alpha_seq[k],dist="partial_wasserstein")
    #Run for different lambda
    for(l in 1:length(lambda_seq)){
      print(paste("Lambda =",lambda_seq[l]))
      #Alignment
      alignment <- align_proteins(match1[[i]],match2[[i]],return_meas = TRUE,lambda = lambda_seq[l],
                                  iter = 10, wstein_matr = wstein_matr,
                                  dist = "partial_wasserstein",print_info = FALSE,alpha = alpha_seq[k])
      so_matr[k,l] <- alignment$measures$so[length(alignment$measures$so)]
      nali_matr[k,l] <- alignment$measures$nali[length(alignment$measures$nali)]
    }
  }
  #Choose optimal alpha and lambda
  opt_so <- which(so_matr == max(so_matr), arr.ind = TRUE)
  nali_vec <- nali_matr[opt_so]
  opt_nali <- max(which(nali_vec == max(nali_vec)))
  opt_par <- opt_so[opt_nali,]
  print(opt_par)
  print(so_matr)
  print(nali_matr)
  alpha_opt <- alpha_seq[opt_par[1]]
  lambda_opt <- lambda_seq[opt_par[2]]
  print(paste(alpha_opt,lambda_opt))
  
  start_time <- Sys.time()
  alignment <- align_proteins(match1[[i]],match2[[i]],return_meas = TRUE,lambda = lambda_opt,iter = 10,
                              dist = "partial_wasserstein",print_info = FALSE,alpha = alpha_opt)
  end_time <- Sys.time()
  #Save alignment measures
  meas_matr_wstein[i,1] <- alignment$measures$so[length(alignment$measures$so)]
  meas_matr_wstein[i,2] <- alignment$measures$rmsd[length(alignment$measures$rmsd)]
  meas_matr_wstein[i,3] <- alignment$measures$nali[length(alignment$measures$nali)]
  meas_matr_wstein[i,4] <- end_time - start_time
  meas_matr_wstein[i,5] <- lambda_opt
  meas_matr_wstein[i,6] <- alpha_opt
}

#Nonsequential Gromov-Wasserstein
meas_matr_gromov <- matrix(nrow = length(match1),ncol = 5)
alpha_seq <- seq(0.3,1,0.05)
#Run on all protein pairs with both structure sizes less than 80
for(i in 20:20){
  set.seed(123512)
  #Get loss tensor for protein pair
  Loss_tensor <- get_loss_tensor(match1[[i]],match2[[i]])
  so_seq <- alpha_seq
  rmsd_seq <- alpha_seq
  nali_seq <- alpha_seq
  runtimes <- 15
  #Run for different alpha
  for (j in 1:length(alpha_seq)) {
    print(j)
    coupl_gromov <- list()
    gw_cost <- vector(length = runtimes)
    learn_rate <- 0.3
    #Run many times to avoid local minima
    k <- 1
    while(k <= runtimes){
      coupl_gromov[[k]] <- partial_gromov_entropic(pr1 =match1[[i]],pr2 =match2[[i]],
                                                   alpha = alpha_seq[j],iter_n = 200,
                                                   Loss_tensor = Loss_tensor,tau = learn_rate)
      #Check whether the computation failed (in case learning rate is too large)
                                                   
      if(anyNA(coupl_gromov[[k]])==FALSE){
        #If succesful, save the Gromov-Wasserstein cost
        gw_cost[k] <- gromov_wstein_cost(coupl_gromov[[k]], Loss_tensor)
        k <- k + 1
      } else{
        #If failed, decrease learning rate
        learn_rate <- learn_rate - 0.02
        print(learn_rate)
      }
    }
    #Choose coupling that mimimizes Gromov-Wasserstein cost
    best_coupl <- coupl_gromov[[min(which(gw_cost==min(gw_cost)))]]
    #Create alignment from coupling
    algn_gromov <- coupling_to_alignNS(best_coupl)
    algn_meas <- get_align_measures(match1[[i]],match2[[i]],algn_gromov)
    so_seq[j] <- algn_meas$so
    rmsd_seq[j] <- algn_meas$rmsd
    nali_seq[j] <- algn_meas$nali
  }
  #Get optimal alpha
  opt_so <- which(so_seq ==max(so_seq))
  nali_seq_soopt <- nali_seq[opt_so]
  opt_index <- min(which(nali_seq ==max(nali_seq_soopt)))
  #Save alignment measures
  meas_matr_gromov[i,1] <- so_seq[opt_index]
  meas_matr_gromov[i,2] <- rmsd_seq[opt_index]
  meas_matr_gromov[i,3] <- nali_seq[opt_index]
  meas_matr_gromov[i,4] <- alpha_seq[opt_index]
}

#Quantized Gromov-Wasserstein
#Run on all protein pairs with at least one structure bigger than 80
for(i in 1:30){
  set.seed(2348)
  so_seq <- alpha_seq
  rmsd_seq <- alpha_seq
  nali_seq <- alpha_seq
  runtimes <- 15
  #Run for different alpha
  for (j in 1:length(alpha_seq)) {
    print(j)
    #Get quantized coupling
    coupl_quant <- quantized_partial_gromov(match1[[i]],match2[[i]],runtimes = 15,m_repr = 80,
                                            alpha = alpha_seq[j],tau = 0.03,iter_n = 200)
    #Get alignment from coupling
    algn_gromov_quant <- coupling_to_alignNS(coupl_quant)
    algn_meas <- get_align_measures(match1[[i]],match2[[i]],algn_gromov_quant)
    so_seq[j] <- algn_meas$so
    rmsd_seq[j] <- algn_meas$rmsd
    nali_seq[j] <- algn_meas$nali
  }
  #Get optimal alpha
  opt_so <- which(so_seq ==max(so_seq))
  nali_seq_soopt <- nali_seq[opt_so]
  opt_index <- min(which(nali_seq ==max(nali_seq_soopt)))
  #save alignment measures
  meas_matr_gromov[i,1] <- so_seq[opt_index]
  meas_matr_gromov[i,2] <- rmsd_seq[opt_index]
  meas_matr_gromov[i,3] <- nali_seq[opt_index]
  meas_matr_gromov[i,4] <- alpha_seq[opt_index]
  
  
}
#Results---------------------------------------------------------------------------------
#Monotone rearrangement median
meas_matr <- matrix(nrow = 30,ncol = 5)
meas_matr[,1] <- c(92.92,96.65,97,76.88,86.66,96.6,83.65,59.12,81.61,83.44,33.58,25.47,95,71.42,94.54,87.72,40.63,20,35.6,38.18,100,97.53,97.53,53.97,22.84,20.31,12.31,7.28,31.76,44.44)
meas_matr[,2] <- c(1.32,0.96,0.97,2.27,2.55,1.53,2.06,1.91,1.78,1.84,3.68,2.79,2.24,3.35,1.31,2.29,2.69,2.70,3.01,2.68,0.73,1.05,1.09,3.26,3.42,3.69,2.8,2.18,3.17,2.52)
meas_matr[,3] <- c(215,234,229,147,146,145,143,97,145,140,68,53,160,61,54,53,29,17,26,24,160,160,160,42,50,40,19,16,59,29)
#meas_matr[,4] <- c(4.68,3.15,34.08,1.16,3.83,1.54)
meas_matr[,5] <- c(1000,1000,1000,1000,2.64,1000,1000,78.46,1000,1000,33.6,14.38,1000,1000,1000,1000,33.6,33.6,6.16,1000,1000,1000,1000,6.16,78.47,1000,14.38,1000,2.64,1000)

#Monotone rearrangement partial Wasserstein
meas_matr_wstein <- matrix(nrow = 30,ncol = 6)
meas_matr_wstein[,1] <- c(92.92,92.65,97,77.45,87.3,96.57,83.64,69.74,81.6,83.5,67.1,77.6,95,72.72,94.54,92.11,82.8,27.5,91.5,100,100,97.3,98.5,61.9,29.01,27.3,26.15,23.79,33.11,32)
meas_matr_wstein[,2] <- c(1.33,2.34,0.97,2.28,2.12,1.533,2.289,2.9,2.2,1.89,2.47,2.41,2.24,2.37,1.31,2.52,2.16,3.12,2.02,1.16,0.73,1.04,1.08,3.48,2.66,2.79,2.47,2.78,3.75,2.73)
meas_matr_wstein[,3] <- c(215,234,229,147,142,145,146,113,154,141,121,137,160,59,54,54,56,26,58,55,160,160,160,49,53,38,36,54,63,32)
meas_matr_wstein[,5] <- c(1000,1000,1000,1.13,0.07,1000,0.016,0.0069,1000,2.63,0.089,0.038,1000,1000,1000,0.038,0.48,6.16,1000,1000,1000,1000,1000,1000,0.038,0.086,183.3,0.016,0.016,1000)
meas_matr_wstein[,6] <- c(1,1,1,1,0.5,1,0.7,0.6,0.6,0.9,0.7,0.8,1,1,1,0.8,0.9,1,0.7,1,1,1,1,0.8,0.6,0.6,1,0.6,0.6,0.5)
  
#(Quantized) Gromov-Wasserstein
meas_matr_gromov <- matrix(nrow = 30,ncol = 4)
meas_matr_gromov[,1] <- c(29.64,28.87,28.2,28.9,30.67,48.6,45.28,13.83,32.18,26.11,18.97,29.19,40.65,61.03,92.72,85.96,79.68,20,72.71,67.27,73.12,46.91,45.06,77.7,17.28,15.62,6.15,7.76,6.08,49.2)
meas_matr_gromov[,2] <- c(6.81,5.59,6.2,5.86,5.11,4.24,4.57,8.89,5.3,6.27,7.15,4.87,5.16,2.95,1.91,2.53,2.66,5.62,2.72,2.47,4.34,4.22,4.33,2.55,8.93,8.13,10.22,7.42,10.1,4.88)
meas_matr_gromov[,3] <- c(225,217,198,159,126,126,132,119,150,140,137,102,149,59,52,55,60,40,45,42,151,138,150,55,161,111,85,82,101,61)
meas_matr_gromov[,4] <- c(1,0.9,0.85,0.95,0.9,0.85,0.85,0.8,0.85,0.9,0.9,0.75,0.95,0.75,0.95,0.95,0.85,0.45,0.65,0.6,0.95,0.85,0.95,0.8,0.8,0.8,0.75,0.6,0.7,0.95)

#SPalign 
meas_matr_sp <- matrix(nrow = 30,ncol = 3)
meas_matr_sp[,1] <- c(92.04,96.65,96.58,76.3,86,95.8,83.02,59.12,79.31,82.8,76.64,76.4,95,72.7,
                      94.55,82.46,78.12,16.25,89.83,100,100,96.91,97.53,55.56,26.54,24.22,
                      20.77,24.76,33.11,49.21)
meas_matr_sp[,2] <- c(1.39,0.83,1.13,2.36,2.14,1.49,2.23,2.57,2.15,2.18,2.51,2.17,1.22,1.52,
                      1.32,2.49,2.12,3.34,1.74,1.16,0.73,1.05,1.09,3.25,3.92,3.89,3.1,3.99,
                      3.27,3.29)
meas_matr_sp[,3] <- c(216,233,231,148,142,144,145,110,153,144,121,135,155,57,54,53,56,19,57,
                      55,160,160,160,46,76,55,35,95,63,43)
#SPalign- NS
meas_matr_spns <- matrix(nrow = 30,ncol = 3)
meas_matr_spns[,1] <- c(94.25,96.65,96.15,76.88,88,96.58,84.91,60.38,91.61,86.62,78.83,81.37,
                        95,72.73,65.5,84.21,85.94,25,89.83,100,100,96.91,97.53,73.02,59.02,
                        48.44,33.85,48.54,39.19,71.43)
meas_matr_spns[,2] <- c(1.2,0.82,0.85,1.74,1.45,1.14,1.76,1.73,1.68,1.49,1.66,1.64,1.01,1.31,
                        0.99,1.94,1.46,2.13,1.27,1.16,0.73,0.9,0.96,2.16,2.12,2.19,2.02,2.24,
                        1.99,2.18)
meas_matr_spns[,3] <- c(213,231,225,133,132,141,135,96,142,136,108,131,152,56,52,48,55,20,53,
                        55,160,157,158,46,94,62,44,100,58,45)

#Plotting-----------------------------------------------------------------------------------
#Create tibbles for boxplot plotting
tib_homestrad <- tibble(method = c(rep("Median",23),rep("Wstein",23),rep("Gromov",23),
                            rep("SP",23),rep("SP-NS",23)),
                    so = c(meas_matr[1:23,1],meas_matr_wstein[1:23,1],
                        meas_matr_gromov[1:23,1],meas_matr_sp[1:23,1],meas_matr_spns[1:23,1]),
                    rmsd = c(meas_matr[1:23,2],meas_matr_wstein[1:23,2],
                        meas_matr_gromov[1:23,2],meas_matr_sp[1:23,2],meas_matr_spns[1:23,2]),
                    nali = c(meas_matr[1:23,3],meas_matr_wstein[1:23,3],
                        meas_matr_gromov[1:23,3],meas_matr_sp[1:23,3],meas_matr_spns[1:23,3]))
tib_ripc <- tibble(method = c(rep("Median",7),rep("Wstein",7),rep("Gromov",7),
                                   rep("SP",7),rep("SP-NS",7)),
                so = c(meas_matr[24:30,1],meas_matr_wstein[24:30,1],
                    meas_matr_gromov[24:30,1],meas_matr_sp[24:30,1],meas_matr_spns[24:30,1]),
                rmsd = c(meas_matr[24:30,2],meas_matr_wstein[24:30,2],
                    meas_matr_gromov[24:30,2],meas_matr_sp[24:30,2],meas_matr_spns[24:30,2]),
                nali = c(meas_matr[24:30,3],meas_matr_wstein[24:30,3],
                    meas_matr_gromov[24:30,3],meas_matr_sp[24:30,3],meas_matr_spns[24:30,3]))

#Homestred boxplots
#SO
basisplot(data = tib_homestrad,mapping = aes(x = method,y = so)) + 
  geom_boxplot(color="black",fill="blue",alpha=0.2,
               outlier.colour="black",outlier.fill="red", outlier.size=2,outlier.alpha = 0.4)+
  theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="darkgreen", fill="darkgreen") +
  ggtitle("Structure overlap (Homestrad)")
#RMSD
basisplot(data = tib_homestrad,mapping = aes(x = method,y = rmsd)) + 
  geom_boxplot(color="black",fill="red",alpha=0.2,
               outlier.colour="red",outlier.fill="red", outlier.size=2,outlier.alpha = 0.4)+
  theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="darkgreen", fill="darkgreen") +
  ggtitle("Rmsd (Homestrad)")
#Nali
basisplot(data = tib_homestrad,mapping = aes(x = method,y = nali)) + 
  geom_boxplot(color="black",fill="lightgreen",alpha=0.2,
               outlier.colour="red",outlier.fill="red", outlier.size=2,outlier.alpha = 0.4)+
  theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="darkgreen", fill="darkgreen") +
  ggtitle("Nali (Homestrad)")

#RIPC Boxplots
#SO
basisplot(data = tib_ripc,mapping = aes(x = method,y = so)) + 
  geom_boxplot(color="black",fill="blue",alpha=0.2,
               outlier.colour="red",outlier.fill="red", outlier.size=3)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="darkgreen", fill="darkgreen") +
  geom_jitter(color="black", size=1.2, alpha=0.9)+
  theme(legend.position="none") +
  ggtitle("Structure overlap (RIPC)")
#RMSD
basisplot(data = tib_ripc,mapping = aes(x = method,y = rmsd)) + 
    geom_boxplot(color="black",fill="red",alpha=0.2,
                 outlier.colour="red",outlier.fill="red", outlier.size=3)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="darkgreen", fill="darkgreen") +
    geom_jitter(color="black", size=1.2, alpha=0.9)+
  theme(legend.position="none") +
    ggtitle("Rmsd (RIPC)")
#Nali
basisplot(data = tib_ripc,mapping = aes(x = method,y = nali)) + 
    geom_boxplot(color="black",fill="lightgreen",alpha=0.2,
                 outlier.colour="red",outlier.fill="red", outlier.size=3)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="darkgreen", fill="darkgreen") +
    geom_jitter(color="black", size=1.2, alpha=0.9)+
  theme(legend.position="none") +
    ggtitle("Nali (RIPC)")

#Histograms of the distributions of lambda parameters for med. difference and part. Wasserstein
#Median difference
basisplot(tibble(x = log(meas_matr[,5])),aes(x = x)) + 
  geom_histogram(binwidth = 0.4,fill = "blue", alpha = 0.5) + xlab("log("~lambda~")")+
  ggtitle("Optimal "~ lambda ~ " (median difference)")
#partial Wasserstein
basisplot(tibble(x = log(meas_matr_wstein[,5])),aes(x = x)) + 
  geom_histogram(binwidth = 0.4,fill = "red", alpha = 0.5) + xlab("log("~lambda~")")+
  ggtitle("Optimal "~ lambda ~ " (partial Wasserstein)")

#Histograms of the distributions of alpha parameters for both datasets
#Homestrad
basisplot(tibble(x = meas_matr_wstein[1:23,6]),aes(x = x)) + 
  geom_histogram(binwidth = 0.05,fill = "blue", alpha = 0.5) + xlab(expression(alpha))+
  ggtitle("Optimal partial Wstein"~ alpha ~ " (Homestrad)")
#RIPC
basisplot(tibble(x = meas_matr_wstein[24:30,6]),aes(x = x)) + 
  geom_histogram(binwidth = 0.05,fill = "red", alpha = 0.5) + xlab(expression(alpha))+
  ggtitle("Optimal partial Wstein"~ alpha ~ " (RIPC)")

#Experiment: Applying quantized Gromov-Wasserstein with different numbers of representatives------
#Fix protein pair
p1 <- match1[[15]]
p2 <- match2[[15]]
#List of numbers of representatives
repr_seq <- seq(5,55,5)
#Matrix to save alignment measures
meas_matr_quant <- matrix(nrow = length(repr_seq),ncol = 4)
#List of alpha parameters
alpha_seq <- seq(0.6,1,0.05)

#Run for all numbers of representatives
for(i in 1:length(repr_seq)){
  set.seed(2348)
  so_seq <- alpha_seq
  rmsd_seq <- alpha_seq
  nali_seq <- alpha_seq
  runtimes <- 15
  #Run for different alpha
  for (j in 1:length(alpha_seq)) {
    print(j)
    #Create coupling
    coupl_quant <- quantized_partial_gromov(p1,p2,runtimes = 15,m_repr = repr_seq[i],
                                            alpha = alpha_seq[j],tau = 0.03,iter_n = 200)
    #Create alignment
    algn_gromov_quant <- coupling_to_alignNS(coupl_quant)
    algn_meas <- get_align_measures(p1,p2,algn_gromov_quant)
    so_seq[j] <- algn_meas$so
    rmsd_seq[j] <- algn_meas$rmsd
    nali_seq[j] <- algn_meas$nali
  }
  #Get optimal alpha
  opt_so <- which(so_seq ==max(so_seq))
  nali_seq_soopt <- nali_seq[opt_so]
  opt_index <- min(which(nali_seq ==max(nali_seq_soopt)))
  #Save alignment measures
  meas_matr_quant[i,1] <- so_seq[opt_index]
  meas_matr_quant[i,2] <- rmsd_seq[opt_index]
  meas_matr_quant[i,3] <- nali_seq[opt_index]
  meas_matr_quant[i,4] <- alpha_seq[opt_index]
}

#Plot alignment measures over number of representatives
basisplot() + 
  geom_point(stat="identity",size = 3,data = tibble(x = repr_seq,y = meas_matr_quant[,1]),mapping = aes(x = x,y = y,color = "blue"))+
  geom_point(stat="identity",size = 3,data = tibble(x = repr_seq,y = meas_matr_quant[,2]),mapping = aes(x = x,y = y,color = "red"))+
  geom_point(stat="identity",size = 3,data = tibble(x = repr_seq,y = meas_matr_quant[,3]),mapping = aes(x = x,y = y,color= "green"))+
  ylab("value")+  xlab("Number of representatives")+
  ggtitle("Align. meas. over number of representatives")+
  scale_color_identity(name = "Align. measure",
                       breaks = c("blue","red","green"),
                       labels = c("SO","RMSD","Nali"),
                       guide = "legend")

#Example where Gromov-Wasserstein fails----------------------------------------------------
#Fix protein pair
p1 <- match1[[20]]
p2<- match2[[20]]

set.seed(31434)
coupl_gromov <- list()
gw_cost <- vector(length = 15)
#Compute loss tensor
Loss_tensor <- get_loss_tensor(p1,p2)
#Run many times to avoid local minima
k <- 1
while(k <= 15){
  #Get coupling
  coupl_gromov[[k]] <- partial_gromov_entropic(pr1 =p1,pr2 =p2,
                                               alpha = 1,iter_n = 200,
                                               Loss_tensor = Loss_tensor,tau = 0.3)
  #Get Gromov-Wasserstein cost
  gw_cost[k] <- gromov_wstein_cost(coupl_gromov[[k]], Loss_tensor)
  k <- k + 1
}
#Choose coupling minimizing the Gromov-Wasserstein cost
best_coupl <- coupl_gromov[[min(which(gw_cost==min(gw_cost)))]]
#Get coupling from alignment
align_gromov <- coupling_to_alignNS(best_coupl)
get_align_measures(p1,p2,align_gromov)
#3-D plot of alignment
plot_algn_3D(p1,p2,align_gromov)

#Monotone rearrangement with partial Wasserstein distance (works)
align_wstein <- align_proteins(p1,p2,return_meas = FALSE,lambda = 1000,iter = 10,
                               dist = "partial_wasserstein",print_info = FALSE,alpha = 1)
get_align_measures(p1,p2,align_wstein)
plot_algn_3D(match1[[20]],match2[[20]],align = align_wstein)

#Gromov-Wasserstein approach works if we cut of remaining residues
set.seed(32634)
coupl_gromov_cut <- list()
gw_cost_cut <- vector(length = 15)
#Loss tensor when cutting off
Loss_tensor_cut <- get_loss_tensor(p1,p2[1:57,])
#Run many times to avoid local minima
k <- 1
while(k <= 15){
  #Compute coupling
  #Choose a larger number than 55 to investigate the behaviour when cutting less residues 
  coupl_gromov_cut[[k]] <- partial_gromov_entropic(pr1 =p1,pr2 =p2[1:57,],
                                               alpha = 1,iter_n = 200,
                                               Loss_tensor = Loss_tensor_cut,tau = 0.3)
  #Get Gromov-Wasserstein cost
  gw_cost_cut[k] <- gromov_wstein_cost(coupl_gromov_cut[[k]], Loss_tensor_cut)
  k <- k + 1

}
#Coupling minimizing Gromov-Wasserstein cost
best_coupl_cut <- coupl_gromov_cut[[min(which(gw_cost_cut==min(gw_cost_cut)))]]
#Get alignment from coupling
align_gromov_cut <- coupling_to_alignNS(best_coupl_cut)
get_align_measures(p1,p2,align_gromov_cut)
#Plot coupling
plot_wstein_matrix(wstein_matrix = best_coupl_cut, title = "Gromov-Wasserstein coupling (1)")
#3-D plot of alignment
plot_algn_3D(p1,p2,align_gromov_cut)

#Set alpha = 55/106 and compare Gromov-coupling with actual optimal coupling
#First: rerun above code with alpha = 55/106 to obain best_coupl
#Now run code below to obtain the actual best coupling (best_coupl_2)
part1 <- diag(1/106,ncol = 8, nrow = 8)
zeros_part1 <- matrix(0,nrow = 8,ncol = 98)
part1 <- cbind(part1,zeros_part1)
part2 <- diag(1/106,ncol = 47, nrow = 47)
part2 <- cbind(matrix(0,nrow=47,ncol=10),part2,matrix(0,nrow=47,ncol=49))
best_coupl2 <- rbind(part1,part2)

#Comparing Gromov-Wasserstein costs -> best_coupl2 ist better
gromov_wstein_cost(best_coupl,Loss_tensor)
gromov_wstein_cost(best_coupl2,Loss_tensor)


#Experiment: Perturbing the residues of a fixed proteins with different variance-------------
#Fix structure 1
p_sim1 <- match1[[15]]
set.seed(347284)
#List of standard deviations
sigma_list <- seq(0,10,0.5)
#Get list of structures by perturbing structure 1 with sigma
p_perturb <- list()
for(i in 1:length(sigma_list)){
  p_perturb[[i]] <- perturb_protein(p_sim1,sigma_list[i])
}
#Write structures to file to apply SP-align (NS)
write_protein_to_pdb(p_sim1,"p_sim1.pdb")
write_protein_to_pdb(p_perturb[[21]],"p_sim2.pdb")

#Apply different alignment measures
#
#Monotone rearrangement with median difference
lambda_seq <- exp(seq(log(0.0001), log(1000), length.out = 20))
meas_matr <- matrix(nrow = 21,ncol = 4)
#Run over all perturbed structures
for (i in 1:21) {
  print(i)
  print(Sys.time())
  so_vec <- vector(length = length(lambda_seq))
  #Run for different lambda
  for(l in 1:length(lambda_seq)){
    print(l)
    alignment <- align_proteins(p_sim1,p_perturb[[i]],return_meas = TRUE,lambda = lambda_seq[l],iter = 10,
                                dist = "median_diff",print_info = FALSE)
    so_vec[l] <- alignment$measures$so[length(alignment$measures$so)]
  }
  #Optimal lambda
  lambda_opt <- max(lambda_seq[which(so_vec==max(so_vec))])
  print(lambda_opt)
  alignment <- align_proteins(p_sim1,p_perturb[[i]],return_meas = TRUE,lambda = lambda_opt,iter = 10,
                              dist = "median_diff",print_info = FALSE)
  #Save alignment measures
  meas_matr[i,1] <- alignment$measures$so[length(alignment$measures$so)]
  meas_matr[i,2] <- alignment$measures$rmsd[length(alignment$measures$rmsd)]
  meas_matr[i,3] <- alignment$measures$nali[length(alignment$measures$nali)]
  meas_matr[i,4] <- lambda_opt
}

#Monotone rearrangement with partial Wasserstein distance
lambda_seq <- exp(seq(log(0.0001), log(1000), length.out = 20))
alpha_seq <- seq(0.5,1,0.1)
meas_matr_wstein <- matrix(nrow = 21,ncol = 5)
#Run over all perturbed structures
for (i in 1:21) {
  print(i)
  print(Sys.time())
  so_matr <- matrix(nrow = length(alpha_seq),ncol = length(lambda_seq))
  nali_matr <- matrix(nrow = length(alpha_seq),ncol = length(lambda_seq))
  ird1 <- get_inter_res_dist(p_sim1)
  ird2 <- get_inter_res_dist(p_perturb[[i]])
  #Run for different alpha
  for(k in 1:length(alpha_seq)){
    print(paste("Alpha =",alpha_seq[k]))
    #Compute Wasserstein matrix for alpha
    wstein_matr <- wstein_dist_matr(ird1,ird2,p = 1,alpha = alpha_seq[k],dist="partial_wasserstein")
    #Run for different lambda
    for(l in 1:length(lambda_seq)){
      print(paste("Lambda =",lambda_seq[l]))
      #Alignment
      alignment <- align_proteins(p_sim1,p_perturb[[i]],return_meas = TRUE,lambda = lambda_seq[l],
                                  iter = 10, wstein_matr = wstein_matr,
                                  dist = "partial_wasserstein",print_info = FALSE,alpha = alpha_seq[k])
      so_matr[k,l] <- alignment$measures$so[length(alignment$measures$so)]
      nali_matr[k,l] <- alignment$measures$nali[length(alignment$measures$nali)]
    }
  }
  #Optimal Parameters
  opt_so <- which(so_matr == max(so_matr), arr.ind = TRUE)
  nali_vec <- nali_matr[opt_so]
  opt_nali <- max(which(nali_vec == max(nali_vec)))
  opt_par <- opt_so[opt_nali,]
  alpha_opt <- alpha_seq[opt_par[1]]
  lambda_opt <- lambda_seq[opt_par[2]]
  print(paste(alpha_opt,lambda_opt))
  #Optimal alignment
  alignment <- align_proteins(p_sim1,p_perturb[[i]],return_meas = TRUE,lambda = lambda_opt,iter = 10,
                              dist = "partial_wasserstein",print_info = FALSE,alpha = alpha_opt)
  #Save alignment measures
  meas_matr_wstein[i,1] <- alignment$measures$so[length(alignment$measures$so)]
  meas_matr_wstein[i,2] <- alignment$measures$rmsd[length(alignment$measures$rmsd)]
  meas_matr_wstein[i,3] <- alignment$measures$nali[length(alignment$measures$nali)]
  meas_matr_wstein[i,4] <- lambda_opt
  meas_matr_wstein[i,5] <- alpha_opt
}

#Gromov Wasserstein
meas_matr_gromov <- matrix(nrow = 21,ncol = 4)
alpha_seq <- seq(0.3,0.8,0.05)
##Run over all perturbed structures
for(i in 15:21){
  print(i)
  set.seed(2348)
  #Get loss tensor for protein pair
  Loss_tensor <- get_loss_tensor(p_sim1,p_perturb[[i]])
  so_seq <- alpha_seq
  rmsd_seq <- alpha_seq
  nali_seq <- alpha_seq
  runtimes <- 15
  #Run for different alpha
  for (j in 1:length(alpha_seq)) {
    print(j)
    coupl_gromov <- list()
    gw_cost <- vector(length = runtimes)
    learn_rate <- 0.2
    #Run many times to avoid local minima
    k <- 1
    while(k <= runtimes){
      coupl_gromov[[k]] <- partial_gromov_entropic(pr1 =p_sim1,pr2 =p_perturb[[i]],
                                                   alpha = alpha_seq[j],iter_n = 200,
                                                   Loss_tensor = Loss_tensor,tau = learn_rate)
      #Check whether computation failed (learning rate too large)
      if(anyNA(coupl_gromov[[k]])==FALSE){
        gw_cost[k] <- gromov_wstein_cost(coupl_gromov[[k]], Loss_tensor)
        print(gw_cost[k])
        k <- k + 1
      } else{
        learn_rate <- learn_rate - 0.02
        print(learn_rate)
      }
    }
    #Coupling minimizing Gromov-Wasserstein cost
    best_coupl <- coupl_gromov[[min(which(gw_cost==min(gw_cost)))]]
    #Get alignment from coupling
    algn_gromov <- coupling_to_alignNS(best_coupl)
    algn_meas <- get_align_measures(p_sim1,p_perturb[[i]],algn_gromov)
    so_seq[j] <- algn_meas$so
    rmsd_seq[j] <- algn_meas$rmsd
    nali_seq[j] <- algn_meas$nali
    print(so_seq[j])
  }
  #Get optimal alpha
  opt_so <- which(so_seq ==max(so_seq))
  nali_seq_soopt <- nali_seq[opt_so]
  opt_index <- min(which(nali_seq ==max(nali_seq_soopt)))
  #Save alignment measures
  meas_matr_gromov[i,1] <- so_seq[opt_index]
  meas_matr_gromov[i,2] <- rmsd_seq[opt_index]
  meas_matr_gromov[i,3] <- nali_seq[opt_index]
  meas_matr_gromov[i,4] <- alpha_seq[opt_index]
}

#Results---------------------------------------------------------------------------------
#Median difference monotone rearrangement
meas_matr <- matrix(nrow = 21,ncol = 4)
meas_matr[,1] <- c(100,100,100,89.09,61.82,54.54,54.54,49.1,32.73,25.45,30.9,27.27,25.45,18.18,21.82,20,25.45,14.54,16.36,16.36,16.36)
meas_matr[,2] <- c(0,0.85,1.69,2.46,3.2,3.06,2.7,2.96,2.44,4.96,3.87,2.81,2.96,4.45,2.57,2.88,2.92,2.48,3.33,3.02,3.08)
meas_matr[,3] <- c(55,55,55,53,47,39,35,33,19,38,36,17,16,22,13,13,16,8,12,10,12)
meas_matr[,4] <- c(1000,1000,1000,1000,1000,14.38,1000,1000,1000,1000,6.16,1000,1000,6.16,1000,1000,1000,1000,189.3,33.6,2.64)
#partial Wasserstein distance
meas_matr_wstein <- matrix(nrow = 21,ncol = 5)
meas_matr_wstein[,1] <- c(100,100,100,89.9,63.6,52.7,56.4,49.1,32.7,30.91,30.91,29.1,25.45,23.64,23.64,21.8,25.45,18.18,21.82,18.18,18.18)
meas_matr_wstein[,2] <- c(0,0.85,1.69,2.5,2.9,3.18,2.88,2.96,2.96,3.77,2.96,3.73,3.36,4.43,2.57,2.29,5.22,4.26,2.81,2.65,2.97)
meas_matr_wstein[,3] <- c(55,55,55,53,43,41,36,33,23,27,20,22,17,24,14,12,26,15,13,10,12)
meas_matr_wstein[,4] <- c(1000,1000,1000,1000,1000,0.48,0.089,1000,1000,2.64,1000,0.089,6.16,1000,1000,33.6,2.64,1000,1.13,0.21,14.38)
meas_matr_wstein[,5] <- c(1,1,1,1,1,0.5,0.5,0.8,1,0.8,0.8,0.6,0.8,0.9,1,0.9,0.8,0.5,0.8,0.5,0.8)
#Gromov-Wasserstein
meas_matr_gromov <- matrix(nrow = 21,ncol = 4)
meas_matr_gromov[,1] <- c(100,100,100,83.63,65.45,58.18,60,61.81,52.72,45.45,41.82,49.1,36.36,41.82,40,36.36,41.82,29.09,30.91,23.64,14.54)
meas_matr_gromov[,2] <- c(0,0.85,1.69,2.2,2.62,2.65,2.65,2.96,2.45,3.47,3.96,3.14,4.17,2.54,2.4,2.77,4.22,2.9,2.42,2.42,7.48)
meas_matr_gromov[,3] <- c(55,55,55,49,41,38,38,44,32,36,41,38,36,25,24,25,33,22,18,15,18)
meas_matr_gromov[,4] <- c(1,1,1,0.9,0.75,0.7,0.7,0.8,0.6,0.65,0.75,0.7,0.65,0.45,0.45,0.45,0.6,0.4,0.35,0.3,0.35)
#SPalign
meas_matr_sp <- matrix(nrow = 21,ncol = 4)
meas_matr_sp[,1] <- c(100,100,100,83.64,52.73,47.27,45.45,41.82,29.09,21.82,23.64,21.82,20,21.82,18.18,16.36,18.18,12.73,16.36,16.36,14.55)
meas_matr_sp[,2] <- c(0,0.85,1.69,2.65,3.59,3.6,3.75,3.69,4.56,5.02,3.97,4.37,4.21,4.49,4.25,4.26,4.04,4.06,3.99,3.1,3.68)
meas_matr_sp[,3] <- c(55,55,55,55,52,46,45,48,35,38,28,26,26,24,21,16,21,16,15,13,12)
#SPalign NS
meas_matr_spns <- matrix(nrow = 21,ncol = 4)
meas_matr_spns[,1] <- c(100,100,100,90.91,56.36,61.82,56.36,63.64,43.64,40,43.64,49.09,41.82,34.55,45.45,38.18,36.36,29.09,29.09,21.82,29.09)
meas_matr_spns[,2] <- c(0,0.85,1.69,2.2,2.27,2.31,2.24,2.29,2.34,2.29,2.24,2.28,2.43,2.2,2.59,2.29,2.26,1.83,2.56,1.88,2.38)
meas_matr_spns[,3] <- c(55,55,55,50,31,34,31,35,24,22,24,27,23,19,25,21,20,16,16,12,16)

#Plot alignment measures over sigma
#SO
basisplot() + 
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr[,1]),mapping = aes(x = x,y = y,color = "blue"))+
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr_wstein[,1]),mapping = aes(x = x,y = y,color = "red"))+
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr_gromov[,1]),mapping = aes(x = x,y = y,color= "green"))+
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr_sp[,1]),mapping = aes(x = x,y = y,color= "purple"))+
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr_spns[,1]),mapping = aes(x = x,y = y,color= "black"))+
  ylab("SO")+  xlab(expression(sigma))+
  ggtitle("Structure overlap")+
  scale_color_identity(name = "Align. measure",
                       breaks = c("blue","red","green","purple","black"),
                       labels = c("Median","Wstein","Gromov","SP","SP-NS"),
                       guide = "legend")
#RMSD
basisplot() + 
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr[,2]),mapping = aes(x = x,y = y,color = "blue"))+
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr_wstein[,2]),mapping = aes(x = x,y = y,color = "red"))+
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr_gromov[,2]),mapping = aes(x = x,y = y,color= "green"))+
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr_sp[,2]),mapping = aes(x = x,y = y,color= "purple"))+
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr_spns[,2]),mapping = aes(x = x,y = y,color= "black"))+
  ylab("RMSD")+  xlab(expression(sigma))+
  ggtitle("RMSD")+
  scale_color_identity(name = "Align. measure",
                       breaks = c("blue","red","green","purple","black"),
                       labels = c("Median","Wstein","Gromov","SP","SP-NS"),
                       guide = "legend")
#Nali
basisplot() + 
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr[,3]),mapping = aes(x = x,y = y,color = "blue"))+
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr_wstein[,3]),mapping = aes(x = x,y = y,color = "red"))+
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr_gromov[,3]),mapping = aes(x = x,y = y,color= "green"))+
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr_sp[,3]),mapping = aes(x = x,y = y,color= "purple"))+
  geom_line(stat="identity",size = 1,data = tibble(x = sigma_list,y = meas_matr_spns[,3]),mapping = aes(x = x,y = y,color= "black"))+
  ylab("Nali")+  xlab(expression(sigma))+
  ggtitle("Number of aligned residues")+
  scale_color_identity(name = "Align. measure",
                       breaks = c("blue","red","green","purple","black"),
                       labels = c("Median","Wstein","Gromov","SP","SP-NS"),
                       guide = "legend")



