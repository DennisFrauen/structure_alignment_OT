#Code for all plots in chapter 4
library(bio3d)
source("structure_allignment.R")
source("gromov_wstein.R")
#Import the fixed protein pair of size 63/77 from the RIPC dataset
pdb1 <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d1biaa1.pdb",ATOM.only = TRUE)
pdb2 <-  bio3d::read.pdb(file = "./protein_structures/RIPC/d1qdma1.pdb",ATOM.only = TRUE)
p1 <- pdb1$atom %>% select(x,y,z) %>% slice(which(pdb1$atom$elety=="CA")) %>% as.matrix()
p2 <- pdb2$atom %>% select(x,y,z) %>% slice(which(pdb2$atom$elety=="CA")) %>% as.matrix()

#Create a basisplot for later plots
basisplot <- ggplot() + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=17),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=17))

#Gromov-Wasserstein: Plot of Alpha against kappa----------------------------------------
set.seed(32490324)
#list of alphas
alpha_seq <- seq(0.1,1,0.1)
#list of kappas
nali_seq <- vector(length = length(alpha_seq))
#Compute Loss tensor
Loss_tensor <- get_loss_tensor(p1,p2)
#Run over different alphas
for(i in 1:length(alpha_seq)){
  #Compute Gromov-Wasserstein coupling
  coup_entr <- partial_gromov_entropic(pr1 = p1,pr2 = p2,alpha = 1,
                                       iter_n = 1000,iter_k = 100,tau = 0.15,Loss_tensor = Loss_tensor)
  #Get alignment from coupling
  align <- coupling_to_alignNS(coup_entr)
  #Save resulting kappa
  nali_seq[i] <- length(align$I)
  print(nali_seq[i])
}

#Plot
basisplot + 
  geom_line(stat="identity",size = 1, data = tibble(x= alpha_seq,y = nali_seq),
            mapping = aes(x = x,y = y),color = "blue")+
  ylab(expression(kappa))+  xlab(expression(alpha))+
  ggtitle("Plot of "~ kappa ~ " against" ~alpha)

#Monotone rearrangement: Alignment measures for different choices of kappa--------------------
s = min(dim(p1)[1],dim(p2)[1])
#Compute Wasserstein distance matrix
ird1 <- get_inter_res_dist(p1)
ird2 <- get_inter_res_dist(p2)
wstein_matr <- wstein_dist_matr(ird1,ird2,p = 1,dist="wasserstein")
#Compute monotone rearrangement tensor
H <- monotone_rearr(wstein_matr)
#Lists for saving the alignment measures
so_vec <- vector(length = s)
rmsd_vec <- vector(length = s)
so_vec[1] <- 100/63
#Run over all possible kappa
for(k in 3:(s+1)){
  #Get alignment by traceback procedure for kappa
  algn <- get_optimal_align(H[,,1:k])
  print(get_align_measures(p1,p2,algn,"nali",TRUE))
  #Save alignment measures
  so_vec[k-1] <- get_align_measures(p1,p2,algn,"so",TRUE)
  rmsd_vec[k-1] <- get_align_measures(p1,p2,algn,"rmsd",TRUE)
}
#Plot alignment measures over kappa
basisplot + 
  geom_line(stat="identity",size = 1, data = tibble(x= 1:s,y = so_vec),
            mapping = aes(x = x,y = y,color = "blue"))+
  geom_line(stat="identity",size = 1, data = tibble(x= 1:s,y = rmsd_vec),
            mapping = aes(x = x,y = y,color = "red"))+
  ylab("value")+  xlab(expression(kappa))+
  ggtitle("Plot of SO and RMSD over "~ kappa)+
  scale_color_identity(name = "Align. measure",
                       breaks = c("blue","red"),
                       labels = c("SO","RMSD"),
                       guide = "legend")

#Partial Wstein monotone rearrangement: Alignment measures over different alpha-------------------------------------
#List of possible alpha
alpha_seq <- seq(0.2,1,0.02)
#empty vectors to save alignment measures
so_seq <- vector(length = length(alpha_seq))
rmsd_seq <- vector(length = length(alpha_seq))
nali_seq <- vector(length = length(alpha_seq))

#Run over different alpha
for (i in 1:length(alpha_seq)) {
  #Get monotone rearrangement alignment using the alpha-partial Wasserstein distance
  alignment <- align_proteins(p1 = p1,p2 = p2,alpha = alpha_seq[i],dist = "partial_wasserstein")
  #Save alignment measures
  so_seq[[i]] <- get_align_measures(p1,p2,alignment,type = "so")
  rmsd_seq[[i]] <- get_align_measures(p1,p2,alignment,type = "rmsd")
  nali_seq[[i]] <- get_align_measures(p1,p2,alignment,type = "nali")
}

#Plot alignment measures over alpha
#SO + Nali
basisplot + geom_line(stat="identity",size = 1, data = tibble(x= alpha_seq,y = so_seq),
            mapping = aes(x = x,y = y,color = "blue"))+
  geom_line(stat="identity",size = 1, data = tibble(x= alpha_seq,y = nali_seq),
            mapping = aes(x = x,y = y,color = "green"))+
  ylab("value")+  xlab(expression(alpha))+
  ggtitle("SO and Nali over "~ alpha)+
  scale_color_identity(name = "Align. measure",
                       breaks = c("blue","red","green"),
                       labels = c("SO","RMSD","Nali"),
                       guide = "legend")
#RMSD
basisplot + geom_line(stat="identity",size = 1, data = tibble(x= alpha_seq,y = rmsd_seq),
                      mapping = aes(x = x,y = y,color = "red"))+
  ylab("value")+  xlab(expression(alpha))+
  ggtitle("RMSD over "~ alpha)+
  scale_color_identity(name = "Align. measure",
                       breaks = c("red"),
                       labels = c("RMSD"),
                       guide = "legend")

#Iterative monotone rearrangement-----------------------------------------------------------
#Alignment measures over different lambda
lambda_seq <- exp(seq(log(0.0001), log(1000), length.out = 100))
so_vec <- vector(length = length(lambda_seq))
rmsd_vec <- vector(length = length(lambda_seq))
nali_vec <- vector(length = length(lambda_seq))
#Run for different lambda
for (i in 1:length(lambda_seq)) {
  #Apply iterative monotone rearrangement with lambda and 10 iterations
  alignment <- align_proteins(p1,p2,return_meas = TRUE,lambda = lambda_seq[i],iter = 10,
                              dist = "median_diff",alpha = 0.6,remove_res = TRUE)
  so_vec[i] <- alignment$measures$so[length(alignment$measures$so)]
  #Save alignment measures
  rmsd_vec[i] <- alignment$measures$rmsd[length(alignment$measures$rmsd)]
  nali_vec[i] <- alignment$measures$nali[length(alignment$measures$nali)]
}

#Plot alignment measures over lambda
basisplot + geom_line(stat="identity",size = 1, data = tibble(x= log(lambda_seq),y = log(so_vec)),
                      mapping = aes(x = x,y = y,color = "blue"))+
  geom_line(stat="identity",size = 1, data = tibble(x= log(lambda_seq),y = log(rmsd_vec)),
            mapping = aes(x = x,y = y,color = "red"))+
  geom_line(stat="identity",size = 1, data = tibble(x= log(lambda_seq),y = log(nali_vec)),
            mapping = aes(x = x,y = y,color = "green"))+
  ylab("value (log)")+  xlab(lambda ~ "(log)")+
  ggtitle("Align. meas. over "~ lambda ~ "(median)")+
  scale_color_identity(name = "Align. measure",
                       breaks = c("blue","red","green"),
                       labels = c("SO","RMSD","Nali"),
                       guide = "legend")

#Get Alignment measures for each iteration step (using the optimal lambda)
#Optimal lambda alignment
lambda_opt <- max(lambda_seq[which(so_vec==max(so_vec))])
alignment <- align_proteins(p1,p2,return_meas = TRUE,lambda = lambda_opt,iter = 10,
                            dist = "median_dist",alpha = 0.6)
#alignment measures for each iteration
align_meas_iter <- alignment$measures
iter <-1:length(align_meas_iter$so)

#Plot alignment measures over iterations
basisplot + geom_point(stat="identity",size = 2, data = tibble(x= iter,y = log(align_meas_iter$so)),
                      mapping = aes(x = x,y = y,color = "blue"))+
  geom_point(stat="identity",size = 2, data = tibble(x= iter,y = log(align_meas_iter$rmsd)),
            mapping = aes(x = x,y = y,color = "red"))+
  geom_point(stat="identity",size = 2, data = tibble(x= iter,y = log(align_meas_iter$nali)),
            mapping = aes(x = x,y = y,color = "green"))+
  ylab("value (log)")+  xlab("Iteration")+
  ggtitle("Align. meas. over iterations (median")+
  scale_color_identity(name = "Align. measure",
                       breaks = c("blue","red","green"),
                       labels = c("SO","RMSD","Nali"),
                       guide = "legend")

