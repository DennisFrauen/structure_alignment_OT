#Code for artificial examples discussed in section 5.2

library(bio3d)
source("structure_allignment.R")
source("gromov_wstein.R")
source("sim_fcn.R")


#First artificial example--------------------------------------------------------------------
#Create structure 1 (hexagon)
create_art_structure1 <- function(radius){
  r1 <- c(radius,0,0)
  r2 <- c(radius/2,sqrt(3)*radius/2,0)
  r3 <- c(-radius/2,sqrt(3)*radius/2,0)
  r4 <- c(-radius,0,0)
  r5 <- c(-radius/2,-sqrt(3)*radius/2,0)
  r6 <- c(radius/2,-sqrt(3)*radius/2,0)
  return(t(matrix(c(r1,r2,r3,r4,r5,r6),nrow = 3)))
}
#Create structure 2 (modified hexagon)
create_art_structure2 <- function(radius){
  r1 <- c(radius,sqrt(3)*radius,0)
  r2 <- c(radius,0,0)
  r3 <- c(radius/2,sqrt(3)*radius/2,0)
  r4 <- c(-radius/2,sqrt(3)*radius/2,0)
  r5 <- c(-radius,0,0)
  r6 <- c(radius,sqrt(3)*radius + radius,0)
  return(t(matrix(c(r1,r2,r3,r4,r5,r6),nrow = 3)))
}
p_sim1 <- create_art_structure1(14)
p_sim2 <- create_art_structure2(14)
#Write to pdb file, so SP align can be applied
write_protein_to_pdb(p_sim1,"p_sim1.pdb")
write_protein_to_pdb(p_sim2,"p_sim2.pdb")
#Plot both proteins
plot_algn_3D(p_sim1,superpos = FALSE,plt_algn = FALSE,plt_connections = TRUE,align = list(I = 1:6,J = 1:6),
             x_cap = "",y_cap = "",z_cap = "",startcol1 = "green",size = 0.5)
plot_algn_3D(p_sim2,superpos = FALSE,plt_algn = FALSE,plt_connections = TRUE,align = list(I = 1:6,J = 1:6),
             x_cap = "",y_cap = "",z_cap = "",col1 = "red",startcol1 = "green",size=0.5,cex = 2)
#Perfect alignment/ SP align/ SP-NS
align_perfect <- list(I = 1:4, J = 2:5)
get_align_measures(p_sim1,p_sim2,align_perfect)
plot_algn_3D(p_sim1,p_sim2,align = align_perfect,
      superpos = TRUE,plt_algn = TRUE, plt_connections = TRUE,x_cap = "",y_cap = "",z_cap = "",
       size = 0.6,startcol1 = "green",startcol2 = "green")

#Gromov-Wasserstein
Loss_tensor <- get_loss_tensor(p_sim1,p_sim2)
coupl_gromov <- matrix(0,nrow = 6,ncol = 6)
coupl_gromov[1,4] <- 1/6
coupl_gromov[2,5] <- 1/6
coupl_gromov[5,2] <- 1/6
coupl_gromov[6,3] <- 1/6
gromov_wstein_cost(coupl_gromov,Loss_tensor)
align_gromov <- coupling_to_alignNS(coupl_gromov)
get_align_measures(p_sim1,p_sim2,align_gromov)
plot_algn_3D(p_sim1,p_sim2,align = align_gromov,
             superpos = TRUE,plt_algn = TRUE, plt_connections = TRUE,x_cap = "",y_cap = "",z_cap = "",
             size = 0.6,startcol1 = "green",startcol2 = "green")

#Monotone rearrangement (Median) (doesnt get perfect alignment)
align_med <- align_proteins(p_sim1,p_sim2,return_meas = FALSE,lambda = 1,iter = 1,
                              dist = "median_diff",print_info = FALSE)
get_align_measures(p_sim1,p_sim2,align_med)



#Wasserstein alignment (doesnt get perfect alignment - the same as for median difference)
align_wstein <- align_proteins(p_sim1,p_sim2,return_meas = FALSE,lambda = 1,
                            iter = 1,
                            dist = "partial_wasserstein",print_info = FALSE,alpha = 1)

get_align_measures(p_sim1,p_sim2,align_wstein)
#Plot
plot_algn_3D(p_sim1,p_sim2,align = align_wstein,
             superpos = TRUE,plt_algn = TRUE, plt_connections = TRUE,x_cap = "",y_cap = "",z_cap = "",
             size = 0.6,startcol1 = "green",startcol2 = "green")


#IRD plotting
ird1 = get_inter_res_dist(p_sim1)
ird2 = get_inter_res_dist(p_sim2)
wmatr <- wstein_dist_matr(ird1,ird2,alpha = 1,dist = "partial_wasserstein")
plot_wstein_matrix(wmatr,plot_align1 = TRUE,align1 = align_wstein,
                   title = "Wasserstein distance matrix (1)")
#Alpha = 5/6
wmatr <- wstein_dist_matr(ird1,ird2,alpha = 5/6,dist = "partial_wasserstein")
plot_wstein_matrix(wmatr,plot_align1 = TRUE,align1 = list(I = 1:4,J = 1:4),
                   title = "Wasserstein distance matrix (5/6)")
#Alpha = 4/6
wmatr <- wstein_dist_matr(ird1,ird2,alpha = 4/6,dist = "partial_wasserstein")
plot_wstein_matrix(wmatr,plot_align1 = TRUE,align1 = list(I = 1:4,J = 1:4),
                   title = "Wasserstein distance matrix (4/6)")

#Median
wmatr <- wstein_dist_matr(ird1,ird2,alpha = 5/6,dist = "median_diff")
plot_wstein_matrix(wmatr,plot_align1 = TRUE,align1 = align_med,
                   title = "Median distance matrix")


#Artificial Gromov-Wasserstein example-------------------------------------------------------

create_shape1 <- function(size){
  r1 <- c(0,0,0)
  r2 <- c(size,0,size*2)
  r3 <- c(size,0,0)
  r4 <- c(size,4*size,0)
  r5 <- c(size,4*size,-size)
  r6 <- c(size,-size,0)
  return(t(matrix(c(r1,r2,r3,r4,r5),nrow = 3)))
}

p_sim1 <- create_shape1(4)
mirror_matr <- matrix(nrow = 3,ncol = 3,data = c(1,0,0,0,1,0,0,0,-1))
p_sim2 <- p_sim1 %*% mirror_matr

#Write to pdb file, so SP align can be applied
write_protein_to_pdb(p_sim1,"p_sim1.pdb")
write_protein_to_pdb(p_sim2,"p_sim2.pdb")
#Plot both proteins
plot_algn_3D(p_sim1,superpos = FALSE,plt_algn = FALSE,plt_connections = TRUE,align = list(I = 1:6,J = 1:6),
             x_cap = "",y_cap = "",z_cap = "",startcol1 = "green",size = 0.5)
plot_algn_3D(p_sim2,superpos = FALSE,plt_algn = FALSE,plt_connections = TRUE,align = list(I = 1:6,J = 1:6),
             x_cap = "",y_cap = "",z_cap = "",col1 = "red",startcol1 = "green",size=0.5)

#Gromov Wasserstein
lt <- get_loss_tensor(p_sim1,p_sim2)

#Run many times to avoid local minima
k <- 1
coupl_gromov <- list()
gw_cost <- vector(length = 15)
set.seed(38376)
while(k <= 15){
  coupl_gromov[[k]] <- partial_gromov_entropic(pr1 =p_sim1,pr2 =p_sim2,
                                               alpha = 4/5,iter_n = 200,
                                               Loss_tensor = lt,tau = 0.02)
  if(anyNA(coupl_gromov[[k]])==FALSE){
    gw_cost[k] <- gromov_wstein_cost(coupl_gromov[[k]], lt)
    k <- k + 1
  } else{
    learn_rate <- learn_rate - 0.02
    print(learn_rate)
  }
}
coupl_gromov <- coupl_gromov[[min(which(gw_cost==min(gw_cost)))]]
align_gromov <- coupling_to_alignNS(coupl_gromov)

get_align_measures(p_sim1,p_sim2,align_gromov)
plot_algn_3D(p_sim1,p_sim2,align = align_gromov,
             superpos = TRUE,plt_algn = TRUE, plt_connections = TRUE,x_cap = "",y_cap = "",z_cap = "",
             size = 0.5,startcol1 = "green",startcol2 = "green")

align_wstein <- align_proteins(p_sim1,p_sim2,return_meas = FALSE,lambda = 0.001,
                               iter = 2,remove_res = FALSE,
                        dist = "partial_wasserstein",print_info = FALSE,alpha = 1)
get_align_measures(p_sim1,p_sim2,align_wstein)

plot_algn_3D(p_sim1,p_sim2,align = align_wstein,
             superpos = TRUE,plt_algn = TRUE, plt_connections = TRUE,x_cap = "",y_cap = "",z_cap = "",
             size = 0.5,startcol1 = "green",startcol2 = "green")



