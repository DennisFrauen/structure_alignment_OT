#Code for several plot layouts

#Basis settings for ggplot
basisplot <- function(data = c(),mapping = c()){
  if(length(data)!= 0 & length(mapping)!= 0){
    basisplt <- ggplot(data,mapping)
  } else{
    basisplt <- ggplot()
  }
  basisplt <- basisplt + 
    theme(axis.text=element_text(size=16), axis.title=element_text(size=17),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          plot.title = element_text(size=17))
  return(basisplt)
}
#3d Structure alignment plot
#Input: pr1,pr2:  protein structures 
#       align:    alignment to be plotted
#       superpos: indicates whether plotted in superimposition
#       size:     size of plotted residues
#       plt_algn: indicates whether green lines are plotted bewteen aligned residues
#       plt_connections: indicates whether protein structure residues are connected in plot
#       width, height: size of the plot
#       x_cap, y_cap, z_cap: Caption of the axes
#       col1, col2: colors of the protein structures
#       startcol1, startcol2: colors of the first residues
#Output: Saves an interactive html file directly to the working space
plot_algn_3D <- function(pr1,pr2 = c(),align = c(),superpos=TRUE,size = 0.2,plt_algn = TRUE,
                         plt_connections=TRUE,width = 1200,height=700,
                         x_cap = "x-axis (Angstr)",y_cap ="y-axis (Angstr)",
                         z_cap = "z-axis (Angstr)",col1 = "blue",col2 = "red",
                         startcol1 = c(),startcol2 = c(),title = "",cex = 1){
  #First apply kabsch algorithm if necissary, to perform superimposition
  if(superpos==TRUE){
    #apply katsch on aligned residues to get rotation matrix and translation vector
    kabsch <- pracma::kabsch(A = t(pr1[align$I,]),B = t(pr2[align$J,]))
    #translate/ rotate first protein accordingly
    Trans_mat <- matrix(kabsch$R, nrow=dim(pr1)[1], ncol=3, byrow=TRUE)
    pr1 <- pr1 %*% t(kabsch$U) +Trans_mat
  }
  #Protein lenghts
  N <- dim(pr1)[1]
  if(length(pr2) > 0){
    M <- dim(pr2)[1]
  } else{
    M <- 0
  }
  #Create dataframes for plotting that contain coordinate and color information
  data_1 <- tibble(x = pr1[,1],y = pr1[,2],z = pr1[,3])
  data_2 <- tibble(x = pr2[,1],y = pr2[,2],z = pr2[,3])
  data <- bind_rows(data_1,data_2)
  #Create color vectors
  if(length(startcol1)>0){
    col_1 <- c(startcol1,ifelse(2:N %in% align$I,rep(col1,N),rep("lightblue",N)))
  } else{
    col_1 <- ifelse(1:N %in% align$I,rep(col1,N),rep("lightblue",N))
  }
  if(M > 0){
    if(length(startcol2)>0){
      col_2 <- c(startcol2,ifelse(2:M %in% align$J,rep(col2,M),rep("orange",M)))  
    } else{
      col_2 <- ifelse(1:M %in% align$J,rep(col2,M),rep("orange",M))
    }
  } else{
    col_2 <- c()
  }
  
  #Plot residues
  par(mar=c(0,0,0,0))
  plot3d( 
    x=data$x, y=data$y, z=data$z, 
    type = 's', 
    col = c(col_1,col_2),
    radius = size,
    xlab=x_cap, ylab=y_cap, zlab=z_cap,cex=cex)
  #Connect residuals
  if(plt_connections == TRUE){
    lines3d(x = data_1$x,y = data_1$y,z = data_1$z,col = col1)
    lines3d(x = data_2$x,y = data_2$y,z = data_2$z,col = col2)
  }
  #Plot alignment
  if(plt_algn==TRUE){
    res1 <- pr1[align$I,]
    res2 <- pr2[align$J,]
    segments3d(x = as.vector(rbind(res1[,1],res2[,1])),
               y = as.vector(rbind(res1[,2],res2[,2])),
               z = as.vector(rbind(res1[,3],res2[,3])),
               col = "green")
  }
  #Write to interactive html file
  htmlwidgets::saveWidget(rglwidget(width = width,height = height), 
                          file.path(getwd(), "align_plot.html"))
  #close rgl device
  rgl.close()
}

#Contour plot of a matrix with possibility to add alignments
#(e.g. Wasserstein distance matrix)
plot_wstein_matrix <- function(wstein_matrix, title = "Wasserstein distance matrix",
                               align1=0,align2=0,plot_align1=FALSE,plot_align2=FALSE){
  #wstein_matrix <- eukl_dist_matrix(pr1,pr2,align1)
  rownames(wstein_matrix) <- paste(1:dim(wstein_matrix)[1])
  colnames(wstein_matrix) <- paste(1:dim(wstein_matrix)[2])
  #Format to tibble for plotting
  wstein_tibble <- wstein_matrix %>% as_tibble() %>% rowid_to_column(var = "Protein1") %>%
    gather(key = "Protein2",value = "value",-1) %>% 
    mutate(p1 = as.integer(Protein1), p2 = as.integer(Protein2))
  #Plotting
  plt <- basisplot(wstein_tibble,aes(p1,p2,fill = value))+
    geom_tile() +
    ggtitle(title) + xlab("Structure 1") + ylab("Structure 2")
  #Plot first alignment in heatmap
  if(plot_align1==TRUE){
    algn1 <- tibble(p1 = align1$I,p2 = align1$J)
    plt <- plt + geom_rect(data=algn1, size=1, fill=NA, colour="red",
                           aes(xmin=p1 - 0.5, xmax=p1 + 0.5, ymin=p2 - 0.5, ymax=p2 + 0.5))
  }
  #Plot second alignment in heatmap
  if(plot_align2==TRUE){
    algn1 <- tibble(p1 = align1$I,p2 = align1$J)
    algn2 <- tibble(p1 = align2$I,p2 = align2$J)
    plt <- plt + geom_rect(data=algn2, size=1, fill=NA, colour="green",
                           aes(xmin=p1 - 0.5, xmax=p1 + 0.5, ymin=p2 - 0.5, ymax=p2 + 0.5))
    #Color overlapping matchings
    #align_inter <- inner_join(algn1,algn2)
    # plt <- plt + geom_rect(data=align_inter, size=1, fill=NA, colour="orange",
    # aes(xmin=p1 - 0.5, xmax=p1 + 0.5, ymin=p2 - 0.5, ymax=p2 + 0.5))
  }
  ggplotly(plt)
}

#Plot the inter-residue distributions of chosen residues as histograms
hist_ird_dist <- function(pr1,pr2=0,h,res1,res2=0,plotp2 = FALSE,dist_cutoff=0,
                          cutoff_by_dist=TRUE,seq_cutoff=0, remove_res1 = c(), 
                          remove_res2 = c()){
  #Extract IRD of first protein
  ird1 <- get_inter_res_dist(pr1,dist_cutoff,cutoff_by_dist,seq_cutoff,remove_res1)[res1]
  #Create tibble for plotting
  plot_tib_list <- list()
  
  for(i in 1:length(res1)){
    plot_tib_list[[i]] <- tibble(Residue = paste(res1[i]),IRD = ird1[[i]])
  }
  
  plot_tib1 <- do.call(rbind,plot_tib_list)
  #Save plot
  plot1 <- ggplot(data = plot_tib1,aes(x=IRD, color=Residue, fill=Residue)) +
    geom_histogram(alpha=0.6, binwidth = h) +
    scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) +
    theme_ipsum() +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    ) +
    xlab("Inter residue distances") +
    facet_wrap(~Residue)
  #Second protein
  if(plotp2 == TRUE){
    ird2 <- get_inter_res_dist(pr2,dist_cutoff,cutoff_by_dist,seq_cutoff,remove_res2)[res2]
    #Create tibble for plotting
    plot_tib_list <- list()
    for(i in 1:length(res2)){
      plot_tib_list[[i]] <- tibble(Residue = paste(res2[i]),IRD = ird2[[i]])
    }
    plot_tib2 <- do.call(rbind,plot_tib_list)
    plot2 <- ggplot(data = plot_tib2,aes(x=IRD, color=Residue, fill=Residue)) +
      geom_histogram(alpha=0.6, binwidth = h) +
      scale_fill_viridis(discrete=TRUE) +
      scale_color_viridis(discrete=TRUE) +
      theme_ipsum() +
      theme(
        legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)
      ) +
      xlab("Inter residue distances") +
      facet_wrap(~Residue)
    gridExtra::grid.arrange(plot1,plot2,ncol=2)
  } else{
    plot1
  }
}

