#' #--------------------------------------------------------------------
#' Provides the clique membership for a one-mode network
#'
#' Using an binary, undirected one-mode network (NET1), this function gives the maximum clique solutions.
#' Clique membership is captured by a matrix, with nodes as rows and each clique as a column.
#' The minimum size of the cliques being considered are defined in (Min).
#' The network should be a binary and undirected one-mode network.
#'
#' @param NET1 A binary, undirected one-mode network stored as a 'matrix' object.
#' @param Min The minimum size for cliques being reported. By default a clique has at least 3 nodes.
#'
#' @return A matrix with membership for each node (row) by clique (column).
#' @importFrom igraph graph_from_adjacency_matrix max_cliques
#' @export
#'
#' @references
#' Chapter 11. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' igraph package - https://igraph.org
#'
#' @seealso [xUCINET::xCliquesCoMembership()], [xUCINET::xCliquesOverlap()], [igraph::max_cliques()]
#'
#' @examples
#' ## Examples of undirected networks (See chapter 11 of Borgatti et al., 2022):
#' xCliquesMembership(Hawthorne_BankWiring$Game)

xCliquesMembership<-function(NET1,Min=3){
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  if(NR1!=NC1){warning(' .   ### Network file NET1 is directed, undirected version used. ')}
  NET1i<-igraph::graph_from_adjacency_matrix(NET1, mode="undirected")
  CLMEMB<-igraph::max_cliques(NET1i,min=Min)
  OUTPUT1<- matrix(0,NR1,length(CLMEMB))
  for(k in 1:length(CLMEMB)){
    OUTPUT1[CLMEMB[[k]],k]<-1
  }
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-paste("CL_",c(1:length(CLMEMB)),sep="")
  OUTPUT1
}

#--------------------------------------------------------------------
#' Provides the clique comembership among nodes for a one-mode network
#'
#' Using an binary, undirected one-mode network (NET1), this function calculates the comembership matrix among nodes.
#' Clique comembership is the number of times two nodes are part of a clique.
#' The network should be a binary and undirected one-mode network.
#'
#' @param NET1 A binary, undirected one-mode network stored as a 'matrix' object.
#' @param Min The minimum size for cliques being considered. By default a clique has at least 3 nodes.
#'
#' @return A matrix with comembership among nodes.
#' @export
#'
#' @references
#' Chapter 11. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' igraph package - https://igraph.org
#'
#' @seealso [xUCINET::xCliquesMembership()], [xUCINET::xCliquesOverlap()], [igraph::max_cliques()]
#'
#' @examples
#' ## Examples of undirected networks (See chapter 11 of Borgatti et al., 2022):
#' xCliquesCoMembership(Hawthorne_BankWiring$Game)

xCliquesCoMembership<-function(NET1,Min=3)
{
  MC1<-xCliquesMembership(NET1,Min=Min)
  OUTPUT1<-MC1%*%t(MC1)
  OUTPUT1
}

#--------------------------------------------------------------------
#' Provides the cliques-overlap based on a one-mode network
#'
#' Using an binary, undirected one-mode network (NET1), this function calculates all the maximum cliques in a network and then defines for the each pair of cliques how many nodes they have in common.
#' The network should be a binary and undirected one-mode network.
#'
#' @param NET1 A binary, undirected one-mode network stored as a 'matrix' object.
#' @param Min The minimum size for cliques being considered. By default a clique has at least 3 nodes.
#'
#' @return A matrix where each row and column are unique cliques found, and the values are the number of nodes two cliques have in common.
#' @export
#'
#' @references
#' Chapter 11. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' igraph package - https://igraph.org
#'
#' @seealso [xUCINET::xCliquesMembership()], [xUCINET::xCliquesMembership()], [igraph::max_cliques()]
#'
#' @examples
#' ## Examples of undirected networks (See chapter 11 of Borgatti et al., 2022):
#' xCliquesOverlap(Hawthorne_BankWiring$Game)

xCliquesOverlap<-function(NET1,Min=3)
{
  MC1<-xCliquesMembership(NET1)
  OUTPUT1<-t(MC1)%*%MC1
  OUTPUT1
}

#--------------------------------------------------------------------
#' Provides the membership of nodes for the community detection using the Girvan-Newman algorithm.
#'
#' Using an undirected one-mode network (NET1), this function calculates the solution(s) for community detection membership using the Girvan-Newman algorithm and graphically shows the result.
#' The output is a matrix providing membership for communities with (k) groups.
#' This number of communities can be identified beforehand in the optional argument (NOC). This can be a single value or a series of valued (as a vector).
#' Values should be at least as large as the number of components in a network.
#' If no values are defined in (NOC), then the solution with the largest value is given (again, where the number of communities is at least as large as the number of components).
#' This function relies on igraph's cluster_edge_betweenness() function.
#'
#' @param NET1 An undirected one-mode network stored as a 'matrix' object.
#' @param NOC A number or vector containing the number of communities (NOC) to be defined. This needs to be at least as large as the number of components in the network.
#' @param vertex.shape An option argument to specify the shapes of the nodes. Needs to be a vector of length equal to the number of nodes.
#' @param vertex.label.cex THe size of the font used for the labels (default when NULL is cex=1).
#' @param layout The layout option used. See igraph. (default when NULL is layout_nicely, which tries to find the optimal layout).
#'
#' @return A vector or matrix with membership for each node (row).
#'
#' @importFrom igraph graph_from_adjacency_matrix cluster_edge_betweenness vcount cut_at
#' @importFrom sna component.dist
#' @importFrom graphics lines
#'
#' @export
#'
#' @references
#' Chapter 11. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' igraph package - https://igraph.org
#'
#' @seealso [xUCINET::xFastGreedy()], [xUCINET::xWalkTrap()], [xUCINET::xLouvainMethod()], [xUCINET::xLabelPropagation()], [igraph::cluster_edge_betweenness()]
#'
#' @examples
#' ## Examples of undirected networks (See chapter 11 of Borgatti et al., 2022):
#' xGirvanNewman(Zachary_KarateClub$Connection, vertex.label.cex=.7)
#' ## Add attribute file related to the club split:
#' ## Define the shapes to be used
#' UseShapes<-c("circle","rectangle")
#' ## Define the vector for nodes -
#' ZKAttrClub<-(Zachary_KarateClub$Attributes$Club>1)+1
#' xGirvanNewman(Zachary_KarateClub$Connection, vertex.shape=UseShapes[ZKAttrClub],
#'               vertex.label.cex=.7)
#' xGirvanNewman(Zachary_KarateClub$Connection, NOC=c(3,4,5),vertex.shape=UseShapes[ZKAttrClub],
#'               vertex.label.cex=.7)
#'
#' ## A second dataset:
#' xGirvanNewman(Kapferer_TailorShop$SociationalT1, vertex.label.cex=.7)

xGirvanNewman<-function(NET1, NOC=NULL, vertex.shape=NULL, vertex.label.cex=NULL, layout=NULL)
{
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  if(sum(is.na(NET1)>0)){stop(' .   ### Network file NET1 cannot contain missing data. ')}
  if(NR1!=NC1){stop(' .   ### Network file NET1 is directed, undirected network needed. PLease symmetrize. ')}
  if(!is.null(NOC)){if(!is.vector(NOC)){stop(' .   ### Argument (NOC) needs to be a single value or vector (or NULL)')}}
  NET1i<-igraph::graph_from_adjacency_matrix(NET1, mode="undirected")
  #RUN THE MODULARITY
  CEB<-igraph::cluster_edge_betweenness(NET1i)
  #GET THE MODULARITY SCORES BY NUMBER OF COMMUNITIES
  MODSCORES<-matrix(c(c(length(CEB$modularity):1),CEB$modularity,c(length(CEB$modularity):1)),length(CEB$modularity),3)
  rownames(MODSCORES)<-paste("Q_",c(dim(MODSCORES)[1]:1),sep="")
  #PLOT MODULARITY SOLUTIONS AND WHICH ARE BEING CONSIDERED
  plot(MODSCORES[,1],MODSCORES[,2], type="b", xlab="Number of Communities", ylab=paste("Modularity ( max =",round(max(CEB$modularity), digits=5),")"))
  #IF NONE SPECIFIED GET MAXIMUM - THE MAX CANNOT BE SMALLER THAN THE NUMBER OF COMPONENTS
  NOCOMPONENTS<-length(sna::component.dist(NET1)$csize)
  MODSCORES[,3]<-(MODSCORES[,3]>=NOCOMPONENTS)

  if(is.null(NOC)){NOC<-max(MODSCORES[MODSCORES[,2]==max(MODSCORES[,2]*MODSCORES[,3],na.rm=T),1])}
  #PLOT MINIMUM NUMBER OF GROUPS GIVEN COMPONENTS
  lines(c(NOCOMPONENTS-.5,NOCOMPONENTS-.5), c(min(CEB$modularity),max(CEB$modularity)), col = "red",lty=2,lwd=3)
  #PLOT WHICH SOLUTIONS ARE REQUESTED
  for (k in 1:length(NOC))
  {
    lines(c(NOC[k],NOC[k]), c(min(CEB$modularity),CEB$modularity[length(CEB$modularity)-NOC[k]+1]), col = "blue",lty=3, lwd=2)
  }
  #NOC - number of communities cannot be smaller than number of components
  if(sum(NOCOMPONENTS>NOC)>0){stop(paste(' .   ### At least one of the requested number of communities is lower than the number of components in the network:',NOCOMPONENTS,'. See the plot with modularity values for more information.'))}
  #GET THE MEMBERSHIP SOLUTION FOR REQUESTED NUMBER OF COMMUNITIES (SINGLE OR MULTIPLE)
  OUTPUT1<-matrix(NA,igraph::vcount(NET1i),length(NOC))
  #ENSURE ALL PLOTS FOR DIFFERENT SOLUTIONS HAVE SAME POSITION
  if(is.null(layout)){LAYOUTSOL1<-layout_nicely(NET1i)}else{if(layout=="layout_with_fr"){LAYOUTSOL1<-layout_with_fr(NET1i)}else{LAYOUTSOL1<-layout_nicely(NET1i)}}
  #ITERATE OVER ALL NUMBER OF CLUSTERS NOC
  for (k in 1:length(NOC))
  {
    SOL1k<-igraph::cut_at(CEB,no=NOC[k])
    CEB$membership<-SOL1k
    OUTPUT1[,k]<-SOL1k
    plot(CEB,NET1i,layout=LAYOUTSOL1, vertex.label.family="", vertex.label.color="black", vertex.label.font=2, vertex.shape=vertex.shape, vertex.label.cex=vertex.label.cex, main=paste("Girvan-Newman with",NOC[k],"communities (Modularity", round(MODSCORES[length(CEB$modularity)-NOC[k],2],digits=5),")"))
  }
  #NOW PREPARE OUTPUT
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-paste("CL_",NOC,sep="")
  OUTPUT1
}


#--------------------------------------------------------------------
#' Provides the membership of nodes for the community detection using the fast greedy algorithm.
#'
#' Using an undirected one-mode network (NET1), this function calculates the solution(s) for community detection membership using the fast greedy algorithm and graphically shows the result.
#' The output is a matrix providing membership for communities with (k) groups.
#' This number of communities can be identified beforehand in the optional argument (NOC). This can be a single value or a series of valued (as a vector).
#' Values should be at least as large as the number of components in a network.
#' If no values are defined in (NOC), then the solution with the largest value is given (again, where the number of communities is at least as large as the number of components).
#' This function relies on igraph's cluster_fast_greedy() function.
#'
#' @param NET1 An undirected one-mode network stored as a 'matrix' object.
#' @param NOC A number or vector containing the number of communities (NOC) to be defined. This needs to be at least as large as the number of components in the network.
#' @param vertex.shape An option argument to specify the shapes of the nodes. Needs to be a vector of length equal to the number of nodes.
#' @param vertex.label.cex The size of the font used for the labels (default when NULL is cex=1).
#' @param layout The layout option used. See igraph. (default when NULL is layout_nicely, which tries to find the optimal layout).
#'
#' @return A vector or matrix with membership for each node (row).
#' @importFrom igraph graph_from_adjacency_matrix cluster_fast_greedy vcount cut_at
#' @importFrom sna component.dist
#' @importFrom graphics lines
#' @export
#'
#' @references
#' Chapter 11. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' igraph package - https://igraph.org
#'
#' @seealso [xUCINET::xGirvanNewman()], [xUCINET::xWalkTrap()], [xUCINET::xLouvainMethod()], [xUCINET::xLabelPropagation()], [igraph::cluster_fast_greedy()]
#'
#' @examples
#' ## Examples of undirected networks (See chapter 11 of Borgatti et al., 2022):
#' xFastGreedy(Zachary_KarateClub$Connection, vertex.label.cex=.7)
#' ## Add attribute file related to the club split:
#' ## Define the shapes to be used
#' UseShapes<-c("circle","rectangle")
#' ## Define the vector for nodes -
#' ZKAttrClub<-(Zachary_KarateClub$Attributes$Club>1)+1
#' xFastGreedy(Zachary_KarateClub$Connection, vertex.shape=UseShapes[ZKAttrClub],
#'             vertex.label.cex=.7)
#' xFastGreedy(Zachary_KarateClub$Connection, NOC=c(3,4,5),vertex.shape=UseShapes[ZKAttrClub],
#'             vertex.label.cex=.7)
#'
#' ## A second dataset:
#' xFastGreedy(Kapferer_TailorShop$SociationalT1, vertex.label.cex=.7, layout="layout_with_fr")

xFastGreedy<-function(NET1, NOC=NULL, vertex.shape=NULL, vertex.label.cex=NULL, layout=NULL)
{
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  if(sum(is.na(NET1)>0)){stop(' .   ### Network file NET1 cannot contain missing data. ')}
  if(NR1!=NC1){stop(' .   ### Network file NET1 is directed, undirected network needed. PLease symmetrize. ')}
  if(!is.null(NOC)){if(!is.vector(NOC)){stop(' .   ### Argument (NOC) needs to be a single value or vector (or NULL)')}}
  NET1i<-igraph::graph_from_adjacency_matrix(NET1, mode="undirected")
  #RUN THE MODULARITY
  CFG<-igraph::cluster_fast_greedy(NET1i)
  #GET THE MODULARITY SCORES BY NUMBER OF COMMUNITIES
  MODSCORES<-matrix(c(c(length(CFG$modularity):1),CFG$modularity,c(length(CFG$modularity):1)),length(CFG$modularity),3)
  rownames(MODSCORES)<-paste("Q_",c(dim(MODSCORES)[1]:1),sep="")
  #PLOT MODULARITY SOLUTIONS AND WHICH ARE BEING CONSIDERED
  plot(MODSCORES[,1],MODSCORES[,2], type="b", xlab="Number of Communities", ylab=paste("Modularity ( max =",round(max(CFG$modularity), digits=5),")"))
  #IF NONE SPECIFIED GET MAXIMUM - THE MAX CANNOT BE SMALLER THAN THE NUMBER OF COMPONENTS
  NOCOMPONENTS<-length(sna::component.dist(NET1)$csize)
  MODSCORES[,3]<-(MODSCORES[,3]>=NOCOMPONENTS)

  if(is.null(NOC)){NOC<-max(MODSCORES[MODSCORES[,2]==max(MODSCORES[,2]*MODSCORES[,3],na.rm=T),1])}
  #PLOT MINIMUM NUMBER OF GROUPS GIVEN COMPONENTS
  lines(c(NOCOMPONENTS-.5,NOCOMPONENTS-.5), c(min(CFG$modularity),max(CFG$modularity)), col = "red",lty=2,lwd=3)
  #PLOT WHICH SOLUTIONS ARE REQUESTED
  for (k in 1:length(NOC))
  {
    lines(c(NOC[k],NOC[k]), c(min(CFG$modularity),CFG$modularity[length(CFG$modularity)-NOC[k]+1]), col = "blue",lty=3, lwd=2)
  }
  #NOC - number of communities cannot be smaller than number of components
  if(sum(NOCOMPONENTS>NOC)>0){stop(paste(' .   ### At least one of the requested number of communities is lower than the number of components in the network:',NOCOMPONENTS,'. See the plot with modularity values for more information.'))}
  #GET THE MEMBERSHIP SOLUTION FOR REQUESTED NUMBER OF COMMUNITIES (SINGLE OR MULTIPLE)
  OUTPUT1<-matrix(NA,igraph::vcount(NET1i),length(NOC))
  #ENSURE ALL PLOTS FOR DIFFERENT SOLUTIONS HAVE SAME POSITION
  if(is.null(layout)){LAYOUTSOL1<-layout_nicely(NET1i)}else{if(layout=="layout_with_fr"){LAYOUTSOL1<-layout_with_fr(NET1i)}else{LAYOUTSOL1<-layout_nicely(NET1i)}}
  #ITERATE OVER ALL NUMBER OF CLUSTERS NOC
  for (k in 1:length(NOC))
  {
    SOL1k<-igraph::cut_at(CFG,no=NOC[k])
    CFG$membership<-SOL1k
    OUTPUT1[,k]<-SOL1k
    plot(CFG,NET1i,layout=LAYOUTSOL1, vertex.label.family="", vertex.label.color="black", vertex.label.font=2, vertex.shape=vertex.shape, vertex.label.cex=vertex.label.cex, main=paste("Fast greedy with",NOC[k],"communities (Modularity", round(MODSCORES[length(CFG$modularity)-NOC[k],2],digits=5),")"))
  }
  #NOW PREPARE OUTPUT
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-paste("CL_",NOC,sep="")
  OUTPUT1
}

#--------------------------------------------------------------------
#' Provides the membership of nodes for the community detection using the walktrap algorithm.
#'
#' Using an undirected one-mode network (NET1), this function calculates the solution(s) for community detection membership using the walktrap algorithm and graphically shows the result.
#' The output is a matrix providing membership for communities with (k) groups.
#' This number of communities can be identified beforehand in the optional argument (NOC). This can be a single value or a series of valued (as a vector).
#' Values should be at least as large as the number of components in a network.
#' If no values are defined in (NOC), then the solution with the largest value is given (again, where the number of communities is at least as large as the number of components).
#' This function relies on igraph's cluster_walktrap() function.
#'
#' @param NET1 An undirected one-mode network stored as a 'matrix' object.
#' @param NOC A number or vector containing the number of communities (NOC) to be defined. This needs to be at least as large as the number of components in the network.
#' @param vertex.shape An option argument to specify the shapes of the nodes. Needs to be a vector of length equal to the number of nodes.
#' @param vertex.label.cex The size of the font used for the labels (default when NULL is cex=1).
#' @param layout The layout option used. See igraph. (default when NULL is layout_nicely, which tries to find the optimal layout).
#'
#' @return A vector or matrix with membership for each node (row).
#'
#' @importFrom igraph graph_from_adjacency_matrix cluster_walktrap vcount cut_at
#' @importFrom sna component.dist
#' @importFrom graphics lines
#'
#' @export
#'
#' @references
#' Chapter 11. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' igraph package - https://igraph.org
#'
#' @seealso [xUCINET::xGirvanNewman()], [xUCINET::xFastGreedy()], [xUCINET::xLouvainMethod()], [xUCINET::xLabelPropagation()], [igraph::cluster_walktrap()]
#'
#' @examples
#' ## Examples of undirected networks (See chapter 11 of Borgatti et al., 2022):
#' xWalkTrap(Zachary_KarateClub$Connection, vertex.label.cex=.7)
#' ## Add attribute file related to the club split:
#' ## Define the shapes to be used
#' UseShapes<-c("circle","rectangle")
#' ## Define the vector for nodes -
#' ZKAttrClub<-(Zachary_KarateClub$Attributes$Club>1)+1
#' xWalkTrap(Zachary_KarateClub$Connection, vertex.shape=UseShapes[ZKAttrClub],
#'           vertex.label.cex=.7)
#' xWalkTrap(Zachary_KarateClub$Connection, NOC=c(3,4,5),vertex.shape=UseShapes[ZKAttrClub],
#'           vertex.label.cex=.7)
#'
#' ## A second dataset:
#' xWalkTrap(Kapferer_TailorShop$SociationalT1, vertex.label.cex=.7, layout="layout_with_fr")

xWalkTrap<-function(NET1, NOC=NULL, vertex.shape=NULL, vertex.label.cex=NULL, layout=NULL)
{
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  if(sum(is.na(NET1)>0)){stop(' .   ### Network file NET1 cannot contain missing data. ')}
  if(NR1!=NC1){stop(' .   ### Network file NET1 is directed, undirected network needed. PLease symmetrize. ')}
  if(!is.null(NOC)){if(!is.vector(NOC)){stop(' .   ### Argument (NOC) needs to be a single value or vector (or NULL)')}}
  NET1i<-igraph::graph_from_adjacency_matrix(NET1, mode="undirected")
  #RUN THE MODULARITY
  CWT<-igraph::cluster_walktrap(NET1i)
  #GET THE MODULARITY SCORES BY NUMBER OF COMMUNITIES
  MODSCORES<-matrix(c(c(length(CWT$modularity):1),CWT$modularity,c(length(CWT$modularity):1)),length(CWT$modularity),3)
  rownames(MODSCORES)<-paste("Q_",c(dim(MODSCORES)[1]:1),sep="")
  #PLOT MODULARITY SOLUTIONS AND WHICH ARE BEING CONSIDERED
  plot(MODSCORES[,1],MODSCORES[,2], type="b", xlab="Number of Communities", ylab=paste("Modularity ( max =",round(max(CWT$modularity), digits=5),")"))
  #IF NONE SPECIFIED GET MAXIMUM - THE MAX CANNOT BE SMALLER THAN THE NUMBER OF COMPONENTS
  NOCOMPONENTS<-length(sna::component.dist(NET1)$csize)
  MODSCORES[,3]<-(MODSCORES[,3]>=NOCOMPONENTS)

  if(is.null(NOC)){NOC<-max(MODSCORES[MODSCORES[,2]==max(MODSCORES[,2]*MODSCORES[,3],na.rm=T),1])}
  #PLOT MINIMUM NUMBER OF GROUPS GIVEN COMPONENTS
  lines(c(NOCOMPONENTS-.5,NOCOMPONENTS-.5), c(min(CWT$modularity),max(CWT$modularity)), col = "red",lty=2,lwd=3)
  #PLOT WHICH SOLUTIONS ARE REQUESTED
  for (k in 1:length(NOC))
  {
    lines(c(NOC[k],NOC[k]), c(min(CWT$modularity),CWT$modularity[length(CWT$modularity)-NOC[k]+1]), col = "blue",lty=3, lwd=2)
  }
  #NOC - number of communities cannot be smaller than number of components
  if(sum(NOCOMPONENTS>NOC)>0){stop(paste(' .   ### At least one of the requested number of communities is lower than the number of components in the network:',NOCOMPONENTS,'. See the plot with modularity values for more information.'))}
  #GET THE MEMBERSHIP SOLUTION FOR REQUESTED NUMBER OF COMMUNITIES (SINGLE OR MULTIPLE)
  OUTPUT1<-matrix(NA,igraph::vcount(NET1i),length(NOC))
  #ENSURE ALL PLOTS FOR DIFFERENT SOLUTIONS HAVE SAME POSITION
  if(is.null(layout)){LAYOUTSOL1<-layout_nicely(NET1i)}else{if(layout=="layout_with_fr"){LAYOUTSOL1<-layout_with_fr(NET1i)}else{LAYOUTSOL1<-layout_nicely(NET1i)}}
  #ITERATE OVER ALL NUMBER OF CLUSTERS NOC
  for (k in 1:length(NOC))
  {
    SOL1k<-igraph::cut_at(CWT,no=NOC[k])
    CWT$membership<-SOL1k
    OUTPUT1[,k]<-SOL1k
    plot(CWT,NET1i,layout=LAYOUTSOL1, vertex.label.family="", vertex.label.color="black", vertex.label.font=2, vertex.shape=vertex.shape, vertex.label.cex=vertex.label.cex, main=paste("Walktrap with",NOC[k],"communities (Modularity", round(MODSCORES[length(CWT$modularity)-NOC[k],2],digits=5),")"))
  }
  #NOW PREPARE OUTPUT
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-paste("CL_",NOC,sep="")
  OUTPUT1
}


#--------------------------------------------------------------------
#' Provides the membership of nodes for the community detection using the Louvain method.
#'
#' Using an undirected one-mode network (NET1), this function calculates the solution for community detection membership using the Louvain method and graphically shows the result.
#' The output is a matrix providing membership for communities with (k) groups.
#' This function relies on igraph's cluster_louvain() function.
#'
#' @param NET1 An undirected one-mode network stored as a 'matrix' object.
#' @param vertex.shape An option argument to specify the shapes of the nodes. Needs to be a vector of length equal to the number of nodes.
#' @param vertex.label.cex The size of the font used for the labels (default when NULL is cex=1).
#' @param layout The layout option used. See igraph. (default when NULL is layout_nicely, which tries to find the optimal layout).
#'
#' @return A vector with membership for each node (row).
#'
#' @importFrom igraph layout_nicely vcount cluster_louvain
#' @importFrom graphics lines
#'
#' @export
#'
#' @references
#' Chapter 11. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' igraph package - https://igraph.org
#'
#' @seealso [xUCINET::xGirvanNewman()], [xUCINET::xFastGreedy()], [xUCINET::xWalkTrap()], [xUCINET::xLabelPropagation()], [igraph::cluster_louvain()]
#'
#' @examples
#' ## Examples of undirected networks (See chapter 11 of Borgatti et al., 2022):
#' xLouvainMethod(Zachary_KarateClub$Connection, vertex.label.cex=.7)
#' ## Add attribute file related to the club split:
#' ## Define the shapes to be used
#' UseShapes<-c("circle","rectangle")
#' ## Define the vector for nodes -
#' ZKAttrClub<-(Zachary_KarateClub$Attributes$Club>1)+1
#' xLouvainMethod(Zachary_KarateClub$Connection, vertex.shape=UseShapes[ZKAttrClub],
#'          vertex.label.cex=.7)
#'
#' ## A second dataset:
#' xLouvainMethod(Kapferer_TailorShop$SociationalT1, vertex.label.cex=.7, layout="layout_with_fr")

xLouvainMethod<-function(NET1, vertex.shape=NULL, vertex.label.cex=NULL, layout=NULL)
{
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  if(sum(is.na(NET1)>0)){stop(' .   ### Network file NET1 cannot contain missing data. ')}
  if(NR1!=NC1){stop(' .   ### Network file NET1 is directed, undirected network needed. PLease symmetrize. ')}
  NET1i<-igraph::graph_from_adjacency_matrix(NET1, mode="undirected")
  #RUN THE MODULARITY
  CLM<-igraph::cluster_louvain(NET1i)
  #GET THE MEMBERSHIP SOLUTION FOR REQUESTED NUMBER OF COMMUNITIES (SINGLE OR MULTIPLE)
  OUTPUT1<-matrix(CLM$membership,igraph::vcount(NET1i),1)
  #USE INPUT PLOT
  if(is.null(layout)){LAYOUTSOL1<-igraph::layout_nicely(NET1i)}else{if(layout=="layout_with_fr"){LAYOUTSOL1<-igraph::layout_with_fr(NET1i)}else{LAYOUTSOL1<-igraph::layout_nicely(NET1i)}}
  #ITERATE OVER ALL NUMBER OF CLUSTERS NOC
  TITLE1<-paste("Louvain method with", max(CLM$membership),"communities ( Modularity", round(as.vector(CLM$modularity),digits=5),")")
  plot(CLM,NET1i,layout=LAYOUTSOL1, vertex.label.family="", vertex.label.color="black", vertex.label.font=2, vertex.shape=vertex.shape, vertex.label.cex=vertex.label.cex, main=TITLE1)
  #NOW PREPARE OUTPUT
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-paste("CL_", max(CLM$membership),sep="")
  OUTPUT1
}

#--------------------------------------------------------------------
#' Provides the membership of nodes for the community detection using label propagation.
#'
#' Using an undirected one-mode network (NET1), this function calculates the solution for community detection membership using label propagation and graphically shows the result.
#' The output is a matrix providing membership for communities with (k) groups.
#' This function relies on igraph's cluster_label_prop() function.
#'
#' @param NET1 An undirected one-mode network stored as a 'matrix' object.
#' @param vertex.shape An option argument to specify the shapes of the nodes. Needs to be a vector of length equal to the number of nodes.
#' @param vertex.label.cex The size of the font used for the labels (default when NULL is cex=1).
#' @param layout The layout option used. See igraph. (default when NULL is layout_nicely, which tries to find the optimal layout).
#'
#' @return A vector with membership for each node (row).
#'
#' @importFrom igraph graph_from_adjacency_matrix cluster_label_prop layout_with_fr layout_nicely
#' @importFrom graphics lines
#'
#' @export
#'
#' @references
#' Chapter 11. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' igraph package - https://igraph.org
#'
#' @seealso [xUCINET::xGirvanNewman()], [xUCINET::xFastGreedy()], [xUCINET::xWalkTrap()], [xUCINET::xLouvainMethod()], [igraph::cluster_label_prop()]
#'
#' @examples
#' ## Examples of undirected networks (See chapter 11 of Borgatti et al., 2022):
#' xLabelPropagation(Zachary_KarateClub$Connection, vertex.label.cex=.7)
#' ## Add attribute file related to the club split:
#' ## Define the shapes to be used
#' UseShapes<-c("circle","rectangle")
#' ## Define the vector for nodes -
#' ZKAttrClub<-(Zachary_KarateClub$Attributes$Club>1)+1
#' xLabelPropagation(Zachary_KarateClub$Connection, vertex.shape=UseShapes[ZKAttrClub],
#'          vertex.label.cex=.7)
#'
#' ## A second dataset:
#' xLabelPropagation(Kapferer_TailorShop$SociationalT1, vertex.label.cex=.7,
#'                   layout="layout_with_fr")

xLabelPropagation<-function(NET1, vertex.shape=NULL, vertex.label.cex=NULL, layout=NULL)
{
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  if(sum(is.na(NET1)>0)){stop(' .   ### Network file NET1 cannot contain missing data. ')}
  if(NR1!=NC1){stop(' .   ### Network file NET1 is directed, undirected network needed. PLease symmetrize. ')}
  NET1i<-igraph::graph_from_adjacency_matrix(NET1, mode="undirected")
  #RUN THE MODULARITY
  CLP<-igraph::cluster_label_prop(NET1i)
  #GET THE MEMBERSHIP SOLUTION FOR REQUESTED NUMBER OF COMMUNITIES (SINGLE OR MULTIPLE)
  OUTPUT1<-matrix(CLP$membership,igraph::vcount(NET1i),1)
  #USE INPUT PLOT
  if(is.null(layout)){LAYOUTSOL1<-igraph::layout_nicely(NET1i)}else{if(layout=="layout_with_fr"){LAYOUTSOL1<-igraph::layout_with_fr(NET1i)}else{LAYOUTSOL1<-igraph::layout_nicely(NET1i)}}
  #ITERATE OVER ALL NUMBER OF CLUSTERS NOC
  TITLE1<-paste("Label Propagation with", max(CLP$membership),"communities ( Modularity", round(as.vector(CLP$modularity),digits=5),")")
  plot(CLP,NET1i,layout=LAYOUTSOL1, vertex.label.family="", vertex.label.color="black", vertex.label.font=2, vertex.shape=vertex.shape, vertex.label.cex=vertex.label.cex, main=TITLE1)
  #NOW PREPARE OUTPUT
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-paste("CL_", max(CLP$membership),sep="")
  OUTPUT1
}
