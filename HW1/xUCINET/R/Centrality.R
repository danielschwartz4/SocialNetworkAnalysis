#' Calculates the outdegree for a network (i.e. for each row).
#'
#' Calculates the (out)degree centrality for a one-mode or two-mode matrix.
#' To calculate the indegree simply transpose the matrix by using t().
#' Values can be binary or valued. For valued networks, the normalized
#' (out)degree ('nDegree') is the average value of the outgoing ties.
#' Any missing data (indicated by values 'NA') are ignored.
#' - Last updated: 1 February 2022.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#' @param Loops Logical statement, whether Loops should be included.
#'
#' @return A matrix with degree and normalized degree.
#' @export
#'
#' @examples
#' ## Directed valued one-mode with missing data and with Loops
#' M1<-matrix(c(6,2,4,6,1,NA,2,NA,3),3,3, byrow=TRUE)
#' xDegreeCentrality(M1, Loops=FALSE)
#'
#' ## To obtain the indegree, use the transposed t()
#' xDegreeCentrality(t(M1), Loops=FALSE)
#'
#' ## Undirected one-mode network
#' xDegreeCentrality(ASNR_Fig09x1)
#'
#' ## Directed one-mode network
#' xDegreeCentrality(ASNR_Fig09x4)
#'
#' ## Undirected two-mode network
#' xDegreeCentrality(Davis_SouthernWomen$Attendance)
#' xDegreeCentrality(t(Davis_SouthernWomen$Attendance))
#'
xDegreeCentrality<-function(NET1, Loops=FALSE)
{
  NN1<-dim(NET1)[1]
  TwoMode1<-dim(NET1)[1]!=dim(NET1)[2]
  if (TwoMode1==TRUE) {cat("TwoMode network\n")}
  if (Loops==FALSE & TwoMode1==FALSE) {diag(NET1)<-NA}
  OUTPUT1<-matrix(NA,NN1,2)
  OUTPUT1[,1]<-rowSums(NET1, na.rm=TRUE)
  OUTPUT1[,2]<-rowMeans(NET1, na.rm=TRUE)
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-c("Degree","nDegree")
  OUTPUT1
}

#' Calculates the eigenvector centrality for a network (i.e. for each row).
#'
#' Calculates the eigenvector centrality for a one-mode network.
#' Tends to be 'problematic' for directed networks.
#' Values can be binary or valued. For valued networks.
#' The normalized version ensures that the maximum value is 1.
#' No missing data are allowed.
#' - Last updated: 1 February 2022.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#' @param Loops Logical statement, whether Loops should be used.
#'
#' @return A matrix with eigenvector centrality and normalized eigenvector centrality (max=1).
#' @importFrom sna evcent
#' @export
#'
#' @examples
#' ##Undirected networks
#' xEigenvectorCentrality(ASNR_Fig09x1)
#' xEigenvectorCentrality(ASNR_Fig09x2)
xEigenvectorCentrality<-function(NET1, Loops=FALSE)
{
  NN1<-dim(NET1)[1]
  if (Loops==FALSE) {diag(NET1)<-0}
  OUTPUT1<-matrix(NA,NN1,2)
  OUTPUT1[,1]<-sna::evcent(NET1)
  OUTPUT1[,2]<-OUTPUT1[,1]/max(OUTPUT1[,1])
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-c("Eigenvector","nEigenvector")
  OUTPUT1
}


#' Calculates Beta centrality (Bonacich's power centrality)
#'
#' Calculates the Beta centrality for a one-mode network, where Beta can be a
#' vector of values or a single value, but must be smaller than 1/highest eigenvalue.
#' Can to be 'problematic' for directed networks.
#' Values can be binary or valued. For valued networks.
#' Different normalizations are available in Option.
#' No missing data are allowed.
#' - Last updated: 1 February 2022.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#' @param Beta A single value or vector contain Beta values which need to be smaller than 1/highest eigenvalue.
#' @param Loops Logical statement, whether Loops should be included.
#' @param Option The normalization chosen, which can be one of the following "None" (no normalization), "MAX1" (set maximum value to 1), "SUM1" (set sum of values to 1), "SSQ1" (set sum of squares equal to 1), "SSQN"  (set sum of squares equal to number of nodes)
#'
#' @return A matrix with Beta centrality values.
#' @export
#'
#' @examples
#' ## Undirected one-mode network
#' xBetaCentrality(ASNR_Fig09x3, c(0,.1,.2,.3,.4,-.1,-.2), Option="SSQN")
#'
#' ## Directed one-mode network
#' xBetaCentrality(ASNR_Fig09x4, Beta=c(0,1,10,100,1000) , Option="SSQN")
#' xBetaCentrality(t(ASNR_Fig09x4), Beta=c(0,1,10,100,1000) , Option="SSQN")
#'
xBetaCentrality<-function(NET1, Beta=0, Loops=FALSE, Option="None")
{
  eigen1<-max(eigen(NET1)$values)
  cat("\n", "Highest eigenvalue =", eigen1)
  cat("\n", "Maximum Beta =", .9999/eigen1)
  cat("\n", "Beta value used =", Beta, "\n", "\n")
  NN1<-dim(NET1)[1]
  I1<-matrix(0, NN1, NN1)
  diag(I1)<-1
  OUTPUT1<-matrix(NA,NN1,length(Beta))
  if (Loops==FALSE) {diag(NET1)<-0}
  kn<-1
  for (Beta.k in Beta)
  {
    OUTPUT0<-rowSums(solve(I1-Beta.k*(NET1))%*%NET1)
    if(Option=="None")
    {OUTPUT1[,kn]<-OUTPUT0}
    if(Option=="MAX1")
    {OUTPUT1[,kn]<-OUTPUT0/max(OUTPUT1)}
    if(Option=="SUM1")
    {OUTPUT1[,kn]<-OUTPUT0/sum(OUTPUT1)}
    if(Option=="SSQ1")
    {OUTPUT1[,kn]<-OUTPUT0*(1/sum(OUTPUT0^2))^.5}
    if(Option=="SSQN")
    {OUTPUT1[,kn]<-OUTPUT0*(NN1/sum(OUTPUT0^2))^.5}
    kn<-kn+1
  }
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-paste("Beta", Beta, sep=".")
  OUTPUT1
}

#' Calculates 3 versions of closeness centrality
#'
#' Calculates the closeness centrality for a one-mode or two-mode network.
#' Three versions of closeness are provided: Freeman's closeness, reciprocal closeness
#' and Valente-Foreman's closeness centrality measure.
#' Each version also is normalized.
#' For directed networks it calculates the out-closeness.
#' To calculate the in-closeness simply transpose the matrix by using t().
#' No missing data are allowed. Loops for the geodesic distance is set to 0.
#' - Last updated: 1 February 2022.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#'
#' @return A matrix with closeness centrality values.
#' @importFrom sna geodist
#' @export
#'
#' @examples
#' ## Undirected one-mode network
#' xClosenessCentrality(Padgett_FlorentineFamilies$Marriage)
#'
#' ## Directed one-mode network (out-closeness followed by in-closeness)
#' xClosenessCentrality(Krackhardt_HighTech$Advice)
#' xClosenessCentrality(t(Krackhardt_HighTech$Advice))
#'
xClosenessCentrality<-function(NET1)
{
  NN1<-dim(NET1)[1]
  #use sna::geodist from sna
  GDIST1<-geodist(NET1)$gdist
  #what is max geodist ignoring inf?
  MAXGDIST<-max(GDIST1[is.finite(GDIST1)])
  #Freeman's version
  GDIST1x<-GDIST1
  GDIST1x[GDIST1==Inf]<-MAXGDIST+1
  OUTPUT1<-rowSums(GDIST1x)
  OUTPUT2<-(NN1-1)/rowSums(GDIST1x)
  #Reciprocal
  RGDIST1<-1/GDIST1
  diag(RGDIST1)<-0
  OUTPUT3<-rowSums(RGDIST1)
  OUTPUT4<-OUTPUT3/(NN1-1)
  #Valente-Foreman
  REVDIST<-(MAXGDIST-GDIST1x+1)
  diag(REVDIST)<-0
  OUTPUT5<-rowSums(REVDIST)
  OUTPUT6<-rowSums(REVDIST)/((NN1-1)*MAXGDIST)
  #Combine all
  OUTPUTALL<-cbind(OUTPUT1,OUTPUT2,OUTPUT3,OUTPUT4, OUTPUT5, OUTPUT6)
  #copy names of nodes
  rownames(OUTPUTALL)<-rownames(NET1)
  #give column names
  colnames(OUTPUTALL)<-c("FreemanCloseness","nFreemanCloseness","ReciprocalCloseness","nReciprocalCloseness","ValenteForemanCloseness","nValenteForemanCloseness")
  #Print
  OUTPUTALL
}

#' Calculate kReach centrality.
#'
#' Calculates the kReach centrality for a one-mode or two-mode network,
#' i.e., the number/proportion of nodes reachable in k steps (if Cumulative=TRUE),
#' or in exactly k steps if Cumulative=FALSE.
#' k is defined with the argument kReach, which can be a single value or a set of values.
#' The Loops for the geodesic distance matrix is by default set to NA.
#' Missing data should be indicated by 'NA'.
#' - Last updated: 1 February 2022.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#' @param kReach A single value or vector contain kReach values.
#' @param Loops Logical statement, whether Loops should be included (default=FALSE).
#' @param Normalize Whether to return the normalized value or raw counts.
#' @param Cumulative Whether to consider nodes reached in exactly k steps (FALSE) or in <=k steps (TRUE).
#'
#' @return A matrix with kReach centrality values.
#' @importFrom sna geodist
#' @export
#'
#' @examples
#' ## Undirected one-mode network
#' xReachCentrality(Padgett_FlorentineFamilies$Marriage, kReach=c(1:6))
#'
#' ## Directed one-mode network (out-reach followed by in-reach)
#' xReachCentrality(Krackhardt_HighTech$Advice, kReach=c(1:6))
#' xReachCentrality(t(Krackhardt_HighTech$Advice), kReach=c(1:6))
#'
xReachCentrality<-function(NET1, kReach=2, Loops=FALSE, Normalize=FALSE, Cumulative=TRUE)
{
  #get unique positive values for k only
  if (Loops==FALSE) {diag(NET1)<-0}
  kReach<-unique(kReach)
  if (sum(kReach%%1!=0)>0)
  {stop(' ERROR: Some values for k are not positive integer values')}
  if (sum(kReach<0)>0)
  {stop(' ERROR: Some values for k are not positive integer values')}
  NN1<-dim(NET1)[1]
  OUTPUT1<-matrix(NA,NN1,length(kReach))
  GDIST<-geodist(NET1)$gdist
  if (Loops==FALSE) {diag(GDIST)<-NA}
  kn<-1
  for (kReach.k in kReach)
  {
    if(Cumulative==F)
    {
      if(Normalize==T)
      {
        OUTPUT1[,kn]<-(rowMeans(GDIST==kReach.k, na.rm = TRUE))
      }
      if(Normalize==F)
      {
        OUTPUT1[,kn]<-(rowSums(GDIST==kReach.k, na.rm = TRUE))
      }
    }
    if(Cumulative==T)
    {
      if(Normalize==T)
      {
        OUTPUT1[,kn]<-(rowMeans(GDIST<=kReach.k, na.rm = TRUE))
      }
      if(Normalize==F)
      {
        OUTPUT1[,kn]<-(rowSums(GDIST<=kReach.k, na.rm = TRUE))
      }
    }
    kn<-kn+1
  }
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-paste("kReach", kReach, sep=".")
  OUTPUT1
}

#' Calculate Beta Reach centrality.
#'
#' Calculates the Beta Reach centrality for a one-mode or two-mode network,
#' i.e., the number/proportion of nodes reachable in k steps weighted by Beta.
#' Beta can be a single value or set of values.(default=.5)
#' The Loops for the geodesic distance matrix is by default set to NA.
#' Missing data should be indicated by 'NA'.
#' - Last updated: 1 February 2022.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#' @param Beta A single value or vector contain Beta values.
#' @param Loops Logical statement, whether Loops should be included (default=FALSE).
#' @param Normalize Whether to return the normalized value (divide by theoretical maximum value) or raw counts.
#'
#' @return A matrix with Beta reach centrality values.
#' @importFrom sna geodist
#' @export
#'
#' @examples
#' ## Undirected one-mode network
#' xBetaReachCentrality(Padgett_FlorentineFamilies$Marriage, Beta=c(.5,.3))
#
#' ## Directed one-mode network (out-beta reach followed by in-beta reach)
#' xBetaReachCentrality((Krackhardt_HighTech$Advice), Beta=c(.4), Normalize=TRUE)
#' xBetaReachCentrality(t(Krackhardt_HighTech$Advice), Beta=c(.4), Normalize=TRUE)
#'
xBetaReachCentrality<-function(NET1, Beta=.5, Loops=FALSE, Normalize=FALSE)
{
  if (Loops==FALSE) {diag(NET1)<-0}
  #use sna::geodist from sna
  NN1<-dim(NET1)[1]
  GDIST1<-geodist(NET1)$gdist
  if (Loops==FALSE) {diag(GDIST1)<-NA}
  OUTPUT1<-matrix(NA,NN1,length(Beta))
  if (Normalize==F)
  {
    kn<-1
    for (Beta.k in Beta)
    {
      OUTPUT1[,kn]<-rowSums(Beta.k^(GDIST1-1), na.rm=T)
      kn<-kn+1
    }
  }
  if (Normalize==T)
  {
    kn<-1
    for (Beta.k in Beta)
    {
      OUTPUT1[,kn]<-rowMeans(Beta.k^(GDIST1-1), na.rm=T)
      kn<-kn+1
    }
  }
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-paste("Beta", Beta, sep=".")
  OUTPUT1
}

#' Calculates betweenness centrality.
#'
#' Calculates betweenness centrality for a one-mode or two-mode network,
#' i.e., the number (proportion) of times a node 'i' is in between others.
#' For directed networks is consider the "in-out" connections around node i,
#' i.e., j->i->k.
#' Note: Undirected counts are counted in both directions j->i->k and k->i->j.
#' Provide both non-normalized and normalized measure (divided by the theoretical maximum).
#' No missing data are allowed.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#'
#' @return A matrix with raw and normalized betweenness centrality values.
#' @importFrom sna geodist
#' @export
#'
#' @examples
#' ## Undirected one-mode network
#' xBetweennessCentrality(Padgett_FlorentineFamilies$Marriage)
#'
#' ## Directed one-mode network
#' xBetweennessCentrality(Krackhardt_HighTech$Advice)
#'
xBetweennessCentrality<-function(NET1)
{
  diag(NET1)<-0
  NN1<-dim(NET1)[1]
  #use sna::geodist from sna
  GDIST1<-geodist(NET1)
  GDIST1
  OUTPUT1<-matrix(NA,NN1,2)
  for (node.k in c(1:NN1))
  {
    MAT2<-NET1
    MAT2[node.k,]<-0
    MAT2[,node.k]<-0
    GDIST2<-geodist(MAT2)
    CHANGE_GDIST<-(GDIST2$gdist!=GDIST1$gdist)
    CHANGE_GDIST[node.k,]<-0
    CHANGE_GDIST[,node.k]<-0
    CHANGE_COUNTS<-(GDIST2$gdist==GDIST1$gdist)*(GDIST1$counts-GDIST2$counts)/(GDIST1$counts+(GDIST1$counts==0))
    CHANGE_COUNTS[node.k,]<-0
    CHANGE_COUNTS[,node.k]<-0
    OUTPUT1[node.k,1]<-sum(CHANGE_GDIST+CHANGE_COUNTS,na.rm=T)
  }
  OUTPUT1[,2]<-OUTPUT1[,1]/((NN1-1)*(NN1-2))
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-c("Betweenness","nBetweenness")
  OUTPUT1
}

#' Calculates reach betweenness centrality.
#'
#' Calculates betweenness centrality for a one-mode or two-mode network,
#' i.e., the number (proportion) of times a node 'i' is in between others
#' for paths of length k or less (defined by 'kReach').
#' kReach has to be a single value.
#' For directed networks is consider the "in-out" connections around node i,
#' i.e., j->i->k, which are of length k or less.
#' Note: Undirected counts are counted in both directions j->i->k and k->i->j.
#' No missing data are allowed.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#' @param kReach A single value indicating the maximum distances being considered.
#'
#' @return A matrix with reach betweenness centrality values for a specific k (kReach).
#' @importFrom sna geodist
#' @export
#'
#' @examples
#' ## Undirected one-mode network
#' xReachBetweennessCentrality(Padgett_FlorentineFamilies$Marriage)
#'
#' ## Directed one-mode network
#' xReachBetweennessCentrality(Krackhardt_HighTech$Advice)
#'
xReachBetweennessCentrality<-function(NET1, kReach=2)
{
  diag(NET1)<-0
  NN1<-dim(NET1)[1]
  #use sna::geodist from sna
  GDIST1<-sna::geodist(NET1)
  OUTPUT1<-matrix(NA,NN1,2)
  for (node.k in c(1:NN1))
  {
    MAT2<-NET1
    MAT2[node.k,]<-0
    MAT2[,node.k]<-0
    GDIST2<-sna::geodist(MAT2)
    CHANGE_GDIST<-(GDIST2$gdist!=GDIST1$gdist)*(GDIST1$gdist<=kReach)
    CHANGE_GDIST[node.k,]<-0
    CHANGE_GDIST[,node.k]<-0
    CHANGE_COUNTS<-(GDIST2$gdist==GDIST1$gdist)*(GDIST1$counts-GDIST2$counts)/(GDIST1$counts+(GDIST1$counts==0))*(GDIST1$gdist<=kReach)
    CHANGE_COUNTS[node.k,]<-0
    CHANGE_COUNTS[,node.k]<-0
    OUTPUT1[node.k,1]<-sum(CHANGE_GDIST+CHANGE_COUNTS,na.rm=T)
  }
  OUTPUT1[,2]<-OUTPUT1[,1]/((NN1-1)*(NN1-2))
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-c("ReachBetweenness","nReachBetweenness")
  OUTPUT1
}

#' Calculate Beta Reach betweenness centrality.
#'
#' Calculates the Beta Reach between centrality for a one-mode or two-mode network,
#' i.e., the number (proportion) of times a node is in between others, where this is weights by its length.
#' Beta needs to be a single value (default=.5)
#' Missing data should be indicated by 'NA'.
#' - Last updated: 1 February 2022.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#' @param Beta A single value for Beta (max = 1).
#'
#' @return A matrix with Beta reach betweenness centrality values.
#' @importFrom sna geodist
#' @export
#'
#' @examples
#' ## Undirected one-mode network
#' xBetaReachBetweennessCentrality(Padgett_FlorentineFamilies$Marriage, Beta=.5)
#
#' ## Directed one-mode network
#' xBetaReachBetweennessCentrality(Krackhardt_HighTech$Advice, Beta=.3)
#'
xBetaReachBetweennessCentrality<-function(NET1, Beta=.5)
{
  diag(NET1)<-0
  NN1<-dim(NET1)[1]
  #use sna::geodist from sna
  GDIST1<-geodist(NET1)
  OUTPUT1<-matrix(NA,NN1,2)
  for (node.k in c(1:NN1))
  {
    MAT2<-NET1
    MAT2[node.k,]<-0
    MAT2[,node.k]<-0
    GDIST2<-geodist(MAT2)
    CHANGE_GDIST<-(GDIST2$gdist!=GDIST1$gdist)*(Beta^(GDIST1$gdist-1))
    CHANGE_GDIST[node.k,]<-0
    CHANGE_GDIST[,node.k]<-0
    CHANGE_COUNTS<-(GDIST2$gdist==GDIST1$gdist)*(GDIST1$counts-GDIST2$counts)/(GDIST1$counts+(GDIST1$counts==0))*(Beta^(GDIST1$gdist-1))
    CHANGE_COUNTS[node.k,]<-0
    CHANGE_COUNTS[,node.k]<-0
    OUTPUT1[node.k,1]<-sum(CHANGE_GDIST+CHANGE_COUNTS,na.rm=T)
  }
  OUTPUT1[,2]<-OUTPUT1[,1]/((NN1-1)*(NN1-2))
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-c("BetaReachBetw","nBetaReachBetw")
  OUTPUT1
}

#' Calculate induced centrality.
#'
#' Induced centrality calculates the change in specific network properties (graph invariants) when a node is removed.
#' A graph invariant can be used to induce a centrality score by measuring the contribution each node makes to the network property (e.g. the number of components).
#' The induced centrality of a node is the value for the graph invariant for the observed network minus the value for the graph invariant for the observed network but with the node deleted.
#' This routine calculates the induced centrality for the following network properties:
#' - "kReach1": Number of nodes reachable within geodesic distance 1 (Degree)
#' - "kReach2": Number of nodes reachable within geodesic distance 2
#' - "kReach3": Number of nodes reachable within geodesic distance 3
#' - "kReach4": Number of nodes reachable within geodesic distance 4
#' - "kReachAll": Number of nodes reachable within a finite number of nodes
#' - "RecipCloseness": Sum of reciprocal geodesic distance to all other nodes
#' - "Betweenness": Number of times a node is between other nodes
#' No missing data are allowed.
#' - Last updated: 1 February 2022.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#'
#' @return A matrix with a series of induced centrality measures, followed by a split between endogeneous and exogeneous.
#' @importFrom sna geodist
#' @export
#'
#' @examples
#' ## Undirected one-mode network
#' xInducedCentrality(Padgett_FlorentineFamilies$Marriage)
#'
#' ## Directed one-mode network
#' xInducedCentrality(Krackhardt_HighTech$Advice)
#' xInducedCentrality(t(Krackhardt_HighTech$Advice))

xInducedCentrality<-function(NET1)
{
  diag(NET1)<-0
  NET1N<-NET1
  diag(NET1N)<-NA
  NN1<-dim(NET1)[1]
  OUTPUT0<-matrix(NA,NN1,7)
  OUTPUTR<-matrix(NA,NN1,7)
  OUTPUT1<-matrix(NA,NN1,7)
  MATQ<-matrix(1,NN1,NN1)
  diag(MATQ)<-0
  GDIST<-geodist(NET1)$gdist
  diag(GDIST)<-NA
  OUTPUT0[,1]<-sum(GDIST<2,na.rm=T)
  OUTPUT0[,2]<-sum(GDIST<3,na.rm=T)
  OUTPUT0[,3]<-sum(GDIST<4,na.rm=T)
  OUTPUT0[,4]<-sum(GDIST<5,na.rm=T)
  OUTPUT0[,5]<-sum(GDIST<NN1,na.rm=T)
  OUTPUT0[,6]<-sum(1/GDIST,na.rm=T)
  OUTPUT0[,7]<-sum(xBetweennessCentrality(NET1)[,1], na.rm=T)
  for (node.k in c(1:NN1))
  {
    MAT2<-NET1
    MAT2[node.k,]<-0
    MAT2[,node.k]<-0
    GDIST2<-geodist(MAT2)$gdist
    diag(GDIST2)<-NA
    OUTPUTR[node.k,1]<-sum(GDIST2<2,na.rm=T)
    OUTPUTR[node.k,2]<-sum(GDIST2<3,na.rm=T)
    OUTPUTR[node.k,3]<-sum(GDIST2<4,na.rm=T)
    OUTPUTR[node.k,4]<-sum(GDIST2<5,na.rm=T)
    OUTPUTR[node.k,5]<-sum(GDIST2<NN1,na.rm=T)
    OUTPUTR[node.k,6]<-sum(1/GDIST2,na.rm=T)
    OUTPUTR[node.k,7]<-sum(xBetweennessCentrality(MAT2)[,1], na.rm=T)
  }
  OUTPUT1[,1]<-OUTPUT0[,1]-OUTPUTR[,1]
  OUTPUT1[,2]<-OUTPUT0[,2]-OUTPUTR[,2]
  OUTPUT1[,3]<-OUTPUT0[,3]-OUTPUTR[,3]
  OUTPUT1[,4]<-OUTPUT0[,4]-OUTPUTR[,4]
  OUTPUT1[,5]<-OUTPUT0[,5]-OUTPUTR[,5]
  OUTPUT1[,6]<-OUTPUT0[,6]-OUTPUTR[,6]
  OUTPUT1[,7]<-OUTPUT0[,7]-OUTPUTR[,7]
  OUTPUT2<-xReachCentrality(NET1,kReach=c(1,2,3,4,NN1))
  OUTPUT3<-xClosenessCentrality(NET1)[,3]
  OUTPUT4<-xBetweennessCentrality(NET1)[,1]
  OUTPUT1<-cbind(OUTPUT1,OUTPUT2,OUTPUT3,OUTPUT4)
  OUTPUT1<-cbind(OUTPUT1,OUTPUT1[,1:7]-OUTPUT1[,8:14])

  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-c("kReach1","kReach2","kReach3","kReach4","kReachAll","RecipCloseness#","Betweenness",
                       "kReach1R","kReach2R","kReach3R","kReach4R","kReachAllR","RecipClosenessR","BetweennessR",
                       "kReach1D","kReach2D","kReach3D","kReach4D","kReachAllD","RecipClosenessD#","BetweennessD")
  OUTPUT1
}

#' Calculates the outdegree for negative relations (i.e. for each row).
#'
#' Calculates the (out)degree centrality for a one-mode containing a negative relation.
#' To calculate the indegree simply transpose the matrix by using t().
#' This function calculates degree centrality for negative relations as explained in Everett and Borgatti (2014).
#' Only for binary data.
#' Any missing data (indicated by values 'NA') are ignored.
#' - Last updated: 1 February 2022.
#'
#' @param NET1N A one-mode or two-mode negative network stored as a 'matrix' object.
#'
#' @return A matrix with outdegree centrality measures.
#' @export
#'
#' @examples
#' ## Undirected one-mode network # TABLE 1
#' xNegativeDegreeCentrality(Read_NewGuinea$Opposition)
#'
#' ## Directed one-mode network # TABLE 3
#' NEG1<-((Sampson_Monastery$Disesteem+Sampson_Monastery$Dislike+Sampson_Monastery$Blame
#'         +Sampson_Monastery$NegativeInfluence)>0)*1
#' xNegativeDegreeCentrality(NEG1)
#' xNegativeDegreeCentrality(t(NEG1))
#'
#' @references
#' Everett, M.G. Borgatti, S.P. (2014) Networks containing negative ties. Social Networks 38, 111-120.

xNegativeDegreeCentrality<-function(NET1N)
{
  diag(NET1N)<-NA
  OUTPUT1<-rowSums(!is.na(NET1N))
  OUTPUT2<-rowSums(-NET1N,na.rm=T)
  OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
  OUTPUT2<-1-rowSums(NET1N,na.rm=T)/OUTPUT1[,1]
  OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
  rownames(OUTPUT1)<-rownames(NET1N)
  colnames(OUTPUT1)<-c("Valid","DegreeNegative","nDegreeNegative")
  OUTPUT1
}

#' Calculates the weighted centrality for negative relations (i.e. for each row).
#'
#' This function calculates a weighted centrality for negative relations as explained in Everett and Borgatti (2014).
#' No missing data allowed.
#' - Last updated: 1 February 2022.
#'
#' @param NET1N A one-mode negative network stored as a 'matrix' object.
#'
#' @return A matrix with weighted centrality measures.
#' @export
#'
#' @examples
#' ## Undirected one-mode network # TABLE 1
#' xNegativeWeightedCentrality(Read_NewGuinea$Opposition)
#'
#' ## Directed one-mode network # TABLE 3
#' NEG1<-((Sampson_Monastery$Disesteem+Sampson_Monastery$Dislike+Sampson_Monastery$Blame
#'         +Sampson_Monastery$NegativeInfluence)>0)*1
#' xNegativeWeightedCentrality(NEG1)
#' xNegativeWeightedCentrality(t(NEG1))
#'
#' @references
#' Everett, M.G. Borgatti, S.P. (2014) Networks containing negative ties. Social Networks 38, 111-120.

xNegativeWeightedCentrality<-function(NET1N)
{
  NN1<-dim(NET1N)[1]
  diag(NET1N)<-0
  IDENTM<-diag(NN1)
  VECT1<-matrix(1,NN1,1)
  OUTPUT1<-solve(IDENTM-((NET1N%*%t(NET1N))/((NN1-1)^2)))%*%(IDENTM-( NET1N/(NN1-1) ))%*%VECT1
  rownames(OUTPUT1)<-rownames(NET1N)
  colnames(OUTPUT1)<-c("WeightedNegative")
  OUTPUT1
}

#' PN centrality for a negative and positive relation.
#'
#' An actor's PN centrality is one plus the weighted sum of the PN centrality of the actors
#' they are positively connected to minus a weighted sum of the PN centrality scores of
#' the actors they are negatively connected to. The default weights are used to normalize the scores.
#' The weights are expressed as a beta value (this is because if there are no negative ties the measure
#' is the same as the Hubbell measure) and a weight w for the negative ties
#' No missing data allowed.
#' - Last updated: 1 February 2022.
#'
#' @param NET1P A one-mode positive network stored as a 'matrix' object.
#' @param NET1N A one-mode negative network stored as a 'matrix' object.
#'
#' @return A matrix with PN centrality measures.
#' @export
#'
#' @examples
#' ## Undirected one-mode network ## TABLE 4
#' NEG1<-((Sampson_Monastery$Disesteem+Sampson_Monastery$Dislike+Sampson_Monastery$Blame
#'         +Sampson_Monastery$NegativeInfluence)>0)*1
#' POS1<-((Sampson_Monastery$Esteem+Sampson_Monastery$Praise+Sampson_Monastery$PositiveInfluence
#'         +Sampson_Monastery$LikeT3)>0)*1
#' POS2<-(POS1+t(POS1))>0
#' NEG2<-(NEG1+t(NEG1))>0
#' xPNCentrality(POS2,NEG2)
#'
#' ## Directed one-mode network ## not the same as TABLE 7
#' xPNCentrality(POS1,NEG1)
#' xPNCentrality(t(POS1),t(NEG1))
#'
#' @references
#' Everett, M.G. Borgatti, S.P. (2014) Networks containing negative ties. Social Networks 38, 111-120.

xPNCentrality<-function(NET1P,NET1N)
{
  NN1<-dim(NET1N)[1]
  PN1<-NET1P-NET1N*2
  diag(PN1)<-0
  IDENTM<-diag(NN1)
  VECT1<-matrix(1,NN1,1)
  OUTPUT1<-solve(IDENTM-((PN1%*%t(PN1))/(4*(NN1-1)^2)))%*%(IDENTM+( PN1/(2*(NN1-1)) ))%*%VECT1
  rownames(OUTPUT1)<-rownames(NET1N)
  colnames(OUTPUT1)<-c("PNcentrality")
  OUTPUT1
}

# ---------------------------- END -------------------------------------
