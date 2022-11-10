#' Transforms a two-mode network into a one-mode network (one-mode projection).
#'
#' This function takes a two-mode network stored as a matrix and transforms it
#' into a one-mode network (i.e., a one-mode projection).
#' The type of one-mode projections for nodes of mode 1 based on two-mode data
#' involving ties between mode 1 and mode 2 nodes are:
#' - "CrossProd": the number of common connections of nodes of the other mode. In the case of a valued network, this give the sum of the product of the tie strength to other nodes of the other modes. For z\[i,j\]=(Sum\[k\](x\[i,k\]*x\[j,k\]))
#' - "CrossProdDMin": the number of common connections of nodes of the other mode divided by the maximum possible given the degree of both nodes. In the case of a valued network, this give the sum of the product of the tie strength to other nodes of the other modes. For z\[i,j\]=(Sum\[k\](x\[i,k\]*x\[j,k\]))/(Sum\[k\](min(x\[i,k\],x\[j,k\])))
#' - "ExactMatch": the relation between two nodes of a mode is calculated as the number of times the two nodes make the same choices and non-choices. In the case of a valued network, the same strength of tie to the other mode. For z\[i,j\]=(Sum\[k\](x\[i,k\]=x\[j,k\])). This can also be seen as the absolute difference between 2.
#' - "Jaccard": the number of common connections to nodes of the other mode, divided by the total number of connections to nodes of mode B by either (or both). For z\[i,j\]=(Sum\[k\](x\[i,k\]*x\[j,k\]))/(Sum\[k\](x\[i,k\]=1 AND/OR xjk=1))
#' - "CrossMin": the strength of a tie between two nodes of a mode is the sum of the minimum of both ties. For z\[i,j\]=(Sum\[k\](min(x\[i,k\],x\[j,k\])))
#' - "Correlation" Pearson's product-moment correlation.
#' - "Covariance" Mean-centered cross products: Sxy/n - SxSy/n^2
#'
#' @param NET2 A single two-mode network stored in an R object as a matrix.
#' @param Measure Character string, indicating which type of transformation to use.
#'
#' @return A one-mode projection
#' @export
#' @importFrom stats cov
#'
#' @examples
#' xTwoModeToOneMode(Davis_SouthernWomen$Attendance)

xTwoModeToOneMode<-function(NET2, Measure="CrossProd")
{
  #Identify the Measure to be included from the full list below:
  MeasureAll<-c("ExactMatch","ExactMatchN0","Jaccard",
                "CrossProd","CrossProdDMin","CrossMin","MaxCrossMin","SumSqDiff",
                "Bonacich","Correlation","Covariance")
  #Give a warning for any unidentified measure
  if (!Measure %in% MeasureAll)
  {
    stop("\nThe following measure could not be identified, and might contain a typo:", "\n", Measure, "\n", "\n")
  }
  #Now go over all options
  if (Measure=="ExactMatch")
  {
    MAT1<-matrix(-Inf,NROW(NET2),NROW(NET2))
    MAT1[lower.tri(MAT1, diag = FALSE)]<-colSums(combn(nrow(NET2),2,FUN=function(x)(NET2[x[1],]==NET2[x[2],])))
    diag(MAT1)<-rowSums(NET2==NET2)
    OUTPUT1<-pmax(MAT1,t(MAT1))
  }
  if (Measure=="ExactMatchN0")
  {
    MAT1<-matrix(-Inf,NROW(NET2),NROW(NET2))
    MAT1[lower.tri(MAT1, diag = FALSE)]<-colSums(combn(nrow(NET2),2,FUN=function(x)((NET2[x[1],]>0)*(NET2[x[2],]>0))))
    diag(MAT1)<-rowSums(NET2==NET2)
    OUTPUT1<-pmax(MAT1,t(MAT1))
  }
  if (Measure=="Jaccard")
  {
    MAT1<-matrix(-Inf,NROW(NET2),NROW(NET2))
    MAT1[lower.tri(MAT1, diag = FALSE)]<-colSums(combn(nrow(NET2),2,FUN=function(x)((NET2[x[1],]>0)*(NET2[x[2],]>0))))
    diag(MAT1)<-rowSums(NET2>0)
    MaxR<-matrix(rowSums(NET2>0),NROW(NET2),NROW(NET2))
    OUTPUT1<-pmax(MAT1,t(MAT1))/(MaxR+t(MaxR)-pmax(MAT1,t(MAT1)))
  }
  if (Measure=="CrossProd")
  {
    OUTPUT1<-crossprod(t(NET2)) #faster than NET2%*%t(NET2)
  }
  if (Measure=="CrossProdDMin")
  {
    SumR<-matrix(rowSums(NET2),NROW(NET2),NROW(NET2))
    OUTPUT1<-(crossprod(t(NET2)))/(pmin(SumR,t(SumR)))
  }
  if (Measure=="CrossMin")
  {
    MAT1<-matrix(-Inf,NROW(NET2),NROW(NET2))
    MAT1[lower.tri(MAT1, diag = FALSE)]<-colSums(combn(nrow(NET2),2,FUN=function(x)pmin(NET2[x[1],],NET2[x[2],])))
    diag(MAT1)<-rowSums(NET2)
    OUTPUT1<-pmax(MAT1,t(MAT1))
  }
  if (Measure=="MaxCrossMin")
  {
    MAT1<-matrix(-Inf,NROW(NET2),NROW(NET2))
    MAT1[lower.tri(MAT1, diag = FALSE)]<-apply(combn(nrow(NET2),2,FUN=function(x)pmin(NET2[x[1],],NET2[x[2],])),2,max)
    diag(MAT1)<-apply(NET2,1,max)
    OUTPUT1<-pmax(MAT1,t(MAT1))
  }
  if (Measure=="SumSqDiff")
  {
    MAT1<-matrix(-Inf,NROW(NET2),NROW(NET2))
    MAT1[lower.tri(MAT1, diag = FALSE)]<-colSums(combn(nrow(NET2),2,FUN=function(x)((NET2[x[1],]-NET2[x[2],])^2)))
    diag(MAT1)<-0
    OUTPUT1<-pmax(MAT1,t(MAT1))
  }
  if (Measure=="Bonacich")
  {
    NET2D<-NET2>0
    X<-(colSums(combn(nrow(NET2),2,FUN=function(x)((NET2D[x[1],]==1)*(NET2D[x[2],]==1)))))*
      (colSums(combn(nrow(NET2),2,FUN=function(x)((NET2D[x[1],]==0)*(NET2D[x[2],]==0)))))
    Y<-(colSums(combn(nrow(NET2),2,FUN=function(x)((NET2D[x[1],]==1)*(NET2D[x[2],]==0)))))*
      (colSums(combn(nrow(NET2),2,FUN=function(x)((NET2D[x[1],]==0)*(NET2D[x[2],]==1)))))
    VAL1<-(X-((X*Y)^.5))/(X-Y)
    VAL1[X==Y]<-.5
    #VAL1[(X==0)*(Y==0)]<-NaN
    MAT1<-matrix(-Inf,NROW(NET2),NROW(NET2))
    MAT1[lower.tri(MAT1, diag = FALSE)]<-VAL1
    diag(MAT1)<-1
    #WHAT TO DO IF SOMEONE HAS NO CONNECTIONS???????
    diag(MAT1)[rowSums(NET2D)==0]<-NaN
    OUTPUT1<-pmax(MAT1,t(MAT1))
  }
  if (Measure=="Correlation")
  {
    OUTPUT1<-cor(t(NET2))
  }
  if (Measure=="Covariance")
  {
    OUTPUT1<-cov(t(NET2))
    #FOR SOME REASON THIS IS TIME 6 IN UCINET
  }
  rownames(OUTPUT1)<-rownames(NET2)
  colnames(OUTPUT1)<-rownames(NET2)
  return(OUTPUT1)
}

#' Transforms a two-mode network into a bipartite network
#'
#' This function takes a two-mode network stored as a matrix and transforms it
#' into a bipartite network where the row and columns consists of both actors and events.
#'
#' @param NET2 A single two-mode network stored in an R object as a matrix.
#' @param First Whether to first put RowNodes or ColumnNodes.
#'
#' @return A one-mode bipartite network
#' @export
#'
#' @examples
#' xTwoModeToBipartite(Davis_SouthernWomen$Attendance)
xTwoModeToBipartite<-function(NET2,First="RowNodes")
{
  if (First=="ColumnNodes")
  {
    NET2<-t(NET2)
  }
  BIPART1<-matrix(0,NROW(NET2)+NCOL(NET2),NROW(NET2)+NCOL(NET2))
  BIPART1[1:NROW(NET2),(NROW(NET2)+1):(NROW(NET2)+NCOL(NET2))]<-NET2
  BIPART1[(NROW(NET2)+1):(NROW(NET2)+NCOL(NET2)),1:NROW(NET2)]<-t(NET2)
  rownames(BIPART1)<-c(rownames(NET2),colnames(NET2))
  colnames(BIPART1)<-rownames(BIPART1)
  return(BIPART1)
}

#' Calculates the cliquemembership for a two-mode network.
#'
#' This function takes a two-mode network stored as a matrix and calculates the bicliques.
#' By default a biclique consists of at least 3 mode A and 3 mode B nodes who are all maximal connected.
#'
#' @param NET2 A single two-mode network stored in an R object as a matrix.
#' @param MinA The minimum number of mode A nodes (row)
#' @param MinB The minimum number of mode B nodes (columns)
#'
#' @return A one-mode bipartite network
#' @export
#'
#' @examples
#' xBiCliques(Davis_SouthernWomen$Attendance)
#'
xBiCliques<-function(NET2, MinA=3, MinB=3)
{
  NR1<-dim(NET2)[1]
  NC1<-dim(NET2)[2]
  BIP<-xTwoModeToBipartite(NET2)
  ONETWO<-((BIP%*%BIP)+BIP)>0
  CL1<-xCliquesMembership(ONETWO, Min=6)
  CL2<-CL1[,colSums(CL1[1:NR1,])>=MinA]
  CL3<-CL2[,colSums(CL2[(NR1+1):(NR1+NC1),])>=MinB]
  rownames(CL3)<-c(rownames(NET2),colnames(NET2))
  CL3
}

#' Performs a dual projection using the Louvain Method.
#'
#' Performs a separate community detection analysis (Louvain) on row and column nodes using one-mode projections.
#'
#' @param NET2 A single two-mode network stored in an R object as a matrix.
#'
#' @return A membership
#' @export
#'
#' @examples
#' xDualLouvainMethod(Davis_SouthernWomen$Attendance)

xDualLouvainMethod<-function(NET2)
{
  AA<-NET2%*%t(NET2)
  EE<-t(NET2)%*%NET2
  ACL<-xLouvainMethod(AA)
  ECL<-xLouvainMethod(EE)
  list(ACL,ECL)
}

#' Performs a dual core-periphery analysis.
#'
#' Performs a separate core-periphery analysis on row and column nodes using one-mode projections.
#'
#' @param NET2 A single two-mode network stored in an R object as a matrix.
#'
#' @return A membership
#' @export
#'
#' @examples
#' xDualCorePeriphery(Davis_SouthernWomen$Attendance)

xDualCorePeriphery<-function(NET2)
{
  AA<-NET2%*%t(NET2)
  EE<-t(NET2)%*%NET2
  ACL<-xCorePeriphery(AA)
  ECL<-xCorePeriphery(EE)
  list(ACL[,5],ECL[,5])
}


#' Performs a dual structural equivalence analysis.
#'
#' Performs a separate structural equivalence analysis on row and column nodes using one-mode projections.
#'
#' @param NET2 A single two-mode network stored in an R object as a matrix.
#'
#' @return A membership
#' @export
#'
#' @examples
#' xDualStructuralEquivalence(Davis_SouthernWomen$Attendance)
xDualStructuralEquivalence<-function(NET2)
{
  AA<-NET2%*%t(NET2)
  EE<-t(NET2)%*%NET2
  ACL<-xStructuralEquivalence(AA,Method="Euclidean",IncludeTransposed = FALSE,Choiceij = "Original")
  ECL<-xStructuralEquivalence(EE,Method="Euclidean",IncludeTransposed = FALSE,Choiceij = "Original")
  list(ACL,ECL)
}

